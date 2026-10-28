//! Uniqueness flags: which of an index's values only genomes of their own species share.
//!
//! The index holds a species' representative genome(s); the pass reads a larger set of genomes, each
//! labelled with its species, and looks up every k-mer the index would have selected. A value of the
//! index loses
//!  - [`UNIQUE`] when a genome of another species has a k-mer that matches it: on the core-mer alone
//!    for a key without a header (a short value), on core-mer and flanks for a key with one (a long
//!    value);
//!  - [`SUPER_UNIQUE`] (keys with a header only) when a genome of another species has a k-mer of the
//!    same core-mer within one flank mismatch of it.
//!
//! The three classes protal reports follow from the flags and `has_header`:
//! short unique = `UNIQUE` without a header, long unique = `UNIQUE` with a header, long super unique
//! = `UNIQUE | SUPER_UNIQUE` (with a header). The genome set should contain the representatives, so
//! that representatives of two species are checked against each other too.
//!
//! The flags live in a value cell's two top bits ([`crate::values::FLAG_SHIFT`]) and in a packed
//! index's `flag_bits`. An index built without the pass has none, and reads exactly as before.

use std::cmp::min;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use kmerrs::consecutive::kmer::KmerIter;
use kmerrs::kmer_utils::mask_coremer;
use kmerrs::minimizer::context_free::Minimizer;
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;

use crate::flexmap::Flexmap;
use crate::packed::{flank_distance, table_keys};
use crate::values::FLAG_SHIFT;
use crate::VD;

/// No genome of another species has this value's k-mer.
pub const UNIQUE: u8 = 1;
/// (Keys with a header) no genome of another species has a k-mer within one flank mismatch of it.
pub const SUPER_UNIQUE: u8 = 2;

/// Species names and their ids.
#[derive(Default, Debug, Clone)]
pub struct SpeciesTable {
    pub names: Vec<String>,
    ids: HashMap<String, u32>,
}

impl SpeciesTable {
    /// The id of species `name`, a new one if it is new.
    pub fn id(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.ids.get(name) {
            return id;
        }
        let id = self.names.len() as u32;
        self.names.push(name.to_string());
        self.ids.insert(name.to_string(), id);
        id
    }
}

/// Reads a two-column TSV (`name<TAB>species`; `#` lines and blank lines skipped) as pairs.
pub fn read_pairs(path: impl AsRef<Path>) -> std::io::Result<Vec<(String, String)>> {
    let mut pairs = Vec::new();
    for line in BufReader::new(File::open(path.as_ref())?).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut cols = line.split('\t');
        match (cols.next(), cols.next()) {
            (Some(a), Some(b)) => pairs.push((a.to_string(), b.to_string())),
            _ => return Err(std::io::Error::other(format!("{}: not two columns: {line:?}", path.as_ref().display()))),
        }
    }
    Ok(pairs)
}

/// The species of every reference id of an index (`id2reference` as the build returns it, id 0 the
/// dummy), from `sequence -> species`. Fails on a reference without a species.
pub fn reference_species(
    id2reference: &[String],
    sequence_species: &HashMap<String, String>,
    species: &mut SpeciesTable,
) -> std::io::Result<Vec<u32>> {
    let mut out = vec![u32::MAX; id2reference.len()];
    for (id, name) in id2reference.iter().enumerate().skip(1) {
        let s = sequence_species
            .get(name)
            .ok_or_else(|| std::io::Error::other(format!("reference {name:?} has no species")))?;
        out[id] = species.id(s);
    }
    Ok(out)
}

/// What the pass found, counted over the index's values.
#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub struct UniquenessStats {
    pub genomes: usize,
    pub kmers_looked_up: u64,
    pub values: u64,
    pub short_values: u64,
    pub short_unique: u64,
    pub long_values: u64,
    pub long_unique: u64,
    pub long_super_unique: u64,
}

/// Reads a FASTA file (gzip or plain) as its records' (name, upper-case sequence).
fn read_fasta(path: &Path) -> std::io::Result<Vec<(String, Vec<u8>)>> {
    let mut raw = Vec::new();
    File::open(path)?.read_to_end(&mut raw)?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let mut plain = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..]).read_to_end(&mut plain)?;
        raw = plain;
    }
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in raw.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if let Some(head) = line.strip_prefix(b">") {
            let name = String::from_utf8_lossy(head).split(' ').next().unwrap_or("").to_string();
            records.push((name, Vec::new()));
        } else if let Some((_, seq)) = records.last_mut() {
            seq.extend(line.iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    Ok(records)
}

/// Sets the uniqueness flags of `flexmap`'s values: every value starts unique (and, with a header,
/// super unique), then every genome of `genomes` (FASTA file, species id) clears the flags its k-mers
/// contradict. `ref_species[id]` is the species of reference id `id`. K, C, F, S, L must be the ones
/// the index was built with (as `build::default_build`). Runs on `threads` threads, a genome file
/// each at a time.
pub fn mark_uniqueness<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CPB: u64,
    const T: usize,
>(
    flexmap: &mut Flexmap<C, F, CPB, T>,
    ref_species: &[u32],
    genomes: &[(PathBuf, u32)],
    threads: usize,
) -> std::io::Result<UniquenessStats> {
    // Start: everything unique.
    {
        let cells = &mut flexmap.values.data;
        table_keys::<C, CPB>(&flexmap.keys.data)(&mut |_key, start, slots| {
            let (h, flags) = if slots > T as u64 { ((slots + 2) / 3, UNIQUE | SUPER_UNIQUE) } else { (0, UNIQUE) };
            for c in &mut cells[(start + h) as usize..(start + slots) as usize] {
                c.0 |= (flags as u64) << FLAG_SHIFT;
            }
        });
    }

    let flexmap: &Flexmap<C, F, CPB, T> = flexmap;
    // The cells as atomics: flags are only ever cleared, so threads may race on them.
    let cells: &[AtomicU64] = unsafe {
        std::slice::from_raw_parts(flexmap.values.data.as_ptr() as *const AtomicU64, flexmap.values.data.len())
    };
    let clear = |cell: &AtomicU64, flags: u8| {
        let bits = (flags as u64) << FLAG_SHIFT;
        if cell.load(Ordering::Relaxed) & bits != 0 {
            cell.fetch_and(!bits, Ordering::Relaxed);
        }
    };

    let next = AtomicUsize::new(0);
    let looked_up = AtomicU64::new(0);
    let threads = threads.max(1).min(genomes.len().max(1));
    std::thread::scope(|scope| -> std::io::Result<()> {
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| -> std::io::Result<()> {
                    let mut cs = ClosedSyncmer::<K, S, L>::new();
                    let mut lookups = 0u64;
                    loop {
                        let g = next.fetch_add(1, Ordering::Relaxed);
                        let Some((path, species)) = genomes.get(g) else { break };
                        let records = read_fasta(path)
                            .map_err(|e| std::io::Error::other(format!("{}: {e}", path.display())))?;
                        for (_, seq) in &records {
                            // The selection and orientation of `build::default_build_map`.
                            for (_, kmer_fwd, kmer_rev) in KmerIter::<K, true>::new(seq) {
                                let cmer_fwd = kmer_fwd.middle::<C>();
                                let cmer_rev = kmer_rev.middle::<C>();
                                let kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
                                let cmer = min(cmer_fwd, cmer_rev);
                                if !cs.is_minimizer(cmer.0) {
                                    continue;
                                }
                                let Some((start, end)) = flexmap.keys.vrange(mask_coremer::<C>(cmer.0)) else { continue };
                                lookups += 1;
                                let slots = end - start;
                                let other = |c: &AtomicU64| {
                                    let (r, _) = VD::get(c.load(Ordering::Relaxed));
                                    ref_species[r as usize] != *species
                                };
                                if slots > T {
                                    let h = (slots + 2) / 3;
                                    let n = slots - h;
                                    let flank = kmer.flanks::<F>().0 as u32;
                                    for i in 0..n {
                                        let header = cells[start + i / 2].load(Ordering::Relaxed);
                                        let d = flank_distance((header >> (32 * (i % 2))) as u32, flank);
                                        if d > 1 {
                                            continue;
                                        }
                                        let c = &cells[start + h + i];
                                        if other(c) {
                                            clear(c, if d == 0 { UNIQUE | SUPER_UNIQUE } else { SUPER_UNIQUE });
                                        }
                                    }
                                } else {
                                    for c in &cells[start..end] {
                                        if other(c) {
                                            clear(c, UNIQUE);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    looked_up.fetch_add(lookups, Ordering::Relaxed);
                    Ok(())
                })
            })
            .collect();
        for w in workers {
            w.join().expect("uniqueness worker panicked")?;
        }
        Ok(())
    })?;

    let mut stats = UniquenessStats { genomes: genomes.len(), kmers_looked_up: looked_up.into_inner(), ..Default::default() };
    table_keys::<C, CPB>(&flexmap.keys.data)(&mut |_key, start, slots| {
        let header = slots > T as u64;
        let h = if header { (slots + 2) / 3 } else { 0 };
        for c in &flexmap.values.data[(start + h) as usize..(start + slots) as usize] {
            let f = c.flags();
            stats.values += 1;
            if header {
                stats.long_values += 1;
                stats.long_unique += (f & UNIQUE != 0) as u64;
                stats.long_super_unique += (f == UNIQUE | SUPER_UNIQUE) as u64;
            } else {
                stats.short_values += 1;
                stats.short_unique += (f & UNIQUE != 0) as u64;
            }
        }
    });
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;
    use crate::build::default_build;
    use crate::flexmap::VRangeGetter;
    use crate::packed::{CompactKeys, PackedFlexmap};
    use crate::range::AnyRange;

    const K: usize = 19;
    const C: usize = 11;
    const F: usize = 16;
    const S: usize = 3;
    const L: usize = C - S + 1;
    const CPB: u64 = 8;
    const T: usize = 2;

    fn random_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..n)
            .map(|_| {
                state ^= state >> 12;
                state ^= state << 25;
                state ^= state >> 27;
                b"ACGT"[(state.wrapping_mul(0x2545F4914F6CDD1D) >> 62) as usize]
            })
            .collect()
    }

    /// `seq` with the base at every `step`-th position from `from` changed.
    fn mutated(seq: &[u8], from: usize, step: usize) -> Vec<u8> {
        let mut s = seq.to_vec();
        for i in (from..s.len()).step_by(step) {
            s[i] = if s[i] == b'A' { b'C' } else { b'A' };
        }
        s
    }

    fn write_fasta(dir: &Path, name: &str, records: &[(&str, &[u8])]) -> PathBuf {
        let path = dir.join(name);
        let mut f = File::create(&path).unwrap();
        for (head, seq) in records {
            writeln!(f, ">{head}").unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
        }
        path
    }

    fn occurs(text: &[u8], s: &[u8]) -> bool {
        let rc = revcomp(s);
        text.windows(s.len()).any(|w| w == s || w == &rc[..])
    }

    /// Species A has three genomes in the index (so its keys have headers: long values), species B
    /// one (short values). A genome of B carries an exact 2 kb stretch of A (those long values lose
    /// UNIQUE) and a 2 kb stretch of A with a mutation every K bases (those whose mutation falls in
    /// a flank lose SUPER_UNIQUE only).
    #[test]
    fn genomes_of_other_species_clear_flags() {
        let dir = std::env::temp_dir().join(format!("flexmap_uniqueness_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let rep_a = random_seq(30_000, 0x1234_5678_9abc_def1);
        let a2 = mutated(&rep_a, 50, 97);
        let a3 = mutated(&rep_a, 20, 89);
        let rep_b = random_seq(30_000, 0x0fed_cba9_8765_4321);
        let exact = rep_a[10_000..12_000].to_vec();
        let near = mutated(&rep_a[20_000..22_000], 0, K);

        let reference = write_fasta(&dir, "ref.fa", &[("repA", &rep_a), ("a2", &a2), ("a3", &a3), ("repB", &rep_b)]);
        let (mut flexmap, _ref2id, id2ref) =
            default_build::<K, C, F, S, L, CPB, T>(&reference, 10_000).unwrap();
        let mut species = SpeciesTable::default();
        let seq_species: HashMap<String, String> = [("repA", "A"), ("a2", "A"), ("a3", "A"), ("repB", "B")]
            .iter()
            .map(|(a, b)| (a.to_string(), b.to_string()))
            .collect();
        let ref_species = reference_species(&id2ref, &seq_species, &mut species).unwrap();
        let (sa, sb) = (species.id("A"), species.id("B"));

        let genomes = vec![
            (write_fasta(&dir, "a1.fa", &[("a1", &rep_a)]), sa),
            (write_fasta(&dir, "a2.fa", &[("a2", &a2)]), sa),
            (write_fasta(&dir, "a3.fa", &[("a3", &a3)]), sa),
            (write_fasta(&dir, "b1.fa", &[("b1", &rep_b)]), sb),
            (write_fasta(&dir, "b2.fa", &[("b2_exact", &exact), ("b2_near", &near)]), sb),
        ];
        let stats = mark_uniqueness::<K, C, F, S, L, CPB, T>(&mut flexmap, &ref_species, &genomes, 3).unwrap();
        assert_eq!(stats.genomes, 5);
        assert!(stats.long_values > 1000 && stats.short_values > 1000, "{stats:?}");
        assert!(stats.long_super_unique < stats.long_unique && stats.long_unique < stats.long_values, "{stats:?}");

        let b_text: Vec<u8> = [&rep_b[..], b"N", &exact[..], b"N", &near[..]].concat();
        let a_text: Vec<u8> = [&rep_a[..], b"N", &a2[..], b"N", &a3[..]].concat();
        let mut cs = ClosedSyncmer::<K, S, L>::new();
        let (mut exact_seen, mut near_seen, mut short_seen) = (0, 0, 0);
        for (rid, seq, other_text) in [(1u64, &rep_a, &b_text), (4u64, &rep_b, &a_text)] {
            for (pos, f, r) in KmerIter::<K, true>::new(seq) {
                let cmer = min(f.middle::<C>(), r.middle::<C>());
                if !cs.is_minimizer(cmer.0) {
                    continue;
                }
                let range = AnyRange::from(flexmap.get_vrange(mask_coremer::<C>(cmer.0)).unwrap());
                let i = (0..range.len()).find(|&i| range.value(i) == (rid, pos as u64)).unwrap();
                let flags = range.flags(i);
                let window = &seq[pos..pos + K];
                if range.has_header() {
                    let in_other = occurs(other_text, window);
                    // Unique exactly when no genome of the other species has the k-mer.
                    assert_eq!(flags & UNIQUE != 0, !in_other, "long value {rid}:{pos}");
                    if rid == 1 && (10_000..=12_000 - K).contains(&pos) {
                        assert_eq!(flags, 0, "exact stretch {pos}");
                        exact_seen += 1;
                    }
                    if rid == 1 && (20_000..=22_000 - K).contains(&pos) {
                        let m = (K - (pos - 20_000) % K) % K; // the mutated base of this window
                        let in_flank = m < (K - C) / 2 || m >= (K + C) / 2;
                        if in_flank && !in_other {
                            assert_eq!(flags, UNIQUE, "near stretch {pos}");
                            near_seen += 1;
                        }
                    }
                } else {
                    // A short value is matched on its core-mer alone.
                    let core = &window[(K - C) / 2..(K + C) / 2];
                    if !occurs(other_text, core) {
                        assert_eq!(flags, UNIQUE, "short value {rid}:{pos}");
                        short_seen += 1;
                    }
                }
            }
        }
        assert!(exact_seen > 50 && near_seen > 50 && short_seen > 1000, "{exact_seen} {near_seen} {short_seen}");

        // The packed index carries the same flags, through a file and a slice; an index without the
        // pass has none.
        let packed = PackedFlexmap::<CompactKeys, C, F, T>::from_flexmap(&flexmap);
        assert!(packed.has_flags());
        for key in 0..(1u64 << (2 * C)) {
            match (flexmap.get_vrange(key), packed.get(key)) {
                (Some(a), Some(b)) => {
                    for i in 0..b.len() {
                        assert_eq!((a.flags(i), VD::get(a.positions[i].0)), (b.flags(i), b.get(i)), "key {key} value {i}");
                    }
                }
                (None, None) => {}
                _ => panic!("key {key}"),
            }
        }
        let path = dir.join("flagged.compact").to_string_lossy().to_string();
        packed.save(&path).unwrap();
        assert!(PackedFlexmap::<CompactKeys, C, F, T>::is_compatible(&path));
        let loaded = PackedFlexmap::<CompactKeys, C, F, T>::load(&path, 2).unwrap();
        assert!(loaded.has_flags());
        let (lo, hi) = (1000u64, 3_000_000u64);
        let slice = packed.slice(lo, hi);
        assert!(slice.has_flags());
        for key in lo..hi {
            if let Some(a) = packed.get(key) {
                let (b, c) = (loaded.get(key).unwrap(), slice.get(key - lo).unwrap());
                for i in 0..a.len() {
                    assert_eq!((a.flags(i), a.get(i)), (b.flags(i), b.get(i)));
                    assert_eq!((a.flags(i), a.get(i)), (c.flags(i), c.get(i)));
                }
            }
        }

        let (plain, _, _) = default_build::<K, C, F, S, L, CPB, T>(&reference, 10_000).unwrap();
        let plain = PackedFlexmap::<CompactKeys, C, F, T>::from_flexmap(&plain);
        assert!(!plain.has_flags());
        assert_eq!(plain.values.layout.value_bits(), packed.values.layout.value_bits());
        assert_eq!(plain.values.layout.entry_bits() + crate::packed::FLAG_BITS, packed.values.layout.entry_bits());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect()
    }
}
