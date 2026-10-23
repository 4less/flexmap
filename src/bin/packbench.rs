//! Benchmarks of the Flexmap blob (as flexalign loads it) against the packed flexmap: memory, load
//! time and lookup speed on a real reference, with reads simulated from it.
//!
//!   packbench build <fasta> <out-prefix> <max-range> [old,compact,direct]
//!       builds the Flexmap (K 31, C 15, F 16, S 4, CPB 16, threshold 2, as flexalign), converts it,
//!       writes <out-prefix>.old.blob / .compact.pk / .direct.pk (the ones listed; all by default)
//!   packbench reads <fasta> <out.q> <reads> <error-rate> <seed> [read-length (150)]
//!       simulates reads and writes the (key, flank) of each of their syncmer k-mers, read by read
//!   packbench query ... with FLEXMAP_LONG=<long.q>: also times long reads (a `reads` file of
//!       long pieces), each read one call
//!   packbench query <old|compact|direct> <heap|mmap> <index> <queries.q> <threads,..> <rounds>
//!       loads the index and times lookups; prints one TSV row per measurement
use std::cmp::min;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::time::Instant;

use bioreader::fasta_byte_reader::FastaByteReader;
use bioreader::fasta_reader::FastaReader;
use bioreader::sequence::fasta_record::OwnedFastaRecord;
use flexmap::build::default_build;
use flexmap::flexmap::{Flexmap, FlexmapBlob, VRangeGetter};
use flexmap::packed::{CompactKeys, DirectKeys, HashedKeys, PackedFlexmap, SlotLocator};
use kmerrs::consecutive::kmer::{Kmer, KmerIter};
use kmerrs::kmer_utils::mask_coremer;
use kmerrs::minimizer::context_free::Minimizer;
use kmerrs::syncmer::closed_syncmer::ClosedSyncmer;
use std::sync::{Arc, Mutex};

const K: usize = 31;
const C: usize = 15;
const F: usize = 16;
const S: usize = 4;
const L: usize = C - S + 1;
const CPB: u64 = 16;
const T: usize = 2;

type Old = Flexmap<C, F, CPB, T>;
type OldBlob = FlexmapBlob<C, F, CPB, T>;
type Compact = PackedFlexmap<CompactKeys, C, F, T>;
type Direct = PackedFlexmap<DirectKeys<C, CPB>, C, F, T>;
type Hashed = PackedFlexmap<HashedKeys, C, F, T>;

fn status(field: &str) -> u64 {
    let s = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    s.lines()
        .find(|l| l.starts_with(field))
        .and_then(|l| l.split_whitespace().nth(1))
        .and_then(|v| v.parse::<u64>().ok())
        .map(|kb| kb * 1024)
        .unwrap_or(0)
}

fn anon_huge() -> u64 {
    let s = std::fs::read_to_string("/proc/self/smaps_rollup").unwrap_or_default();
    s.lines()
        .find(|l| l.starts_with("AnonHugePages:"))
        .and_then(|l| l.split_whitespace().nth(1))
        .and_then(|v| v.parse::<u64>().ok())
        .map(|kb| kb * 1024)
        .unwrap_or(0)
}

fn gb(x: u64) -> f64 {
    x as f64 / 1e9
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(|s| s.as_str()) {
        Some("build") => build(&args[2..]),
        Some("reads") => reads(&args[2..]),
        Some("query") => query(&args[2..]),
        Some("hash") => {
            // packbench hash <compact.pk> <out.hashed.pk>: the hashed key table from a compact index.
            let compact = Compact::load(&args[2], 8).unwrap();
            let t = Instant::now();
            let hashed = Hashed::from_compact(&compact).expect("hashed table fits");
            report_packed("hashed", &hashed.keys, &hashed.values, t.elapsed().as_secs_f64());
            hashed.save(&args[3]).unwrap();
        }
        _ => {
            eprintln!("usage: see the top of src/bin/packbench.rs");
            std::process::exit(2);
        }
    }
}

fn build(args: &[String]) {
    let (fasta, out, max_range) = (&args[0], &args[1], args[2].parse::<usize>().unwrap());
    let which = args.get(3).cloned().unwrap_or_else(|| "old,compact,direct".into());
    let t = Instant::now();
    let (flexmap, _, id2ref) = default_build::<K, C, F, S, L, CPB, T>(fasta, max_range).unwrap();
    let build_s = t.elapsed().as_secs_f64();
    let (kb, vb) = (flexmap.keys.data.len() * 2, flexmap.values.data.len() * 8);
    let used = (0..1u64 << (2 * C)).filter(|&k| flexmap.keys.vrange(k).is_some()).count();
    println!("build\t{fasta}\trefs\t{}\tused_keys\t{used}\tcells\t{}\tbuild_s\t{build_s:.1}\tpeak_rss_gb\t{:.2}", id2ref.len() - 1, flexmap.values.data.len(), gb(status("VmHWM:")));
    println!("old\tkeys_gb\t{:.3}\tvalues_gb\t{:.3}\ttotal_gb\t{:.3}", gb(kb as u64), gb(vb as u64), gb((kb + vb) as u64));
    if which.contains("old") {
        let t = Instant::now();
        flexmap.save_blob(&format!("{out}.old.blob"));
        println!("old\tsave_s\t{:.1}", t.elapsed().as_secs_f64());
    }
    if which.contains("compact") {
        let t = Instant::now();
        let packed = Compact::from_flexmap(&flexmap);
        let convert_s = t.elapsed().as_secs_f64();
        report_packed("compact", &packed.keys, &packed.values, convert_s);
        packed.save(&format!("{out}.compact.pk")).unwrap();
    }
    if which.contains("hashed") {
        let t = Instant::now();
        match Hashed::from_flexmap_hashed(&flexmap) {
            Some(packed) => {
                let convert_s = t.elapsed().as_secs_f64();
                report_packed("hashed", &packed.keys, &packed.values, convert_s);
                packed.save(&format!("{out}.hashed.pk")).unwrap();
            }
            None => println!("hashed\tdoes not fit"),
        }
    }
    if which.contains("direct") {
        let t = Instant::now();
        let packed = Direct::from_flexmap_direct(&flexmap);
        let convert_s = t.elapsed().as_secs_f64();
        report_packed("direct", &packed.keys, &packed.values, convert_s);
        packed.save(&format!("{out}.direct.pk")).unwrap();
    }
}

fn report_packed<KM: SlotLocator>(name: &str, keys: &KM, values: &flexmap::packed::PackedValues, convert_s: f64) {
    let l = values.layout;
    println!(
        "{name}\tkeys_gb\t{:.3}\tvalues_gb\t{:.3}\ttotal_gb\t{:.3}\tref_bits\t{}\tpos_bits\t{}\tbits_per_cell\t{:.2}\tconvert_s\t{convert_s:.1}",
        gb(keys.memory_bytes() as u64),
        gb(values.memory_bytes() as u64),
        gb((keys.memory_bytes() + values.memory_bytes()) as u64),
        l.ref_bits,
        l.pos_bits,
        l.bits_per_slot()
    );
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545F4914F6CDD1D)
    }

    fn below(&mut self, n: u64) -> u64 {
        ((self.next() as u128 * n as u128) >> 64) as u64
    }
}

fn reads(args: &[String]) {
    let (fasta, out) = (&args[0], &args[1]);
    let (n_reads, err, seed) = (args[2].parse::<usize>().unwrap(), args[3].parse::<f64>().unwrap(), args[4].parse::<u64>().unwrap());
    let len = args.get(5).map(|l| l.parse::<usize>().unwrap()).unwrap_or(150);
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let file = File::open(fasta).unwrap();
    let mut byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, 1 << 24).unwrap()));
    let mut reader = FastaReader::with_capacity(1 << 24);
    let mut record = OwnedFastaRecord::new();
    while let Some(()) = reader.load_batch_par(&mut byte_reader).unwrap() {
        while let Some(_) = reader.next(&mut record) {
            if record.seq().len() >= len + 1000 {
                seqs.push(record.seq().to_vec());
            }
        }
    }
    let mut ends = Vec::with_capacity(seqs.len());
    let mut total = 0u64;
    for s in &seqs {
        total += (s.len() - len) as u64;
        ends.push(total);
    }
    let mut rng = Rng(seed | 1);
    let mut cs = ClosedSyncmer::<C, S, L>::new();
    let (mut keys, mut flexes) = (Vec::new(), Vec::new());
    let mut offsets = vec![0u64];
    let mut made = 0;
    let mut read = vec![0u8; len];
    while made < n_reads {
        let at = rng.below(total);
        let r = ends.partition_point(|&e| e <= at);
        let start = (at - if r == 0 { 0 } else { ends[r - 1] }) as usize;
        let src = &seqs[r][start..start + len];
        if src.iter().any(|&c| !matches!(c, b'A' | b'C' | b'G' | b'T')) {
            continue;
        }
        let rc = rng.next() & 1 == 1;
        for i in 0..len {
            let c = if rc {
                match src[len - 1 - i] { b'A' => b'T', b'C' => b'G', b'G' => b'C', _ => b'A' }
            } else {
                src[i]
            };
            read[i] = if (rng.next() as f64 / u64::MAX as f64) < err { b"ACGT"[rng.below(4) as usize] } else { c };
        }
        for (_, kf, kr) in KmerIter::<K, true>::new(&read) {
            let (cf, cr) = (kf.middle::<C>(), kr.middle::<C>());
            let kmer: Kmer<K> = if cf < cr { kf } else { kr };
            let cmer = min(cf, cr);
            if !cs.is_minimizer(cmer.0) {
                continue;
            }
            keys.push(mask_coremer::<C>(cmer.0));
            flexes.push(kmer.flanks::<F>().0 as u32);
        }
        offsets.push(keys.len() as u64);
        made += 1;
    }
    let mut w = BufWriter::new(File::create(out).unwrap());
    w.write_all(&(keys.len() as u64).to_le_bytes()).unwrap();
    w.write_all(bytemuck::cast_slice(&keys)).unwrap();
    w.write_all(bytemuck::cast_slice(&flexes)).unwrap();
    w.write_all(&(offsets.len() as u64).to_le_bytes()).unwrap();
    w.write_all(bytemuck::cast_slice(&offsets)).unwrap();
    println!("reads\t{fasta}\t{n_reads}\tlength\t{len}\terr\t{err}\tqueries\t{}", keys.len());
}

/// A query file: the queries and where each read's start (reads.len() = reads + 1).
fn read_queries(path: &str) -> (Vec<u64>, Vec<u32>, Vec<u64>) {
    let mut f = File::open(path).unwrap();
    let mut n = [0u8; 8];
    f.read_exact(&mut n).unwrap();
    let n = u64::from_le_bytes(n) as usize;
    let mut keys = vec![0u64; n];
    let mut flexes = vec![0u32; n];
    f.read_exact(bytemuck::cast_slice_mut(&mut keys)).unwrap();
    f.read_exact(bytemuck::cast_slice_mut(&mut flexes)).unwrap();
    let mut m = [0u8; 8];
    f.read_exact(&mut m).unwrap();
    let mut reads = vec![0u64; u64::from_le_bytes(m) as usize];
    f.read_exact(bytemuck::cast_slice_mut(&mut reads)).unwrap();
    (keys, flexes, reads)
}

#[inline(always)]
fn mix(i: usize, pos: u64, val: u64, d: Option<(u32, u32)>) -> u64 {
    let (a, b) = d.map(|(a, b)| (a as u64 + 1, b as u64)).unwrap_or((0, 0));
    let mut h = (i as u64).wrapping_mul(0x9E3779B97F4A7C15) ^ pos.wrapping_mul(0xC2B2AE3D27D4EB4F) ^ val.wrapping_mul(0x165667B19E3779F9);
    h ^= (a << 32 | b).wrapping_mul(0xD6E8FEB86659FD93);
    h ^ (h >> 31)
}

/// One way of answering a run of queries: returns (checksum, values seen).
trait Bench: Sync {
    fn run(&self, work: Work, prefetch: bool, base: usize, keys: &[u64], flexes: &[u32]) -> (u64, u64);
}

#[derive(Clone, Copy, PartialEq)]
enum Work {
    /// The key's values, every one decoded.
    Values,
    /// best_flex_match: the values of the nearest flanks.
    Seeds,
}

impl Bench for OldBlob {
    fn run(&self, work: Work, prefetch: bool, base: usize, keys: &[u64], flexes: &[u32]) -> (u64, u64) {
        let (mut sum, mut seen) = (0u64, 0u64);
        let mut one = |i: usize, vr: Option<flexmap::values::VRange>| {
            if let Some(vr) = vr {
                let q = base + i;
                match work {
                    Work::Values => vr.all_matches(|p, v| {
                        sum = sum.wrapping_add(mix(q, p, v, None));
                        seen += 1;
                    }),
                    Work::Seeds => vr.best_flex_match(&Kmer::<F>(flexes[i] as u64), |p, v, d| {
                        sum = sum.wrapping_add(mix(q, p, v, d));
                        seen += 1;
                    }),
                }
            }
        };
        if prefetch {
            let mut i = 0;
            self.for_each_vrange_prefetched(keys, |_, vr| {
                one(i, vr);
                i += 1;
            });
        } else {
            for (i, &k) in keys.iter().enumerate() {
                one(i, self.get_vrange(k));
            }
        }
        (sum, seen)
    }
}

impl<KM: SlotLocator> Bench for PackedFlexmap<KM, C, F, T> {
    fn run(&self, work: Work, prefetch: bool, base: usize, keys: &[u64], flexes: &[u32]) -> (u64, u64) {
        let (mut sum, mut seen) = (0u64, 0u64);
        let mut one = |i: usize, r: Option<flexmap::packed::PackedRange>| {
            if let Some(r) = r {
                let q = base + i;
                match work {
                    Work::Values => r.all_matches(|p, v| {
                        sum = sum.wrapping_add(mix(q, p, v, None));
                        seen += 1;
                    }),
                    Work::Seeds => r.best_flex_match(&Kmer::<F>(flexes[i] as u64), |p, v, d| {
                        sum = sum.wrapping_add(mix(q, p, v, d));
                        seen += 1;
                    }),
                }
            }
        };
        if prefetch {
            let step = std::env::var("FLEXMAP_PF_STEP").ok().and_then(|s| s.parse().ok()).unwrap_or(8);
            self.for_each_prefetched_with(keys, step, |i, r| one(i, r));
        } else {
            for (i, &k) in keys.iter().enumerate() {
                one(i, self.get(k));
            }
        }
        (sum, seen)
    }
}

/// Runs `bench` over the queries on `threads` threads (contiguous parts); wall seconds, checksum, values.
/// With `reads`, each read's queries are one call (an aligner's read at a time) and the threads take
/// whole reads; without, each thread's part is one call.
fn timed(bench: &dyn Bench, work: Work, prefetch: bool, keys: &[u64], flexes: &[u32], reads: Option<&[u64]>, threads: usize) -> (f64, u64, u64) {
    let t = Instant::now();
    let results: Vec<(u64, u64)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|p| {
                scope.spawn(move || match reads {
                    None => {
                        let part = (keys.len() + threads - 1) / threads;
                        let from = (p * part).min(keys.len());
                        let to = ((p + 1) * part).min(keys.len());
                        bench.run(work, prefetch, from, &keys[from..to], &flexes[from..to])
                    }
                    Some(reads) => {
                        let n = reads.len() - 1;
                        let part = (n + threads - 1) / threads;
                        let (mut sum, mut seen) = (0u64, 0u64);
                        for r in (p * part).min(n)..((p + 1) * part).min(n) {
                            let (from, to) = (reads[r] as usize, reads[r + 1] as usize);
                            let (a, b) = bench.run(work, prefetch, from, &keys[from..to], &flexes[from..to]);
                            sum = sum.wrapping_add(a);
                            seen += b;
                        }
                        (sum, seen)
                    }
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let s = t.elapsed().as_secs_f64();
    let sum = results.iter().fold(0u64, |a, r| a.wrapping_add(r.0));
    let seen = results.iter().map(|r| r.1).sum();
    (s, sum, seen)
}

fn query(args: &[String]) {
    let (kind, mode, index, qpath) = (&args[0], &args[1], &args[2], &args[3]);
    let threads: Vec<usize> = args[4].split(',').map(|t| t.parse().unwrap()).collect();
    let rounds: usize = args[5].parse().unwrap();
    // FLEXMAP_QUERIES: use the first reads up to this many queries (all by default).
    let limit = std::env::var("FLEXMAP_QUERIES").ok().and_then(|s| s.parse::<usize>().ok());
    let load_queries = |path: &str| {
        let (mut keys, mut flexes, mut reads) = read_queries(path);
        if let Some(n) = limit {
            let r = reads.partition_point(|&o| o as usize <= n);
            reads.truncate(r.max(1));
            let end = *reads.last().unwrap() as usize;
            keys.truncate(end);
            flexes.truncate(end);
        }
        (keys, flexes, reads)
    };
    let (keys, flexes, reads) = load_queries(qpath);
    let long = std::env::var("FLEXMAP_LONG").ok().map(|p| load_queries(&p));
    let mut rng = Rng(0x5DEECE66D);
    let random: Vec<u64> = (0..keys.len()).map(|_| rng.next() & ((1u64 << (2 * C)) - 1)).collect();
    let rss0 = status("VmRSS:");
    let load_threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
    let t = Instant::now();
    let bench: Box<dyn Bench> = match (kind.as_str(), mode.as_str()) {
        ("old", "heap") => Box::new(Old::load_blob(index)),
        ("old", "mmap") => Box::new(OldBlob::mmap_from_file(index)),
        ("compact", "heap") => Box::new(Compact::load(index, load_threads).unwrap()),
        ("compact", "mmap") => Box::new(Compact::mmap(index).unwrap()),
        ("hashed", "heap") => Box::new(Hashed::load(index, load_threads).unwrap()),
        ("hashed", "mmap") => Box::new(Hashed::mmap(index).unwrap()),
        ("direct", "heap") => Box::new(Direct::load(index, load_threads).unwrap()),
        ("direct", "mmap") => Box::new(Direct::mmap(index).unwrap()),
        _ => panic!("unknown index kind {kind} {mode}"),
    };
    let load_s = t.elapsed().as_secs_f64();
    let thp = std::env::var("FLEXMAP_THP").map(|v| v != "0").unwrap_or(true);
    let file_gb = gb(std::fs::metadata(index).unwrap().len());
    let label = format!("{kind}\t{mode}\tthp={}", thp as u8);
    println!(
        "load\t{label}\tfile_gb\t{file_gb:.3}\tload_s\t{load_s:.2}\trss_gb\t{:.3}\thuge_gb\t{:.3}",
        gb(status("VmRSS:").saturating_sub(rss0)),
        gb(anon_huge())
    );
    let none: (Vec<u64>, Vec<u32>, Vec<u64>) = Default::default();
    let lg = long.as_ref().unwrap_or(&none);
    // (name, keys, flanks, reads if one call per read, work)
    let mut sets: Vec<(&str, &[u64], &[u32], Option<&[u64]>, Work)> = vec![
        ("reads", &keys, &flexes, None, Work::Seeds),
        ("perread", &keys, &flexes, Some(&reads), Work::Seeds),
    ];
    if long.is_some() {
        sets.push(("long", &lg.0, &lg.1, Some(&lg.2), Work::Seeds));
    }
    sets.push(("reads", &keys, &flexes, None, Work::Values));
    sets.push(("random", &random, &flexes, None, Work::Values));
    // FLEXMAP_SETS: only these "set:work" pairs (e.g. perread:seeds,long:seeds).
    let only = std::env::var("FLEXMAP_SETS").ok();
    for (set, ks, fs, rs, work) in sets {
        let work_name = if work == Work::Seeds { "seeds" } else { "values" };
        if let Some(o) = &only {
            if !o.split(',').any(|x| x == format!("{set}:{work_name}")) {
                continue;
            }
        }
        for &th in &threads {
            for prefetch in [false, true] {
                let mut times = Vec::new();
                let mut first = None;
                let mut check = (0, 0);
                for r in 0..rounds + 1 {
                    let (s, sum, seen) = timed(bench.as_ref(), work, prefetch, ks, fs, rs, th);
                    if r == 0 {
                        first = Some(s);
                    } else {
                        times.push(s);
                    }
                    check = (sum, seen);
                }
                times.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let median = times[times.len() / 2];
                println!(
                    "time\t{label}\t{set}\t{work_name}\tthreads\t{th}\tprefetch\t{}\tqueries\t{}\tmedian_s\t{median:.4}\tfirst_s\t{:.4}\tMq_per_s\t{:.2}\tns_per_q_thread\t{:.1}\tvalues\t{}\tchecksum\t{:016x}",
                    prefetch as u8,
                    ks.len(),
                    first.unwrap(),
                    ks.len() as f64 / median / 1e6,
                    median * th as f64 / ks.len() as f64 * 1e9,
                    check.1,
                    check.0
                );
            }
        }
    }
    println!("end\t{label}\trss_gb\t{:.3}\tpeak_gb\t{:.3}", gb(status("VmRSS:")), gb(status("VmHWM:")));
}
