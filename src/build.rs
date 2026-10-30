use std::{cmp::min, collections::HashMap, fs::File, io::BufRead, path::Path, sync::{Arc, Mutex}};

use kmerrs::{consecutive::kmer::KmerIter, kmer_utils::mask_coremer, minimizer::context_free::Minimizer, syncmer::closed_syncmer::ClosedSyncmer};
use bioreader::{fasta_byte_reader::{self, FastaByteReader}, fasta_reader::{self, FastaReader}, fastq_byte_reader, fastq_reader, sequence::fasta_record::OwnedFastaRecord};
use savefile::save;

use crate::{flexmap::Flexmap, keys::FMKeys, VD};

pub fn default_build<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize
>(path: impl AsRef<Path>, max_range_size: usize) -> 
        Result<(Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>, HashMap<String, usize>, Vec<String>), std::io::Error> {
    eprintln!("Build keys");
    let keys = match default_build_keys::<K, C, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>(&path, max_range_size) {
        Ok(keys) => keys,
        Err(_) => panic!("Keys could not be built."),
    };

    eprintln!("Build map");
    let (map, reference2id, id2reference) = match default_build_map::<K, C, F, S, L, CELLS_PER_BODY, HEADER_THRESHOLD>(path, keys) {
        Ok(map) => map,
        Err(_) => panic!("Map could not be built."),
    };

    Ok((map, reference2id, id2reference))
}


fn default_build_keys<
    const K: usize,
    const C: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(path: impl AsRef<Path>, max_range_size: usize) -> Result<FMKeys<C, CELLS_PER_BODY>, std::io::Error> {
    let mut keys = FMKeys::<C, CELLS_PER_BODY>::new();

    let file: File = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.as_ref().display(), why),
        Ok(file) => file,
    };

    let mut record = OwnedFastaRecord::new();
    
    let buffer_size = usize::pow(2, 24);
    let mut byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size)?));
    let mut fasta_reader = FastaReader::with_capacity(buffer_size);
    let mut cs = ClosedSyncmer::<C,S,L>::new();

    eprintln!("read data");
    while let Some(()) = fasta_reader
        .load_batch_par(&mut byte_reader)
        .expect("Batch is invalid")
    {
        while let Some(_) = fasta_reader.next(&mut record) {
            if !record.valid_extended() {
                panic!("Record is not valid {:?}", record.to_string())
            }

            let mut iter = KmerIter::<K, true>::new(record.seq());

            for (_, kmer_fwd, kmer_rev) in iter {
                let cmer_fwd = kmer_fwd.middle::<C>();
                let cmer_rev = kmer_rev.middle::<C>();
                let _kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
                let cmer = min(cmer_fwd, cmer_rev);

                if !cs.is_minimizer(cmer.0) { continue };

                // Selection stays on the raw canonical core-mer (sampling unchanged); only the
                // direct-addressed KEY is mixed, so the populated keys spread uniformly over the key
                // space instead of piling at the low end -- which is what lets shards slice evenly.
                let key = mask_coremer::<C>(cmer.0);
                keys.get_kmer_cell_mut_ref(key).increment();
            }
        }
    }

    eprintln!("Keys build {} {}", HEADER_THRESHOLD, max_range_size);
    keys.build::<HEADER_THRESHOLD>(max_range_size);

    Ok(keys)
}




fn default_build_map<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(path: impl AsRef<Path>, keys: FMKeys<C, CELLS_PER_BODY>) -> 
        Result<(Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>, HashMap<String, usize>, Vec<String>), std::io::Error> {

    let mut cs = ClosedSyncmer::<K,S,L>::new();
    let file: File = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.as_ref().display(), why),
        Ok(file) => file,
    };

    let mut flexmap = Flexmap::<C,F,CELLS_PER_BODY,HEADER_THRESHOLD>::new(keys);

    let mut record = OwnedFastaRecord::new();
    
    let buffer_size = usize::pow(2, 24);
    let mut byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size)?));
    let mut fasta_reader = FastaReader::with_capacity(buffer_size);

    let mut reference_running_id = 1usize;
    let mut reference2id: HashMap<String, usize> = HashMap::<String, usize>::new();
    let mut id2reference: Vec<String> = Vec::<String>::new();
    id2reference.push("dummy".into());

    let mut total_kmers = 0u64;
    let mut total_minimizers = 0u64;
    
    while let Some(()) = fasta_reader
        .load_batch_par(&mut byte_reader)
        .expect("Batch is invalid")
    {

        while let Some(_) = fasta_reader.next(&mut record) {
            if !record.valid_extended() {
                panic!("Record is not valid {:?}", record.to_string())
            }
            let iter = KmerIter::<K, true>::new(record.seq());

            let header = String::from_utf8_lossy(&record.head()[1..]).split(' ').next().unwrap().to_string();

            let reference_id = match reference2id
                    .try_insert(header.clone(), reference_running_id) {
                Ok(reference_id) => *reference_id,
                Err(_) => panic!("Header {:?} has been seen before!", record.head()),
            };
            id2reference.push(header.clone());
            assert_eq!(id2reference[reference_id], header);

            reference_running_id += 1;


            // println!("{:?} -> {}", String::from_utf8_lossy(record.head()), reference_id);

            for (pos, kmer_fwd, kmer_rev) in iter {
                let cmer_fwd = kmer_fwd.middle::<C>();
                let cmer_rev = kmer_rev.middle::<C>();
                let kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
                let cmer = min(cmer_fwd, cmer_rev);

                total_kmers += 1;

                if !cs.is_minimizer(cmer.0) { continue };

                total_minimizers += 1;

                let flanks = kmer.flanks::<F>();

                // Same key mix as the keys pass -- must be identical so values land in the slots
                // counted there.
                let key = mask_coremer::<C>(cmer.0);
                match flexmap.keys.vrange(key) {
                    Some(range) => {
                        let mut vblock = flexmap.values.get_range_mut(range);
                        vblock.insert(VD::set(reference_id as u64, pos as u64), flanks.0 as u32);
                    },
                    None => {},
                };

                // if reference_id == 1 && pos == 578574;
            }
        }

    }
    eprintln!("Minimizer compression rate: {} ({}/{})", total_minimizers as f64/total_kmers as f64, total_minimizers, total_kmers);

    eprintln!("Number of ids: {}", id2reference.len());

    Ok((flexmap, reference2id, id2reference))
}


/// One piece of a batch of records: positions `start..end` of record `rec` (a long record is cut
/// into pieces; a piece reads K - 1 bases past its end so no k-mer starting in it is lost).
struct Piece {
    rec: usize,
    start: usize,
    end: usize,
}

/// Records are read into batches of about this many bases, and cut into pieces of at most `PIECE`.
const BATCH_BASES: usize = 128 << 20;
const PIECE: usize = 4 << 20;

/// Reads the FASTA in batches of whole records: (header, sequence) per record, in file order.
fn for_each_batch(path: impl AsRef<Path>, mut f: impl FnMut(Vec<(String, Vec<u8>)>)) -> Result<(), std::io::Error> {
    let file: File = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.as_ref().display(), why),
        Ok(file) => file,
    };
    let buffer_size = usize::pow(2, 24);
    let mut byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size)?));
    let mut fasta_reader = FastaReader::with_capacity(buffer_size);
    let mut record = OwnedFastaRecord::new();
    let (mut batch, mut bases) = (Vec::new(), 0usize);
    while let Some(()) = fasta_reader.load_batch_par(&mut byte_reader).expect("Batch is invalid") {
        while let Some(_) = fasta_reader.next(&mut record) {
            if !record.valid_extended() {
                panic!("Record is not valid {:?}", record.to_string())
            }
            let header = String::from_utf8_lossy(&record.head()[1..]).split(' ').next().unwrap().to_string();
            bases += record.seq().len();
            batch.push((header, record.seq().to_vec()));
            if bases >= BATCH_BASES {
                f(std::mem::take(&mut batch));
                bases = 0;
            }
        }
    }
    if !batch.is_empty() {
        f(batch);
    }
    Ok(())
}

/// The k-mers of one batch the index keeps (closed-syncmer cores), extracted on `threads` threads
/// and split by key range: `out[piece][range]` holds, in sequence order, (key, record, position,
/// flanks) of each kept k-mer whose key falls in `range`. `range_of(key)` names the range.
fn extract_batch<const K: usize, const C: usize, const F: usize, const S: usize, const L: usize>(
    batch: &[(String, Vec<u8>)], threads: usize, ranges: usize, range_of: &(dyn Fn(u64) -> usize + Sync),
) -> Vec<Vec<Vec<(u64, u32, u32, u32)>>> {
    let mut pieces = Vec::new();
    for (rec, (_, seq)) in batch.iter().enumerate() {
        let mut start = 0;
        while start < seq.len() {
            let end = (start + PIECE).min(seq.len());
            pieces.push(Piece { rec, start, end });
            start = end;
        }
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let slots: Vec<Mutex<Vec<Vec<(u64, u32, u32, u32)>>>> = (0..pieces.len()).map(|_| Mutex::new(Vec::new())).collect();
    std::thread::scope(|scope| {
        for _ in 0..threads.max(1) {
            scope.spawn(|| {
                let mut cs = ClosedSyncmer::<C, S, L>::new();
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if i >= pieces.len() {
                        break;
                    }
                    let p = &pieces[i];
                    let seq = &batch[p.rec].1;
                    let window = &seq[p.start..(p.end + K - 1).min(seq.len())];
                    let mut out: Vec<Vec<(u64, u32, u32, u32)>> = vec![Vec::new(); ranges];
                    for (pos, kmer_fwd, kmer_rev) in KmerIter::<K, true>::new(window) {
                        if pos >= p.end - p.start {
                            break;
                        }
                        let cmer_fwd = kmer_fwd.middle::<C>();
                        let cmer_rev = kmer_rev.middle::<C>();
                        let kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
                        let cmer = min(cmer_fwd, cmer_rev);
                        if !cs.is_minimizer(cmer.0) {
                            continue;
                        }
                        let key = mask_coremer::<C>(cmer.0);
                        out[range_of(key)].push((key, p.rec as u32, (p.start + pos) as u32, kmer.flanks::<F>().0 as u32));
                    }
                    *slots[i].lock().unwrap() = out;
                }
            });
        }
    });
    slots.into_iter().map(|m| m.into_inner().unwrap()).collect()
}

/// `default_build` on `threads` threads, with the same result: every key's positions are stored in
/// the order the FASTA lists them, as the one-thread build stores them. The k-mers of each batch of
/// records are extracted in parallel and split by key range; each thread then counts (keys pass)
/// or places (map pass) the k-mers of one key range, in sequence order, into its own part of the
/// key table and of the values -- disjoint, as a key range's cells and values are contiguous.
pub fn default_build_par<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(path: impl AsRef<Path>, max_range_size: usize, threads: usize) ->
        Result<(Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>, HashMap<String, usize>, Vec<String>), std::io::Error> {
    type Keys<const C: usize, const B: u64> = FMKeys<C, B>;
    let threads = threads.max(1);
    // Key ranges: whole control blocks, so a range's key cells and values are contiguous.
    let n_keys = 1u64 << (2 * C);
    let per = ((n_keys / threads as u64) / CELLS_PER_BODY).max(1) * CELLS_PER_BODY;
    let ranges = ((n_keys + per - 1) / per) as usize;
    let range_of = move |key: u64| ((key / per) as usize).min(ranges - 1);
    let bounds: Vec<u64> = (0..=ranges).map(|r| (r as u64 * per).min(n_keys)).collect();

    // Keys pass.
    eprintln!("Build keys ({threads} threads)");
    let mut keys = Keys::<C, CELLS_PER_BODY>::new();
    {
        let cell_bounds: Vec<usize> = bounds.iter().map(|&k| if k >= n_keys { keys.data.len() - 4 } else { Keys::<C, CELLS_PER_BODY>::kmer_to_ctrl_block_index(k) }).collect();
        for_each_batch(&path, |batch| {
            let parts = extract_batch::<K, C, F, S, L>(&batch, threads, ranges, &range_of);
            let mut rest: &mut [crate::keys::KCell] = &mut keys.data[..];
            let mut slices = Vec::with_capacity(ranges);
            let mut taken = 0usize;
            for r in 0..ranges {
                let (a, b) = std::mem::take(&mut rest).split_at_mut(cell_bounds[r + 1] - taken);
                slices.push((cell_bounds[r], a));
                taken = cell_bounds[r + 1];
                rest = b;
            }
            let parts = &parts;
            std::thread::scope(|scope| {
                for (r, (base, cells)) in slices.into_iter().enumerate() {
                    scope.spawn(move || {
                        for piece in parts {
                            for &(key, _, _, _) in &piece[r] {
                                cells[Keys::<C, CELLS_PER_BODY>::kmer_to_index(key) - base].increment();
                            }
                        }
                    });
                }
            });
        })?;
    }
    eprintln!("Keys build {} {}", HEADER_THRESHOLD, max_range_size);
    keys.build::<HEADER_THRESHOLD>(max_range_size);

    // Map pass.
    eprintln!("Build map ({threads} threads)");
    let mut flexmap = Flexmap::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::new(keys);
    let mut reference2id: HashMap<String, usize> = HashMap::new();
    let mut id2reference: Vec<String> = vec!["dummy".into()];
    let (mut total_kmers, mut total_minimizers) = (0u64, 0u64);
    {
        let keys = &flexmap.keys;
        // Each key range's values: from its first block's offset to the next range's.
        let value_bounds: Vec<usize> = bounds.iter().map(|&k| if k >= n_keys { keys.get_values_size() } else {
            keys.get_control_header_value(Keys::<C, CELLS_PER_BODY>::kmer_to_ctrl_block_index(k)) as usize }).collect();
        let values = &mut flexmap.values.data;
        for_each_batch(&path, |batch| {
            let first_id = id2reference.len();
            for (header, seq) in &batch {
                let id = id2reference.len();
                if reference2id.insert(header.clone(), id).is_some() {
                    panic!("Header {:?} has been seen before!", header);
                }
                id2reference.push(header.clone());
                total_kmers += seq.len().saturating_sub(K - 1) as u64;
            }
            let parts = extract_batch::<K, C, F, S, L>(&batch, threads, ranges, &range_of);
            total_minimizers += parts.iter().flatten().map(|v| v.len() as u64).sum::<u64>();
            let mut rest: &mut [crate::values::VCell] = &mut values[..];
            let mut slices = Vec::with_capacity(ranges);
            let mut taken = 0usize;
            for r in 0..ranges {
                let (a, b) = std::mem::take(&mut rest).split_at_mut(value_bounds[r + 1] - taken);
                slices.push((value_bounds[r], a));
                taken = value_bounds[r + 1];
                rest = b;
            }
            let parts = &parts;
            std::thread::scope(|scope| {
                for (r, (base, cells)) in slices.into_iter().enumerate() {
                    scope.spawn(move || {
                        for piece in parts {
                            for &(key, rec, pos, flanks) in &piece[r] {
                                if let Some((start, end)) = keys.vrange(key) {
                                    let mut block = crate::values::range_mut::<HEADER_THRESHOLD>(cells, start - base, end - base);
                                    block.insert(VD::set((first_id + rec as usize) as u64, pos as u64), flanks);
                                }
                            }
                        }
                    });
                }
            });
        })?;
    }
    eprintln!("Minimizer compression rate: {} ({}/{})", total_minimizers as f64 / total_kmers as f64, total_minimizers, total_kmers);
    eprintln!("Number of ids: {}", id2reference.len());
    Ok((flexmap, reference2id, id2reference))
}
