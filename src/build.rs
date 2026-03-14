use std::{cmp::min, collections::HashMap, fs::File, io::BufRead, path::Path, sync::{Arc, Mutex}};

use kmerrs::{consecutive::kmer::KmerIter, minimizer::context_free::Minimizer, syncmer::closed_syncmer::ClosedSyncmer};
use bioreader::{fasta_byte_reader::{self, FastaByteReader}, fasta_reader::{self, FastaReader}, fastq_byte_reader, fastq_reader, sequence::fasta_record::OwnedFastaRecord};
use savefile::save;

use crate::{flexmap::{Flexmap, FlexmapHash, KeysHashSmall}, keys::{self, FMKeys, FMKeysHash}, values::VData, VD};

fn find_min<'a, I>(vals: I) -> Option<&'a u32>
where
    I: Iterator<Item = &'a u32>,
{
    vals.min()
}

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


pub fn hash_build<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const HEADER_THRESHOLD: usize
>(path: impl AsRef<Path>, max_range_size: usize) -> 
        Result<(FlexmapHash<C, F, HEADER_THRESHOLD>, HashMap<String, usize>, Vec<String>), std::io::Error> {
    eprintln!("Build keys");
    let keys = match hash_build_keys::<K, C, S, L, HEADER_THRESHOLD>(&path, max_range_size) {
        Ok(keys) => keys,
        Err(_) => panic!("Keys could not be built."),
    };

    eprintln!("Build map");
    let (map, reference2id, id2reference) = match hash_build_map::<K, C, F, S, L, HEADER_THRESHOLD>(path, keys) {
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

    let mut kmer_count = 0;
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
                let kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
                let cmer = min(cmer_fwd, cmer_rev);

                if !cs.is_minimizer(cmer.0) { continue };

                // if cmer.is_own_rc() { continue };

                // println!("{:?}", cmer.to_string());

                // let table_size = keys::FMKeys::<C, CELLS_PER_BODY>::table_size();
                // let index = keys::FMKeys::<C, CELLS_PER_BODY>::kmer_to_ctrl_block_index(cmer.0);

                keys.get_kmer_cell_mut_ref(cmer.0).increment();
            }
        }
    }

    eprintln!("Keys build {} {}", HEADER_THRESHOLD, max_range_size);
    keys.build::<HEADER_THRESHOLD>(max_range_size);

    Ok(keys)
}



fn hash_build_keys<
    const K: usize,
    const C: usize,
    const S: usize,
    const L: usize,
    const HEADER_THRESHOLD: usize,
>(path: impl AsRef<Path>, max_range_size: usize) -> Result<FMKeysHash, std::io::Error> {
    let mut keys_counter = HashMap::<u32, u32>::new();

    let file: File = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.as_ref().display(), why),
        Ok(file) => file,
    };

    let mut record = OwnedFastaRecord::new();
    
    let buffer_size = usize::pow(2, 24);
    let mut byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size)?));
    let mut fasta_reader = FastaReader::with_capacity(buffer_size);
    let mut cs = ClosedSyncmer::<C,S,L>::new();

    let mut kmer_count = 0;
    println!("read data");
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
                let kmer = if cmer_fwd < cmer_rev { kmer_fwd } else { kmer_rev };
                let cmer = min(cmer_fwd, cmer_rev);

                if !cs.is_minimizer(cmer.0) { continue };
                // if cmer.is_own_rc() { continue };

                // keys_counter.get
                let values = keys_counter.entry(cmer.0 as u32).or_insert(0);
                *values += 1;
            }
        }
    }

    let mut keys: FMKeysHash = FMKeysHash::with_capacity(keys_counter.len() * 2);
    
    let mut running_v = 0;
    eprintln!("Insert ranges {}", keys_counter.len());
    let mut inserted = 0;
    keys_counter.iter().for_each(|(&cmer, &size)| {
        keys.insert(cmer as u32, running_v as u64, size as u32);
        running_v += size;
    });
    eprintln!("{}", running_v);

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
    
    
    let kmer_count = 0u64;
    let mut own_rc_count = 0u64;

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

                // if cmer.is_own_rc() { 
                //     own_rc_count += 1;
                //     continue 
                // };

                let flanks = kmer.flanks::<F>();
        
                match flexmap.keys.vrange(cmer.0) {
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

    eprintln!("OwnRC {}", own_rc_count);
    eprintln!("Number of ids: {}", id2reference.len());

    Ok((flexmap, reference2id, id2reference))
}



fn hash_build_map<
    const K: usize,
    const C: usize,
    const F: usize,
    const S: usize,
    const L: usize,
    const HEADER_THRESHOLD: usize,
>(path: impl AsRef<Path>, keys: FMKeysHash) -> 
        Result<(FlexmapHash<C, F, HEADER_THRESHOLD>, HashMap<String, usize>, Vec<String>), std::io::Error> {

    let mut cs = ClosedSyncmer::<K,S,L>::new();
    let file: File = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.as_ref().display(), why),
        Ok(file) => file,
    };

    let mut flexmap = FlexmapHash::<C,F,HEADER_THRESHOLD>::new(keys);

    let mut record = OwnedFastaRecord::new();
    
    let buffer_size = usize::pow(2, 24);
    let mut byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size)?));
    let mut fasta_reader = FastaReader::with_capacity(buffer_size);

    let kmer_count = 0;
    let mut own_rc_count = 0;

    let mut reference_running_id = 1;
    let mut reference2id: HashMap<String, usize> = HashMap::<String, usize>::new();
    let mut id2reference = Vec::<String>::new();
    id2reference.push("dummy".into());

    let mut total_kmers = 0;
    let mut total_minimizers = 0;
    
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

                // if cmer.is_own_rc() { 
                //     own_rc_count += 1;
                //     continue 
                // };

                let flanks = kmer.flanks::<F>();
        

                match flexmap.keys.vrange(cmer.0 as u64) {
                    Some(range) => {
                        let mut vblock = flexmap.values.get_range_mut((range.0 as usize, range.1 as usize));
                        vblock.insert(VD::set(reference_id as u64, pos as u64), flanks.0 as u32);
                    },
                    None => {},
                };

                // if reference_id == 1 && pos == 578574;
            }
        }

    }
    eprintln!("Minimizer compression rate: {}", total_minimizers as f64/total_kmers as f64);

    eprintln!("OwnRC {}", own_rc_count);

    Ok((flexmap, reference2id, id2reference))
}
