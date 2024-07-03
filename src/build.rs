use std::{cmp::min, collections::HashMap, fs::File, path::Path, sync::{Arc, Mutex}};

use kmerrs::consecutive::kmer::KmerIter;
use bioreader::{fasta_byte_reader::{self, FastaByteReader}, fasta_reader::{self, FastaReader}, fastq_byte_reader, fastq_reader, sequence::fasta_record::OwnedFastaRecord};

use crate::{flexmap::Flexmap, keys::{self, FMKeys}, values::VData};

fn find_min<'a, I>(vals: I) -> Option<&'a u32>
where
    I: Iterator<Item = &'a u32>,
{
    vals.min()
}

fn default_build<
    const K: usize,
    const C: usize,
    const F: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize
>(path: impl AsRef<Path>) -> Result<Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>, std::io::Error> {
    let keys = match default_build_keys::<K, C, CELLS_PER_BODY, HEADER_THRESHOLD>(&path) {
        Ok(keys) => keys,
        Err(_) => panic!("Keys could not be built."),
    };


    let mut reference2id = HashMap::<String, usize>::new();
    let mut id2reference = Vec::<String>::new();

    let map = match default_build_map::<K, C, F, CELLS_PER_BODY, HEADER_THRESHOLD>(path, keys) {
        Ok(map) => map,
        Err(_) => panic!("Map could not be built."),
    };
    Ok(map)
}

fn default_build_keys<
    const K: usize,
    const C: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(path: impl AsRef<Path>) -> Result<FMKeys<C, CELLS_PER_BODY>, std::io::Error> {
    let mut keys = FMKeys::<C, CELLS_PER_BODY>::new();

    let file: File = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.as_ref().display(), why),
        Ok(file) => file,
    };

    let mut record = OwnedFastaRecord::new();
    
    let buffer_size = usize::pow(2, 24);
    let mut byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size)?));
    let mut fasta_reader = FastaReader::with_capacity(buffer_size);


    let mut kmer_count = 0;
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
                let kmer = min(kmer_fwd, kmer_rev);
                keys.get_kmer_cell_mut_ref(kmer.middle::<C>().0).increment();
            }
        }
    }

    keys.build::<HEADER_THRESHOLD>();

    Ok(keys)
}


fn default_build_map<
    const K: usize,
    const C: usize,
    const F: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
>(path: impl AsRef<Path>, keys: FMKeys<C, CELLS_PER_BODY>) -> Result<Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>, std::io::Error> {
    let file: File = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.as_ref().display(), why),
        Ok(file) => file,
    };

    let mut flexmap = Flexmap::<C,F,CELLS_PER_BODY,HEADER_THRESHOLD>::new(keys);

    let mut record = OwnedFastaRecord::new();
    
    let buffer_size = usize::pow(2, 24);
    let mut byte_reader = Arc::new(Mutex::new(FastaByteReader::new(file, buffer_size)?));
    let mut fasta_reader = FastaReader::with_capacity(buffer_size);
    type VD = VData<20, 40>;

    let mut kmer_count = 0;
    while let Some(()) = fasta_reader
        .load_batch_par(&mut byte_reader)
        .expect("Batch is invalid")
    {
        while let Some(_) = fasta_reader.next(&mut record) {
            if !record.valid_extended() {
                panic!("Record is not valid {:?}", record.to_string())
            }
            let mut iter = KmerIter::<K, true>::new(record.seq());


            for (pos, kmer_fwd, kmer_rev) in iter {
                let kmer = min(kmer_fwd, kmer_rev);
                let core = kmer.middle::<C>();
                let flanks = kmer.flanks::<F>();
        
                let range = flexmap.keys.kmer_to_value_range(core.0);
                let mut vblock = flexmap.values.get_range_mut(range);
        
                vblock.insert(VD::set(1, pos as u64), flanks.0 as u32);
                println!("Insert {}\n{}", core.to_string().expect("Error"), vblock);
            }
        }
    }

    Ok(flexmap)
}
