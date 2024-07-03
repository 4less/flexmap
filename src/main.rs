#![feature(exposed_provenance)]
use std::{cmp::min, collections::HashMap, fs::File, mem::transmute, path::{Path, PathBuf}};

use flexmap::{example::{build_flexmap, test_flexmap}, flexmap::Flexmap, keys::FMKeys};
use kmerrs::consecutive::kmer::{Kmer, KmerIter};
use bioreader::{fasta_byte_reader, fastq_byte_reader, fasta_reader, fastq_reader};
use savefile::{load, save};


fn test_simple() {
    let seq = "CATCGATCGTACGTGACTGCGTCGTCCTGCGTCGTCGTCGTGCTGCTGCTGTCGTCGTCGTCGTGCTGTCGTCGTA";
    const K: usize = 3;
    let kiter: KmerIter<3, true> = KmerIter::<K, true>::new(seq.as_bytes());
    let mut keys = FMKeys::<K, 8>::new();
    for (pos, kmer_fwd, kmer_rev) in kiter.clone() {
        let kmer = min(kmer_fwd, kmer_rev);
        keys.get_kmer_cell_mut_ref(kmer.0).increment();
    }

    let mut map = HashMap::<u64, u64>::new();
    for (pos, kmer_fwd, kmer_rev) in kiter.clone() {
        let kmer = min(kmer_fwd, kmer_rev);
        let entry = map.entry(kmer.0).or_insert_with(|| 0);
        *entry += 1;
    }

    for (pos, kmer_fwd, kmer_rev) in kiter {
        let kmer = min(kmer_fwd, kmer_rev);
        let cell = keys.get_kmer_cell(kmer.0);
        let ctrl_value = keys.get_control_head_value_from_kmer(kmer.0);
        println!("{} -> {} -> {} ({})", kmer.to_string().unwrap(), cell.0, map[&kmer.0], ctrl_value);
    }

    println!("{}", keys.get_control_head_value_from_kmer(0));
    keys.set_control_header_value_from_kmer(0, 42);
    println!("{}", keys.get_control_head_value_from_kmer(0));

    keys.build::<2>();

}



fn main() {
    // println!("Main");
    // let testvec = vec![0u16, 1u16, 2u16, 3u16];
    // let ptr_to_u64: *const u64 = unsafe { transmute((&testvec[0] as *const u16).expose_provenance()) };
    // println!("{:?}", testvec);
    // let res = unsafe { *ptr_to_u64 };
    // dbg!(res);


    // let result: u64 = (testvec[0] as u64) | (testvec[1] as u64) << 16 | (testvec[2] as u64) << 32 | (testvec[3] as u64) << 48;
    // println!("res : {}", result);

    test_simple();

    
    let flexmap = build_flexmap();
    test_flexmap(&flexmap);

    
    const GLOBAL_VERSION: u32 = 1;
    let path = PathBuf::from("result/flexmap.save");
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };
    let _ = save(&mut file, GLOBAL_VERSION, &flexmap);

    let mut file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };
    let flexmap: Flexmap<3, 8, 8, 2> = load(&mut file, GLOBAL_VERSION).expect("Loading did not work");

}