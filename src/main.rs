#![feature(exposed_provenance)]
use std::{collections::HashMap, fs::File, mem::transmute, path::{Path, PathBuf}};

use flexmap::{example::{build_flexmap, test_flexmap}, flexmap::Flexmap, keys::FMKeys};
use kmerrs::consecutive::kmer::{Kmer, KmerIter};
use savefile::{load, save};


fn test_simple() {
    let seq = "CATCGATCGTACGTGACTGCGTCGTCCTGCGTCGTCGTCGTGCTGCTGCTGTCGTCGTCGTCGTGCTGTCGTCGTA";
    const K: usize = 3;
    let kiter = KmerIter::<K>::new(seq.as_bytes());
    let mut keys = FMKeys::<K, 8, 8>::new();
    for (pos, kmer) in kiter.clone() {
        keys.get_kmer_cell_mut_ref(kmer.0).Increment();
    }

    let mut map = HashMap::<u64, u64>::new();
    for (pos, kmer) in kiter.clone() {
        let mut entry = map.entry(kmer.0).or_insert_with(|| 0);
        *entry += 1;
    }

    for (pos, kmer) in kiter {
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

    dbg!(test_flexmap(&flexmap));
}