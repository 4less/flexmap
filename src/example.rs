use std::{cmp::min, collections::HashMap};

use kmerrs::consecutive::kmer::{Kmer, KmerIter};

use crate::{flexmap::{self, FMKeysSmall, Flexmap, VRangeGetter}, keys::FMKeys, values::VData, VD};


pub fn build_keys() -> FMKeys::<3, 8> {
    let seq = "CATCGATCGTACGTGACTGCGTCGTCCTGCGTCGTCGTCGTGCTGCTGCTGTCGTCGTCGTCGTGCTGTCGTCGTA";
    const K: usize = 11;
    const C: usize = 3;
    let kiter = KmerIter::<K, true>::new(seq.as_bytes());
    let mut keys = FMKeys::<C, 8>::new();
    for (_, kmer_fwd, kmer_rev) in kiter.clone() {
        let kmer = min(kmer_fwd, kmer_rev);
        keys.get_kmer_cell_mut_ref(kmer.middle::<C>().0).increment();
    }

    keys.build::<2>(1000);

    keys
}

pub fn build_flexmap() -> Flexmap<3,8,8,2> {

    const C: usize = 3;
    const F: usize = 8;
    const K: usize = C + F;

    let mut keys = build_keys();

    let mut flexmap = Flexmap::<C,F,8,2>::new(keys);

    let seq = "CATCGATCGTACGTGACTGCGTCGTCCTGCGTCGTCGTCGTGCTGCTGCTGTCGTCGTCGTCGTGCTGTCGTCGTA";

    let kiter = KmerIter::<K, true>::new(seq.as_bytes());


    for (pos, kmer_fwd, kmer_rev) in kiter.clone() {
        let kmer = min(kmer_fwd, kmer_rev);
        let core = kmer.middle::<C>();
        let flanks = kmer.flanks::<F>();

        match flexmap.keys.vrange(core.0) {
            Some(range) => {
                let mut vblock = flexmap.values.get_range_mut(range);
                vblock.insert(VD::set(1, pos as u64), flanks.0 as u32);
                println!("Insert {}\n{}", core.to_string().expect("Error"), vblock);
            },
            None => todo!(),
        }
        
    }

    flexmap
}

pub fn test_flexmap(flexmap: &Flexmap<3,8,8,2>)  {
    // let keys = build_keys();
    let kmer1 = Kmer::<3>::from_slice("ACG".as_bytes()).expect("Error").get_smallest_rc();
    
    // let range = flexmap.get(kmer1.0);
    let range = flexmap.get_vrange(kmer1.0);

    println!("Value range (optional header + stored values)");
    println!("Key: {}\n{}", kmer1.to_string().expect(""), range.expect("Value"));

}
