#![feature(map_try_insert)]
#![feature(exposed_provenance)]
#![feature(const_trait_impl)]
#![feature(test)]
#![feature(generic_const_exprs)]
#![feature(core_intrinsics)]
// #![feature(effects)]

use keys::KCell;
use values::VData;
#[cfg(target_pointer_width = "64")]

const GLOBAL_VERSION: u32 = 1;

pub mod keys;
pub mod values;
pub mod flexmap;
pub mod build;
pub mod example;


extern crate test;

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

pub type VD = VData<28, 34>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}


pub fn get_value(data: &Vec<u16>) -> u64 {
    (data[0] as u64) | 
    (data[0+1] as u64) << 16 | 
    (data[0+2] as u64) << 32 | 
    (data[0+3] as u64) << 48
}

pub fn get_value2(data: &Vec<KCell>) -> u64 {
    (data[0].0 as u64) | 
    (data[0+1].0 as u64) << 16 | 
    (data[0+2].0 as u64) << 32 | 
    (data[0+3].0 as u64) << 48
}


#[cfg(test)]
mod kmer_tests {
    use test::Bencher;
    use crate::{get_value, get_value2, keys::KCell};

    #[bench]
    fn get_control_header_value(b: &mut Bencher) {
        let data = vec![1u16, 2u16, 3u16, 4u16];

        b.iter(|| get_value(&data));
    }

    #[bench]
    fn get_control_header_value2(b: &mut Bencher) {
        let mut data = Vec::new();
        data.push(KCell(1));
        data.push(KCell(2));
        data.push(KCell(3));
        data.push(KCell(4));

        b.iter(|| get_value2(&data));
    }
}
