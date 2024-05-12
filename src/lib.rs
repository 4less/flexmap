#![feature(exposed_provenance)]
#![feature(const_trait_impl)]
#![feature(effects)]
#[cfg(target_pointer_width = "64")] 


const GLOBAL_VERSION: u32 = 1;

pub mod keys;
pub mod values;
pub mod flexmap;
pub mod example;

#[macro_use]
extern crate savefile_derive;

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
