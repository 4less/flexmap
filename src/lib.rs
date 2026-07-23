#![feature(map_try_insert)]
#![feature(exposed_provenance)]
#![feature(const_trait_impl)]
#![feature(test)]
#![feature(generic_const_exprs)]
#![feature(core_intrinsics)]
// #![feature(effects)]

use values::VData;
#[cfg(target_pointer_width = "64")]

const GLOBAL_VERSION: u32 = 1;

pub mod keys;
pub mod values;
pub mod flexmap;
pub mod build;
pub mod example;


#[macro_use]
extern crate savefile_derive;
extern crate test;

pub type VD = VData<28, 34>;