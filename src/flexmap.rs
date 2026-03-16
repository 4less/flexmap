use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};
use std::mem::size_of;
use std::os::unix::fs::MetadataExt;
use std::process::exit;
use std::ptr;

use crate::keys::{FMKeys, FMKeysHash};
use crate::values::{FMValues, VRange};

pub type FlexmapStd = Flexmap<15, 16, 16, 2>;
pub type FMKeysStd = FMKeys<15, 16>;

pub type FlexmapSmall = Flexmap<3, 10, 16, 2>;
pub type FMKeysSmall = FMKeys<3, 16>;


pub type KeysHashSmall = HashMap<u32, (u32, u32)>;

use bincode::{Decode, Encode};
use bioreader::utils::{time, time_noerr};
use savefile::prelude::*;
use ser_raw::{Serialize, Serializer};
// use savefile_derive::Savefile;

const FLEXMAP_BIN_MAGIC: u64 = 0x464c45584d415031; // "FLEXMAP1"
const FLEXMAP_BIN_VERSION: u32 = 1;

#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
#[repr(C)]
struct FlexmapBinHeader {
    magic: u64,
    version: u32,
    _reserved: u32,
    keys_len: u64,
    values_len: u64,
}


pub trait FlexOptions {
    
}

pub trait VRangeGetter {
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange>;
}

pub trait DBBuilder {
    fn build(options: impl FlexOptions) -> Self;
}

// / Explanation Flexmap
// /
// /       8           +      8    = 16
// /       v                  v
// / FFFFFFFFKKKKKKKKKKKKKKKFFFFFFFF
// /              ^
// /              16 
// / 
// / Example:
// / 
// / ACGTAGCTAGCTCTGTCGTCGTCTACATCGTACTACTATCGATCGATCGC
// /       FFFFFFFFKKKKKKKKKKKKKKKFFFFFFFF
// / ->    K-mer: GTCGTCGTCTACATC (length K)
// / ->    F-mer: CTAGCTCT GTACTACT (length F)
// / 
// / 
// / K:                 First key (stored in FMKeys) with exact matching. Good value is 15. 
// /                    Do not work with values > 15 as the memory requirement will explode.
// / F:                 Flanking region around K to store. 
// / CELLS_PER_BODY:    CELLS_PER_BODY determines how many K-sized keys are grouped together in 
// /                    one block in FMKeys.
// / HEADER_THRESHOLD:  if number of occurences per k-mer is strictly larger than HEADER_THRESHOLD
// /                    a header is introduced containing the Flex side regions.
// /
// / 
// / Indexed by: flexmap.get(kmer: u64);
// / kmer needs to follow the 2-bit representation of nucleotides (A: 0, C: 1, G: 2, T: 3)
// /                                                                              A C G T
// / e.g. ACGT -> 00000000 00000000 00000000 00000000 00000000 00000000 00000000 00011011
// / The interface for this is provided with the crate kmerrs


#[derive(Clone, Savefile, Encode, Decode, ser_raw::Serialize)]
#[repr(C)]
pub struct Flexmap<
    const C: usize,
    const F: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
> {
    pub keys: FMKeys<C, CELLS_PER_BODY>,
    pub values: FMValues<F, HEADER_THRESHOLD>,
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    pub fn new(
        keys: FMKeys<C, CELLS_PER_BODY>,
    ) -> Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD> {
        let size = keys.get_values_size();
        Flexmap {
            keys: keys,
            values: FMValues::new(size),
        }
    }

    pub fn save(&self, filename: &String) {
        let mut file = File::create(filename).expect("Failed to create flexmap file");

        let header = FlexmapBinHeader {
            magic: FLEXMAP_BIN_MAGIC,
            version: FLEXMAP_BIN_VERSION,
            _reserved: 0,
            keys_len: self.keys.data.len() as u64,
            values_len: self.values.data.len() as u64,
        };

        file.write_all(bytemuck::bytes_of(&header))
            .expect("Failed to write flexmap header");
        file.write_all(bytemuck::cast_slice(&self.keys.data))
            .expect("Failed to write flexmap keys");
        file.write_all(bytemuck::cast_slice(&self.values.data))
            .expect("Failed to write flexmap values");
    }

    pub fn load(filename: &String) -> Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD> {
        let mut file = File::open(filename).expect("Failed to open flexmap file");

        let mut header_bytes = [0u8; size_of::<FlexmapBinHeader>()];
        file.read_exact(&mut header_bytes)
            .expect("Failed to read flexmap header");
        let header = bytemuck::pod_read_unaligned::<FlexmapBinHeader>(&header_bytes);

        if header.magic != FLEXMAP_BIN_MAGIC {
            panic!("Invalid flexmap binary: magic mismatch");
        }
        if header.version != FLEXMAP_BIN_VERSION {
            panic!("Invalid flexmap binary: unsupported version");
        }

        let mut keys_data = vec![crate::keys::KCell(0); header.keys_len as usize];
        let mut values_data = vec![crate::values::VCell(0); header.values_len as usize];

        file.read_exact(bytemuck::cast_slice_mut(&mut keys_data))
            .expect("Failed to read flexmap keys");
        file.read_exact(bytemuck::cast_slice_mut(&mut values_data))
            .expect("Failed to read flexmap values");

        let mut trailing = [0u8; 1];
        if file.read(&mut trailing).expect("Failed to read trailing byte") != 0 {
            panic!("Invalid flexmap binary: trailing bytes detected");
        }

        Flexmap {
            keys: FMKeys { data: keys_data },
            values: FMValues { data: values_data },
        }
    }
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize> VRangeGetter for
    Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange> {
        let range = self.keys.vrange(canonical_kmer)?;
        Some(self.values.get_range(range))
    }
}


#[derive(Clone, Savefile, Encode, Decode)]
#[repr(C)]
pub struct FlexmapHash<
    const C: usize,
    const F: usize,
    const HEADER_THRESHOLD: usize,
> {
    pub keys: FMKeysHash,
    pub values: FMValues<F, HEADER_THRESHOLD>,
}

impl<const C: usize, const F: usize, const HEADER_THRESHOLD: usize>
FlexmapHash<C, F, HEADER_THRESHOLD>
{
    pub fn new(
        keys: FMKeysHash,
    ) -> FlexmapHash<C, F, HEADER_THRESHOLD> {
        let size = keys.data.iter().fold(0, |acc, entry| {
            acc + entry.range_len
        });
        FlexmapHash {
            keys,
            values: FMValues::new(size as usize),
        }
    }

    // pub unsafe fn load(file: &mut File) -> Self {

    //     let size = file.metadata().unwrap().len();
    //     let mut buffer = vec![0u8; size as usize];
    
    //     file.read_exact(&mut buffer).expect("File read works");
    
    //     let ptr = buffer.as_ptr() as *const Self;
    //     let obj = ptr::read(ptr);
    
    //     obj
    // }

    // pub unsafe fn save<W>(&self, file: W) -> () where W: Write {
    //     // let size = file.metadata().unwrap().len();
    //     // let mut buffer = vec![0u8; size as usize];
    
    //     // file.read_exact(&mut buffer).expect("File read works");
    
    //     // let ptr = buffer.as_ptr() as *const Self;
    //     // let obj = ptr::read(ptr);
    
    //     ()
    // }
}


impl<const C: usize, const F: usize, const HEADER_THRESHOLD: usize> VRangeGetter for 
FlexmapHash<C, F, HEADER_THRESHOLD> {
    /// Gets the VRange for a given k-mer (represented as u64). A VRange has an optional header section and a value section. 
    /// The if there is more than HEADER_THRESHOLD items in the value section, there will be a header, otherwise not. The
    /// header contains additional information about the flanking regions of the k-mer (parameter F). Returns None if 
    /// No such key is stored in the flexmap.
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange> {
        let range = self.keys.get(canonical_kmer as u32)?;
        Some(self.values.get_range((range.0, range.0 + range.1)))
    }
}
