
use std::{array, borrow::Borrow, cell::Cell, collections::HashMap, default, error::Error, fs::{self, File}, hash::{BuildHasher, Hash}, intrinsics::size_of, io::{Read, Write}, mem::{self, transmute}, num::Wrapping, process::exit};
use bincode::{Decode, Encode};
use fxhash::FxBuildHasher;
use bioreader::utils::time_noerr;
use kmerrs::consecutive::kmer::Kmer;

#[derive(Debug)]

// ctrl_block_keys_index
//  │                       Values array
//  │   │          │       │             │
//  │   ├──────────┤       ├─────────────┤
//  └──►│CTRL-HEAD ├──────►│             │<- ctrl_block_values_begin
//      ├──────────┤       │             │
//      │Key1      │       │             │
//      ├──────────┤       │             │
//      │Key2      │       │             │
//      ├──────────┤       │             │
//      │Key3      │       │             │
//      ├──────────┤       │             │    m_keys_per_ctrl_block
//      │Key4      │       │             │
//      ├──────────┤       │             │    #keys in this range
//      │Key5      │       │             │
//      ├──────────┤       │             │
//      │Key6      │       │             │
//      ├──────────┤       │             │
//      │Key7      │       │             │
//      ├──────────┤       │             │
//      │Key8      │       │             │
//      ├──────────┤       │             │
//      │CTRL-Block├───┐   │             │
//      ├──────────┤   │   │             │
//      │          │   │   ├─────────────┤
//      │          │   └──►│             │<- ctrl_block_values_end
//      ├──────────┤       │             │
//      │          │

#[derive(Clone, Copy, Encode, Decode, ser_raw::Serialize, bytemuck::Pod, bytemuck::Zeroable)]
#[repr(C)]
pub struct KCell(pub u16);

impl KCell {
    pub fn increment(&mut self) {
        self.0 += 1;
    }

    pub fn set(&mut self, value: u16) {
        self.0 = value;
    }
}

pub const fn table_size<const C: usize, const CELLS_PER_BODY: u64>() -> usize {
    let number_of_keys: u64 = usize::pow(2, (C*2) as u32) as u64;
    let KEY_TO_CTRL_BLOCK_SHIFT: u64 = CELLS_PER_BODY.ilog2() as u64;
    let CELLS_PER_HEAD = 4;
    ((number_of_keys) + ((number_of_keys >> KEY_TO_CTRL_BLOCK_SHIFT) * CELLS_PER_HEAD) + CELLS_PER_HEAD) as usize
}




#[derive(Clone, Encode, Decode, ser_raw::Serialize)]
#[repr(C)]
pub struct FMKeys<const C: usize, const CELLS_PER_BODY: u64> { //where [(); table_size::<C,CELLS_PER_BODY>()]: 
    pub data: Vec<KCell>,
    // data: [KCell; ],
}

/// One control block contains of HEAD and BODY, the HEAD contains the offset position for all 
/// Keys under this CTRL-BLOCK in VALUES.
/// CELLS PER HE
impl<const C: usize, const CELLS_PER_BODY: u64>
    FMKeys<C, CELLS_PER_BODY> //where [(); table_size::<C,CELLS_PER_BODY>()]:
{
    const K_BITS: usize = C * 2;
    // const F_BITS: usize = F * 2;
    const KEY_TO_CTRL_BLOCK_SHIFT: u64 = CELLS_PER_BODY.ilog2() as u64; // BITSHIFT to get the ctrl block number for key
    const KEY_BLOCK_MASK: u64 = CELLS_PER_BODY - 1;
    const CELLS_PER_HEAD: u64 = 4; // Given this implementation, u16 is fixed as cell type. A head is a u64 so takes 4 cells
    const MAX_BLOCK_VALUESSIZE: usize = 2usize.pow(size_of::<KCell>() as u32 * 8);
    const MAX_KEY_VALUESSIZE: usize = Self::MAX_BLOCK_VALUESSIZE / CELLS_PER_BODY as usize;

    pub const fn table_size() -> u64 { 
        let number_of_keys = usize::pow(2, (C*2) as u32) as u64;
        (number_of_keys as u64) + (Self::kmer_to_ctrl_block(number_of_keys) * Self::CELLS_PER_HEAD) + Self::CELLS_PER_HEAD
    }

    const fn kmer_to_ctrl_block(canonical_kmer: u64) -> u64 {
        canonical_kmer >> Self::KEY_TO_CTRL_BLOCK_SHIFT
    }

    pub const fn kmer_to_ctrl_block_index(canonical_kmer: u64) -> usize {
        return (Self::kmer_to_ctrl_block(canonical_kmer) * (Self::CELLS_PER_HEAD + CELLS_PER_BODY)) as usize;
    }

    pub fn kmer_to_index(canonical_kmer: u64) -> usize {
        let value = (Self::kmer_to_ctrl_block_index(canonical_kmer) as u64) + Self::CELLS_PER_HEAD + (canonical_kmer & Self::KEY_BLOCK_MASK);
        value as usize
    }

    const fn kmer_to_indexes(canonical_kmer: u64) -> (u64, u64) {
        let block_index = Self::kmer_to_ctrl_block_index(canonical_kmer) as u64;
        (block_index, block_index + Self::CELLS_PER_HEAD + (canonical_kmer & Self::KEY_BLOCK_MASK))
    }

    pub fn get_kmer_cell_mut_ref(&mut self, canonical_kmer: u64) -> &mut KCell {
        &mut self.data[Self::kmer_to_index(canonical_kmer)]
    }

    pub fn get_kmer_cell(&self, canonical_kmer: u64) -> KCell {
        self.data[Self::kmer_to_index(canonical_kmer)]
    }

    pub fn set_kmer_cell(&mut self, canonical_kmer: u64, value: u16) {
        self.data[Self::kmer_to_index(canonical_kmer)].set(value);
    }

    pub fn vrange(&self, canonical_kmer: u64) -> Option<(usize, usize)> {
        let (block_index, key_index) = Self::kmer_to_indexes(canonical_kmer);
        let ctrl_block_value = self.get_control_header_value(block_index as usize);

        let value_start = ctrl_block_value as usize + self.data[key_index as usize].0 as usize;
        let value_end: usize = if (canonical_kmer & Self::KEY_BLOCK_MASK) != Self::KEY_BLOCK_MASK {
            ctrl_block_value as usize + self.data[key_index as usize + 1].0 as usize
        } else {
            self.get_control_header_value((block_index + Self::CELLS_PER_HEAD + CELLS_PER_BODY) as usize) as usize
        };
        assert!(value_start <= value_end);

        if value_end - value_start == 0 {
            return None
        }

        Some((value_start, value_end))
    }


    pub fn get_value(data: &[u16]) -> u64 {
        (data[0] as u64) | 
        (data[1] as u64) << 16 | 
        (data[2] as u64) << 32 | 
        (data[3] as u64) << 48
    }

    // #[inline(always)]
    fn get_control_header_value(&self, index: usize) -> u64 {
        unsafe {
            *self.data.as_ptr().add(index).cast::<u64>()
        }
    }

    // #[inline(always)]
    fn get_control_header_value3(&self, index: usize) -> u64 {
        let (duration, result) = time_noerr(|| self.data[index]);
        println!("Data Access:  {:?}, {}", duration, result.0);

        let (duration, result) = time_noerr(|| unsafe {
            *self.data.as_ptr().add(index).cast::<u64>()
        });
        println!("Unsafe:       {:?}, {}", duration, result);

        unsafe {
            *self.data.as_ptr().add(index).cast::<u64>()
        }
    }

    // #[inline(always)]
    fn get_control_header_value2(&self, index: usize) -> u64 {
        (self.data[index].0 as u64) +
        (self.data[index+1].0 as u64) << 16 + 
        (self.data[index+2].0 as u64) << 32 + 
        (self.data[index+3].0 as u64) << 48
    }

    pub fn set_control_header_value(&mut self, index: usize, value: u64) {
        assert!(index + 3 < self.data.len());

        self.data[index + 0].0 = ((0x000000000000FFFF & value) >> 0*16) as u16;
        self.data[index + 1].0 = ((0x00000000FFFF0000 & value) >> 1*16) as u16;
        self.data[index + 2].0 = ((0x0000FFFF00000000 & value) >> 2*16) as u16;
        self.data[index + 3].0 = ((0xFFFF000000000000 & value) >> 3*16) as u16;
    }

    pub fn set_control_header_value_from_kmer(&mut self, canonical_kmer: u64, value: u64) {
        self.set_control_header_value(Self::kmer_to_ctrl_block_index(canonical_kmer), value);
    }

    pub fn get_control_head_value_from_kmer(&self, canonical_kmer: u64) -> u64 {
        self.get_control_header_value(Self::kmer_to_ctrl_block_index(canonical_kmer))
    }

    pub fn new() -> FMKeys<C, CELLS_PER_BODY> {
        FMKeys {
            data: vec![KCell(0); Self::table_size().try_into().unwrap()],
        }
    }

    pub fn with_capacity(capacity: usize) -> FMKeys<C, CELLS_PER_BODY> {
        FMKeys {
            data: Vec::with_capacity(capacity),
        }
    }

    pub const fn calc_header_size(size: usize) -> usize {
        (size as u64 + 1 >> 1) as usize
    }

    pub fn get_values_size(&self) -> usize {
        let block_index = self.data.len() - Self::CELLS_PER_HEAD as usize;
        self.get_control_header_value(block_index) as usize
    }

    pub fn save(&mut self, filename: &String) -> () {

        let mut f = File::create(&filename).expect("no file found");
        // Convert Vec<u16> to raw bytes
        let bytes: &[u8] = unsafe {
            // Get a raw pointer to the vector's data
            let ptr = self.data.as_ptr();

            // Calculate the length of the data in bytes
            let len = self.data.len() * mem::size_of::<u16>();

            // Create a slice of u8 from the raw pointer and length
            std::slice::from_raw_parts(ptr as *const u8, len)
        };
        
        f.write_all(bytes);
    }

    pub fn load(filename: &String) -> FMKeys<C, CELLS_PER_BODY> {
        use std::io::Read;
        let mut f = File::open(filename).expect("Failed to open file");

        let mut buffer = Vec::new();
        f.read_to_end(&mut buffer).expect("Failed to read file");

        // Cast the buffer to &[KCell] using bytemuck
        let data: &[KCell] = bytemuck::try_cast_slice(&buffer).expect("Buffer not properly aligned or sized for KCell");

        FMKeys { data: data.to_vec() }
    }

    pub fn build<const HEADER_THRESHOLD: usize>(&mut self, max_range_size: usize) {
        let mut value_index = 0;

        let mut block_index = usize::MAX;
        let mut block_value = 0;
        let mut running_vindex = 0;
        let mut block_vindex = 0;
        
        let mut skip = 0;

        let size = u64::pow(2, C as u32*2);

        let mut set_keys = 0;
        let last_block_index = 0;
        for ckmer in 0..size {
            let kmer = Kmer::<C>(ckmer);
            let ckmer_block_index = Self::kmer_to_ctrl_block_index(ckmer);
            if block_index != ckmer_block_index {
                block_index = ckmer_block_index;
                running_vindex += block_vindex;
                self.set_control_header_value(block_index, running_vindex);

                assert!(self.get_control_header_value(block_index) == running_vindex);
                block_vindex = 0;
            }
            let mut ckmer_count = self.get_kmer_cell(ckmer).0 as u64;

            if ckmer_count > 0 { set_keys += 1 };
            if ckmer_count as usize > max_range_size || ckmer_count as usize > Self::MAX_KEY_VALUESSIZE {
                // eprintln!("{} {}  (user_max_range: {} , data_structure_max_range {})", kmer.to_string().unwrap(), ckmer_count, max_range_size, Self::MAX_KEY_VALUESSIZE);
                skip += 1;
                ckmer_count = 0;
            }

            // print!("Kmer {} {} ({}): {} -> ", ckmer, kmer.to_string().unwrap(), kmer.is_smallest_rc(), ckmer_count);
            self.set_kmer_cell(ckmer, block_vindex as u16);
            // println!("{} -> {}", block_vindex, running_vindex + block_vindex);
            let key_vsize = ckmer_count + (((ckmer_count > HEADER_THRESHOLD as u64) as u64) * (Self::calc_header_size(ckmer_count as usize) as u64));
            block_vindex += key_vsize;
        }
        block_index = self.data.len() - Self::CELLS_PER_HEAD as usize;
        running_vindex += block_vindex;
        self.set_control_header_value(block_index, running_vindex);

        assert!(self.get_control_header_value(block_index) == running_vindex);

        eprintln!("Non null k-mers {}", set_keys);

        eprintln!("Skipped {}", skip);
    }


}  


#[derive(Clone, Encode, Decode)]
#[repr(C)]
pub struct KHashEntry {
    pub key: u32,
    pub range_len: u32,
    pub range_start: u64,
}

impl Default for KHashEntry {
    fn default() -> Self {
        Self { key: 0, range_len: 0, range_start: 0 }
    }
}

#[derive(Clone, Encode, Decode)]
#[repr(C)]
pub struct FMKeysHash {
    pub data: Vec<KHashEntry>,
    pub load_factor: f64,
}

impl FMKeysHash {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            data: vec![KHashEntry::default(); capacity],
            load_factor: 0.6,
        }
    }
}

impl KHashEntry {
    pub fn is_empty(&self) -> bool {
        return self.range_len == 0;
    }
}

impl FMKeysHash {
    #[inline(always)]
    pub fn hash(key: u64) -> u64 {
        let mut k = Wrapping(key);
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccd;
        k ^= k >> 33;
        k *= 0xc4ceb9fe1a85ec53;
        k ^= k >> 33;
        return k.0;
    }

    pub fn insert(&mut self, mut key: u32, mut range_start: u64, mut range_len: u32) -> Option<()> {
        let mut index = Self::hash(key as u64) as usize % self.data.len();
        let mut distance = 0;

        loop {
            if unsafe { self.data.get_unchecked(index).is_empty() } {
                // If the cell is empty, fill it and return
                let cell = unsafe { self.data.get_unchecked_mut(index) };
                cell.key = key;
                cell.range_start = range_start;
                cell.range_len = range_len;

                return Some(());
            }

            let cell_hash_distance = {
                let cell = unsafe { self.data.get_unchecked(index) };
                let cell_hash = Self::hash(cell.key as u64) as usize % self.data.len();
                let cell_hash_distance = if index > cell_hash {
                    index - cell_hash
                } else {
                    index + self.data.len() - cell_hash
                };
    
                cell_hash_distance
            };
    
            // After the block, the mutable borrow ends, and you can safely borrow `self.data` again
            if cell_hash_distance < distance {
                // Reborrow the cell mutably just for the swap operation
                let cell = unsafe { self.data.get_unchecked_mut(index) };
                std::mem::swap(&mut cell.key, &mut key);
                std::mem::swap(&mut cell.range_start, &mut range_start);
                std::mem::swap(&mut cell.range_len, &mut range_len);
                distance = cell_hash_distance;
            }
    
            // Increment the distance and move to the next cell
            distance += 1;
            index += 1;
            if index >= self.data.len() { index -= self.data.len() };

            if distance > self.data.len() {
                panic!("Insert failed due to insufficient coverage");
            }
        }
    }

    pub fn get(&self, key: u32) -> Option<(usize, usize)> {
        let mut index = Self::hash(key as u64) as usize % self.data.len();

        let mut distance = 0;
        let mut cell_hash_distance = {
            let cell = unsafe { self.data.get_unchecked(index) };
            let cell_hash = Self::hash(cell.key as u64) as usize % self.data.len();
            let cell_hash_distance = if index > cell_hash {
                index - cell_hash
            } else {
                index + self.data.len() - cell_hash
            };

            cell_hash_distance
        };

        while distance <= cell_hash_distance {
            let cell = unsafe { self.data.get_unchecked(index) };

            if key == cell.key {
                return Some((cell.range_start as usize, cell.range_len as usize));
            }
            distance += 1;
            index += 1;
            if index >= self.data.len() { index -= self.data.len() };

            cell_hash_distance = {
                let cell = unsafe { self.data.get_unchecked(index) };
                let cell_hash = Self::hash(cell.key as u64) as usize % self.data.len();
                let cell_hash_distance = if index > cell_hash {
                    index - cell_hash
                } else {
                    index + self.data.len() - cell_hash
                };
    
                cell_hash_distance
            };

            if distance > self.data.len() {
                panic!("Get failed due to insufficient coverage");
            }
        }
        None
    }
}



// pub struct HashKeys {
//     pub data: KeysHashSmall,
// }

impl FMKeysHash {
    pub fn vrange(&self, canonical_kmer: u64) -> Option<(usize, usize)> {
        match self.get(canonical_kmer as u32) {
            Some(entry) => Some((entry.0, entry.0 + entry.1)),
            None => None,
        }
    }
}


#[cfg(test)]
mod tests {

    use std::{cmp::min, collections::HashMap};

    use kmerrs::consecutive::kmer::KmerIter;

    use super::*;
    use test::Bencher;

    #[test]
    fn test_kmer_to_indexes_1() {
        assert_eq!(FMKeys::<15, 8>::kmer_to_indexes(7), (0, 1*4 + 7));
        assert_eq!(FMKeys::<15, 8>::kmer_to_indexes(8), (12, 2*4 + 8));
    }

    #[test]
    fn test_table_size_1() {
        assert_eq!(FMKeys::<15, 8>::table_size(), 1610612740);
        assert_eq!(FMKeys::<4, 8>::table_size(), 256 + (32 * 4) + 4);
    }

    #[test]
    fn whole_table_1() {
        let seq = "CATCGATCGTACGTGACTGCGTCGTCCTGCGTCGTCGTCGTGCTGCTGCTGTCGTCGTCGTCGTGCTGTCGTCGTA";
        const K: usize = 10;
        let kiter = KmerIter::<K, true>::new(seq.as_bytes());
        let mut keys = FMKeys::<K, 8>::new();
        for (pos, kmer_fwd, kmer_rev) in kiter.clone() {
            let kmer = min(kmer_fwd, kmer_rev);
            keys.get_kmer_cell_mut_ref(kmer.0).increment();
        }
    
        let mut map = HashMap::<u64, u64>::new();
        for (pos, kmer_fwd, kmer_rev) in kiter.clone() {
            let kmer = min(kmer_fwd, kmer_rev);
            let mut entry = map.entry(kmer.0).or_insert_with(|| 0);
            *entry += 1;
        }
    
    
        for (pos, kmer_fwd, kmer_rev) in kiter {
            let kmer = min(kmer_fwd, kmer_rev);
            let cell = keys.get_kmer_cell(kmer.0);
            assert_eq!(cell.0 as u64, map[&kmer.0]);
            let ctrl_value = keys.get_control_head_value_from_kmer(kmer.0); 
            assert_eq!(ctrl_value, 0)
        }

    }

    #[test]
    fn test_ctrl_head() {
        const K: usize = 10;
        let mut keys = FMKeys::<K, 8>::new();
        keys.set_control_header_value(0, 42);
        assert_eq!(keys.get_control_header_value(0), 42);
    }

    #[test]
    fn test_fm_keys_hash() {
        let capa = 100_000;
        let mut hashmap = FMKeysHash::with_capacity(capa);

        for i in 0..capa {
            let key = i as u32;
            let range_start = i as u64;
            hashmap.insert(key, range_start, 1);
            assert!(hashmap.get(key).unwrap() == (range_start as usize, 1))
        }

        for i in 0..capa {
            let key = i as u32;
            let range_start = i as u64;
            assert!(hashmap.get(key).unwrap() == (range_start as usize, 1))
        }
    }

    #[bench]
    fn bench_fm_keys_hash(b: &mut Bencher) {
        let size = 100_000;
        let load_factor = 0.6;
        let capa = (size as f64 * (1.0/load_factor)) as usize;
        let mut hashmap = FMKeysHash::with_capacity(capa);

        for i in 0..size {
            let key = i as u32;
            let range_start = i as u64;
            hashmap.insert(key, range_start, 1);
            assert!(hashmap.get(key).unwrap() == (range_start as usize, 1))
        }

        b.iter(|| {
            for i in 0..size {
                let key = i as u32;
                let range_start = i as u64;
                assert!(hashmap.get(key).unwrap() == (range_start as usize, 1))
            }
        });
    }

    #[bench]
    fn bench_std_hashmap(b: &mut Bencher) {
        let size = 100_000;
        let load_factor = 0.6;
        let capa = (size as f64 * (1.0/load_factor)) as usize;
        let mut hashmap = HashMap::<u32, (u32, u64)>::with_capacity(capa);

        for i in 0..size {
            let key = i as u32;
            let range_start = i as u64;
            hashmap.insert(key, (1,  range_start));
            assert!(hashmap.get(&key).unwrap() == &(1, range_start))
        }

        b.iter(|| {
            for i in 0..size {
                let key = i as u32;
                let range_start = i as u64;
                assert!(hashmap.get(&key).unwrap() == &(1, range_start))
            }
        });
    }
}
