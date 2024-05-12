#![feature(exposed_provenance)]
use std::mem::transmute;

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

#[derive(Clone, Copy, Savefile)]
pub struct KCell(pub u16);

impl KCell {
    pub fn Increment(&mut self) {
        self.0 += 1;
    }

    pub fn Set(&mut self, value: u16) {
        self.0 = value;
    }
}

#[derive(Clone, Savefile)]
pub struct FMKeys<const K: usize, const F: usize, const CELLS_PER_BODY: u64> {
    data: Vec<KCell>,
}

/// One control block contains of HEAD and BODY, the HEAD contains the offset position for all 
/// Keys under this CTRL-BLOCK in VALUES.
/// CELLS PER HE
impl<const K: usize, const F: usize, const KEYS_PER_BODY: u64>
    FMKeys<K, F, KEYS_PER_BODY>
{
    const K_BITS: usize = K * 2;
    const F_BITS: usize = F * 2;
    const KEY_TO_CTRL_BLOCK_SHIFT: u64 = KEYS_PER_BODY.ilog2() as u64; // BITSHIFT to get the ctrl block number for key
    const KEY_BLOCK_MASK: u64 = KEYS_PER_BODY - 1;
    const CELLS_PER_HEAD: u64 = 4; // Given this implementation, u16 is fixed as cell type. A head is a u64 so takes 4 cells

    const fn table_size() -> u64 {
        let number_of_keys = usize::pow(2, (K*2) as u32) as u64;
        (number_of_keys as u64) + (Self::kmer_to_ctrl_block(number_of_keys) * Self::CELLS_PER_HEAD) + Self::CELLS_PER_HEAD
    }

    const fn kmer_to_ctrl_block(canonical_kmer: u64) -> u64 {
        canonical_kmer >> Self::KEY_TO_CTRL_BLOCK_SHIFT
    }

    const fn kmer_to_ctrl_block_index(canonical_kmer: u64) -> usize {
        (Self::kmer_to_ctrl_block(canonical_kmer) * (Self::CELLS_PER_HEAD + KEYS_PER_BODY)) as usize
    }

    const fn kmer_to_index(canonical_kmer: u64) -> usize {
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
        self.data[Self::kmer_to_index(canonical_kmer)].Set(value);
    }

    pub fn kmer_to_value_range(&self, canonical_kmer: u64) -> (usize, usize) {
        let (block_index, key_index) = Self::kmer_to_indexes(canonical_kmer);
        let ctrl_block_value = self.get_control_header_value(block_index as usize);
        let value_start = ctrl_block_value as usize + self.data[key_index as usize].0 as usize;
        let value_end: usize = if (canonical_kmer & Self::KEY_BLOCK_MASK) != Self::KEY_BLOCK_MASK {
            ctrl_block_value as usize + self.data[key_index as usize + 1].0 as usize
        } else {
            self.get_control_header_value((block_index + Self::CELLS_PER_HEAD + KEYS_PER_BODY) as usize) as usize
        };
        assert!(value_start <= value_end);
        (value_start, value_end)
    }

    fn get_control_header_value(&self, index: usize) -> u64 {
        let ptr_to_u64: *const u64 = unsafe { transmute((&self.data[index].0 as *const u16).expose_provenance()) };
        unsafe { *ptr_to_u64 }
    }

    pub fn set_control_header_value(&mut self, index: usize, value: u64) {
        let mut ptr_to_u64: *mut u64 = unsafe { transmute((&self.data[index].0 as *const u16).expose_provenance()) };
        unsafe { *ptr_to_u64 = value };
    }

    pub fn set_control_header_value_from_kmer(&mut self, canonical_kmer: u64, value: u64) {
        self.set_control_header_value(Self::kmer_to_ctrl_block_index(canonical_kmer), value);
    }

    pub fn get_control_head_value_from_kmer(&self, canonical_kmer: u64) -> u64 {
        self.get_control_header_value(Self::kmer_to_ctrl_block_index(canonical_kmer))
    }

    pub fn new() -> FMKeys<K, F, KEYS_PER_BODY> {
        FMKeys {
            data: vec![KCell(0); Self::table_size().try_into().unwrap()],
        }
    }

    pub const fn calc_header_size(size: usize) -> usize {
        (size as u64 + 1 >> 1) as usize
    }

    pub fn get_values_size(&self) -> usize {
        let block_index = self.data.len() - Self::CELLS_PER_HEAD as usize;
        self.get_control_header_value(block_index) as usize
    }


    pub fn build<const HEADER_THRESHOLD: usize>(&mut self) {
        let mut value_index = 0;

        let mut block_index = usize::MAX;
        let mut block_value = 0;
        let mut running_vindex = 0;
        let mut block_vindex = 0;
        for ckmer in 0..u64::pow(2, K as u32*2) {
            let kmer = Kmer::<K>(ckmer);
            let ckmer_block_index = Self::kmer_to_ctrl_block_index(ckmer);
            if block_index != ckmer_block_index {
                block_index = ckmer_block_index;
                running_vindex += block_vindex;
                self.set_control_header_value(block_index, running_vindex);
                // println!("Block: {}: {}", block_index, running_vindex);
                block_vindex = 0;
            }
            let ckmer_count = self.get_kmer_cell(ckmer).0 as u64;
            // print!("Kmer {} {} ({}): {} -> ", ckmer, kmer.to_string().unwrap(), kmer.is_smallest_rc(), ckmer_count);
            self.set_kmer_cell(ckmer, block_vindex as u16);
            // println!("{} -> {}", block_vindex, running_vindex + block_vindex);
            block_vindex += ckmer_count + (((ckmer_count > HEADER_THRESHOLD as u64) as u64) * (Self::calc_header_size(ckmer_count as usize) as u64));
        }
        block_index = self.data.len() - Self::CELLS_PER_HEAD as usize;
        running_vindex += block_vindex;
        self.set_control_header_value(block_index, running_vindex);
        // println!("Last Block: {}: {}", block_index, running_vindex);
        

        // println!("-------------------------------");
        // for ckmer in 0..u64::pow(2, K as u32*2) {
        //     let kmer = Kmer::<K>(ckmer);
        //     let (start, end) = self.kmer_to_value_range(ckmer);
        //     println!("Range {}: {}-{}", kmer.to_string().unwrap(), start, end);
        // }
        // println!("-------------------------------");
    }


}  

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use kmerrs::consecutive::kmer::KmerIter;

    use super::*;

    #[test]
    fn test_kmer_to_indexes_1() {
        assert_eq!(FMKeys::<15, 8, 8>::kmer_to_indexes(7), (0, 1*4 + 7));
        assert_eq!(FMKeys::<15, 8, 8>::kmer_to_indexes(8), (12, 2*4 + 8));
    }

    #[test]
    fn test_table_size_1() {
        assert_eq!(FMKeys::<15, 8, 8>::table_size(), 1610612740);
        assert_eq!(FMKeys::<4, 8, 8>::table_size(), 256 + (32 * 4) + 4);
    }

    #[test]
    fn whole_table_1() {
        let seq = "CATCGATCGTACGTGACTGCGTCGTCCTGCGTCGTCGTCGTGCTGCTGCTGTCGTCGTCGTCGTGCTGTCGTCGTA";
        const K: usize = 10;
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
            assert_eq!(cell.0 as u64, map[&kmer.0]);
            let ctrl_value = keys.get_control_head_value_from_kmer(kmer.0); 
            assert_eq!(ctrl_value, 0)
        }

    }

    #[test]
    fn test_ctrl_head() {
        const K: usize = 10;
        let mut keys = FMKeys::<K, 8, 8>::new();
        keys.set_control_header_value(0, 42);
        assert_eq!(keys.get_control_header_value(0), 42);
    }
}
