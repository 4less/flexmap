use std::cmp::Ordering::{Equal, Greater, Less};
use std::fs::File;
use std::io::{Read, Write};
use std::iter::zip;
use std::{cmp::Ordering, fmt::Display, slice};

use bincode::{Decode, Encode};
use bytemuck::{Pod, Zeroable};
use kmerrs::consecutive::kmer::Kmer;

use crate::VD;

/// Values holds the sequence positions a kmer occurs in
/// Each key in the keys points to a region in values
/// where all positions of that particular kmer are stored
/// Additionally, every value block above a certain size has a header
/// and the header contains sequences from the flanking region of each k-mer

pub struct HeaderSeq(u32);

impl HeaderSeq {
    pub fn to_string(&self) -> String {
        Kmer::<16>(self.0 as u64).to_string().expect("String")
    }

    pub fn set(&mut self, flank: u32) {
        self.0 = flank;
    }

    pub fn get(&self) -> u32 {
        self.0
    }

    pub fn dist(&self, flex: u32) -> u32 {
        let a = (self.0 ^ flex) & 0x55555555;
        let b = ((self.0 ^ flex) & 0xAAAAAAAA) >> 1;

        // if flex != self.0 {
        //     let self_str = (Kmer::<16> {0: flex as u64}).to_string().unwrap();
        //     let flex_str = (Kmer::<16> {0: self.0 as u64}).to_string().unwrap();

        //     let dist = zip(self_str.chars(), flex_str.chars()).into_iter()
        //         .map(|(a,b)| { (a != b) as u32 })
        //         .sum::<u32>();

        //     println!("Dist: {} ({}) between {} {}",
        //         (a | b).count_ones(),
        //         dist,
        //         self_str,
        //         flex_str);

        //     println!("self: {:#032b}", self.0);
        //     println!("flex: {:#032b}", flex);

        //     println!("a:    {:#032b}", a);
        //     println!("b:    {:#032b}", b);
        // }

        (a | b).count_ones()
    }
}

pub struct VData<const VAL_BITS: usize, const POS_BITS: usize>();

impl<const VAL_BITS: usize, const POS_BITS: usize> VData<VAL_BITS, POS_BITS> {
    const POS_MASK: u64 = (1 << POS_BITS) - 1;
    const VAL_MASK: u64 = (1 << VAL_BITS) - 1;

    pub const fn get(data: u64) -> (u64, u64) {
        let val = (data >> POS_BITS) & Self::VAL_MASK;
        let pos = data & Self::POS_MASK;
        (val, pos)
    }

    pub const fn set(val: u64, pos: u64) -> u64 {
        let mut res = 0;
        res |= (val & Self::VAL_MASK) << POS_BITS;
        res |= pos & Self::POS_MASK;
        res
    }
}

#[derive(Clone, Copy, Savefile, ser_raw::Serialize, Encode, Decode, Zeroable, Pod)]
#[repr(C)]
pub struct VCell(pub u64);

impl VCell {
    const MASK: u64 = (1 << 60) - 1;

    pub fn set_raw(&mut self, value: u64) {
        self.0 = value;
    }

    pub fn set(&mut self, value: u64) {
        self.0 |= value & Self::MASK;
    }

    pub fn get(&self) -> u64 {
        self.0
    }

    pub fn empty(&self) -> bool {
        self.0 == 0
    }
}

#[derive(Clone)]
pub struct VRange<'a> {
    pub header: Option<&'a [HeaderSeq]>,
    pub positions: &'a [VCell],
}

pub struct VRangeMut<'a> {
    pub header: Option<&'a mut [HeaderSeq]>,
    pub positions: &'a mut [VCell],
}

impl<'a> VRange<'a> {
    pub fn new(header: Option<&'a [HeaderSeq]>, positions: &'a [VCell]) -> Self {
        Self { header, positions }
    }
}

impl<'a> Display for VRangeMut<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.header {
            Some(header) => {
                assert_eq!(header.len(), self.positions.len());
                for idx in 0..header.len() {
                    let _ = write!(
                        f,
                        "{}: {}\n",
                        header[idx].to_string(),
                        self.positions[idx].0
                    );
                }
                Ok(())
            }
            None => {
                for idx in 0..self.positions.len() {
                    let _ = write!(f, ".. {}\n", self.positions[idx].0);
                }
                Ok(())
            }
        }
    }
}

impl<'a> PartialEq for VRange<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.positions.len() == other.positions.len()
    }
}

impl<'a> Eq for VRange<'a> {}

impl<'a> PartialOrd for VRange<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.positions.len().partial_cmp(&other.positions.len())
    }

    fn lt(&self, other: &Self) -> bool {
        std::matches!(self.partial_cmp(other), Some(Less))
    }

    fn le(&self, other: &Self) -> bool {
        std::matches!(self.partial_cmp(other), Some(Less | Equal))
    }

    fn gt(&self, other: &Self) -> bool {
        std::matches!(self.partial_cmp(other), Some(Greater))
    }

    fn ge(&self, other: &Self) -> bool {
        std::matches!(self.partial_cmp(other), Some(Greater | Equal))
    }
}

impl<'a> Ord for VRange<'a> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.positions.len().cmp(&other.positions.len())
    }
}

impl<'a> VRange<'a> {
    pub fn to_verbose_string(&self) -> String { //<const V: usize, const P: usize>
        let mut str = String::new();
        match &self.header {
            Some(header) => {
                assert_eq!(header.len(), self.positions.len());
                for idx in 0..header.len() {
                    let (val, pos) = VD::get(self.positions[idx].0);
                    str.push_str(&format!("{}: {} {}\n", header[idx].to_string(), val, pos));
                }
                return str;
            }
            None => {
                for idx in 0..self.positions.len() {
                    let (val, pos) = VD::get(self.positions[idx].0);
                    str.push_str(&format!(".............. : {} {}\n", val, pos));
                }
                return str;
            }
        }
    }

    pub fn best_flex_match<const F: usize, L>(&self, flex: &Kmer<F>, mut lambda: L)
    where
        L: FnMut(u64, u64, Option<(u32, u32)>) -> (), // Put in struct: rpos, rval, Option(distance, count)
    {

        match self.header {
            Some(headers) => {
                let mut count = 0;
                let mut min_dist = u32::MAX;
                for header in headers {
                    let dist = header.dist(flex.0 as u32);
                    if dist < min_dist {
                        min_dist = dist;
                        count = 0;
                    }
                    if dist == min_dist {
                        count += 1
                    };
                }
                
                // eprintln!("Header------");
                for (index, header) in headers.iter().enumerate() {
                    let dist = header.dist(flex.0 as u32);
                    if dist == min_dist {
                        let (value, rpos) = VD::get(self.positions[index].0);

                        lambda(rpos, value, Some((dist, count)));
                    }
                }
            }
            None => {
                for cell in self.positions {
                    // self.seeds.push((*pos, cell.clone()));
                    let (value, rpos) = VD::get(cell.0);
                    lambda(rpos, value, None);
                }
            }
        };
    }


    pub fn all_matches<L>(&self, mut lambda: L)
    where
        L: FnMut(u64, u64) -> (), // Put in struct: rpos, rval, Option(distance, count)
    {
        for cell in self.positions {
            // self.seeds.push((*pos, cell.clone()));
            let (value, rpos) = VD::get(cell.0);
            lambda(rpos, value);
        }
    }

    pub fn len(&self) -> usize {
        self.positions.len()
    }
}

impl<'a> VRangeMut<'a> {
    pub fn insert(&mut self, value: u64, flanks: u32) -> () {
        match &mut self.header {
            Some(header) => {
                assert_eq!(header.len(), self.positions.len());
                for idx in 0..self.positions.len() {
                    if self.positions[idx].empty() {
                        self.positions[idx].set(value);
                        header[idx].set(flanks);
                        break;
                    }
                }
            }
            None => {
                for idx in 0..self.positions.len() {
                    if self.positions[idx].empty() {
                        self.positions[idx].set(value);
                        break;
                    }
                }
            }
        }
    }
}

impl<'a> Display for VRange<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.header {
            Some(header) => {
                assert_eq!(header.len(), self.positions.len());
                for idx in 0..header.len() {
                    let _ = write!(
                        f,
                        "{}: {}\n",
                        header[idx].to_string(),
                        self.positions[idx].0
                    );
                }
                Ok(())
            }
            None => Ok(()),
        }
    }
}

impl<'a> VRangeMut<'a> {
    pub fn new(header: Option<&'a mut [HeaderSeq]>, positions: &'a mut [VCell]) -> Self {
        Self { header, positions }
    }

    fn to_verbose_string<const V: usize, const P: usize>(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> String {
        match &self.header {
            Some(header) => {
                assert_eq!(header.len(), self.positions.len());
                for idx in 0..header.len() {
                    let (val, pos) = VD::get(self.positions[idx].0);
                    let _ = write!(f, "{}: {} {}\n", header[idx].to_string(), val, pos);
                }
                let mut string = String::new();
                f.write_str(&string);
                return string;
            }
            None => {
                for idx in 0..self.positions.len() {
                    let (val, pos) = VD::get(self.positions[idx].0);
                    let _ = write!(f, "............. {} {}\n", val, pos);
                }
                let mut string = String::new();
                f.write_str(&string);
                return string;
            }
        }
    }
}

#[derive(Clone, Savefile, ser_raw::Serialize, Encode, Decode)]
#[repr(C)]
pub struct FMValues<const F: usize, const HEADER_THRESHOLD: usize> {
    pub data: Vec<VCell>,
}

impl<const F: usize, const HEADER_THRESHOLD: usize> FMValues<F, HEADER_THRESHOLD> {
    pub fn new(size: usize) -> Self {
        FMValues {
            data: vec![VCell(0); size],
        }
    }

    pub fn with_capacity(size: usize) -> Self {
        FMValues {
            data: Vec::with_capacity(size),
        }
    }

    pub fn get_header_size(vblock_size: usize) -> usize {
        (vblock_size + 2) / 3
    }

    pub fn get_range(&self, range: (usize, usize)) -> VRange {
        let (start, end) = range;
        let size = end - start;

        if size > HEADER_THRESHOLD {
            let header_size = Self::get_header_size(size);
            let values_size = size - header_size;
            let header_slice = &self.data[start..start + header_size];
            let header = unsafe {
                slice::from_raw_parts(header_slice.as_ptr() as *const HeaderSeq, values_size)
            };
            let vr = VRange::new(Some(header), &self.data[start + header_size..end]);
            vr
        } else {
            let vr = VRange::new(None, &self.data[start..end]);
            vr
        }
        // let v = unsafe { slice::from_raw_parts(value.as_ptr() as *const i8, value.len()) };
    }

    pub fn get_range_mut(&mut self, range: (usize, usize)) -> VRangeMut {
        let (start, end) = range;
        let size: usize = end - start;

        if size > HEADER_THRESHOLD {
            let header_size = Self::get_header_size(size);
            let values_size = size - header_size;
            let header_slice = &mut self.data[start..start + header_size];
            let header: &mut [HeaderSeq] = unsafe {
                slice::from_raw_parts_mut(header_slice.as_mut_ptr() as *mut HeaderSeq, values_size)
            };
            let vr = VRangeMut::new(Some(header), &mut self.data[start + header_size..end]);
            vr
        } else {
            // println!("{} {} -> {}, HT {} HAS HEADER {} SLICESIZE {} len data {}", start, end, size, HEADER_THRESHOLD, size > HEADER_THRESHOLD, end - start, self.data.len());
            let vr = VRangeMut::new(None, &mut self.data[start..end]);

            // let slice = &mut self.data[start..end];
            vr
        }
        // let v = unsafe { slice::from_raw_parts(value.as_ptr() as *const i8, value.len()) };
    }

    pub fn save(&self, filename: &String) -> () {
        let mut f = File::create(filename).expect("no file found");
        let bytes: &[u8] = bytemuck::cast_slice(&self.data);
        f.write_all(bytes).expect("write failed");
    }

    pub fn load(filename: &String) -> FMValues<F, HEADER_THRESHOLD> {
        let mut f = File::open(filename).expect("no file found");
        let mut bytes = Vec::<u8>::new();
        f.read_to_end(&mut bytes).expect("read failed");
        let cell_size = std::mem::size_of::<VCell>();
        assert!(bytes.len() % cell_size == 0, "invalid values file size");
        let len = bytes.len() / cell_size;
        let mut data = vec![VCell(0); len];
        bytemuck::cast_slice_mut::<VCell, u8>(&mut data).copy_from_slice(&bytes);
        FMValues { data }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use kmerrs::consecutive::kmer::KmerIter;

    use super::*;

    #[test]
    fn test_kmer_to_indexes_1() {
        assert_eq!(FMValues::<16, 2>::get_header_size(5), 2); // 2+3
        assert_eq!(FMValues::<16, 2>::get_header_size(6), 2); // 2+4
        assert_eq!(FMValues::<16, 2>::get_header_size(8), 3); // 3+5
        assert_eq!(FMValues::<16, 2>::get_header_size(9), 3); // 3+6
        assert_eq!(FMValues::<16, 2>::get_header_size(11), 4); // 4+7
        assert_eq!(FMValues::<16, 2>::get_header_size(12), 4); // 4+8
    }
}
