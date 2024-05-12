use std::{fmt::Display, slice};

use kmerrs::consecutive::kmer::Kmer;


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
}

pub struct VData<const POS_BITS: usize, const VAL_BITS: usize>();

impl<const POS_BITS: usize, const VAL_BITS: usize> VData<POS_BITS, VAL_BITS> {
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

#[derive(Clone, Savefile)]
pub struct VCell(u64);

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
        Self {
            header,
            positions,
        }
    }

    fn to_string(&self) -> String {
        let result: String = String::new();
        // match self.header {
        //     Some(header) => {
        //         assert_eq!(header.len(), self.positions.len());
        //         for idx in 0..header.len() {
        //             let _ = write!(result., "{}: {}\n", header[idx].to_string(), self.positions[idx].0);
        //         }
        //         Ok(())
        //     },
        //     None => {
        //         for idx in 0..header.len() {
        //             let _ = write!(f, "X: {}\n", self.positions[idx].0);
        //         }
        //         Ok(())
        //     }
        // }
        result
    }
}

impl<'a> Display for VRangeMut<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.header {
            Some(header) => {
                assert_eq!(header.len(), self.positions.len());
                for idx in 0..header.len() {
                    let _  = write!(f, "{}: {}\n", header[idx].to_string(), self.positions[idx].0);
                }
                Ok(())
            },
            None => {
                Ok(())
            }
        }
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
                        break
                    }
                }
            },
            None => {
                for idx in 0..self.positions.len() {
                    if self.positions[idx].empty() {
                        self.positions[idx].set(value);
                        break
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
                    let _ = write!(f, "{}: {}\n", header[idx].to_string(), self.positions[idx].0);
                }
                Ok(())
            },
            None => {

                Ok(())
            }
        }

    }
}


impl<'a> VRangeMut<'a> {
    pub fn new(header: Option<&'a mut [HeaderSeq]>, positions: &'a mut [VCell]) -> Self {
        Self {
            header,
            positions,
        }
    }
}


#[derive(Clone, Savefile)]
pub struct FMValues<const F: usize, const HEADER_THRESHOLD: usize> {
    data: Vec<VCell>,
}

impl<const F: usize, const HEADER_THRESHOLD: usize> FMValues<F, HEADER_THRESHOLD> {
    pub fn new(size: usize) -> Self {
        FMValues {
            data: vec![VCell(0); size],
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
            let header_slice = &self.data[start..header_size];
            let header = unsafe { slice::from_raw_parts(header_slice.as_ptr() as *const HeaderSeq, values_size) };
            let vr = VRange::new(Some(header), &self.data[start+header_size..end]);
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
            let header_slice = &mut self.data[start..start+header_size];
            let header: &mut [HeaderSeq] = unsafe { slice::from_raw_parts_mut(header_slice.as_mut_ptr() as *mut HeaderSeq, values_size) };
            let vr = VRangeMut::new(Some(header), &mut self.data[start+header_size..end]);
            vr
        } else {
            let vr = VRangeMut::new(None, &mut self.data[start..end]);
            vr
        }
        // let v = unsafe { slice::from_raw_parts(value.as_ptr() as *const i8, value.len()) };
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