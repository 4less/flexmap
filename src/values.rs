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

#[repr(transparent)]
pub struct HeaderSeq(u32);

/// 2-bit-per-base Hamming distance encoded in one `u32` (16 bases): number of bases that differ
/// between the two packed core-mers whose XOR is `xored`.
#[inline(always)]
fn coremer_bit_dist(xored: u32) -> u32 {
    ((xored & 0x5555_5555) | ((xored & 0xAAAA_AAAA) >> 1)).count_ones()
}

/// Minimum flank distance from `flex` to any header in `headers`, and how many headers attain it.
/// Returns `(u32::MAX, 0)` for an empty slice. Semantically identical to the scalar loop in
/// [`VRange::best_flex_match`]; dispatches to the widest SIMD the CPU supports at runtime.
/// Below this many headers, the AVX-512 path's fixed cost (setup + horizontal lane-combine)
/// outweighs its throughput and it loses to the auto-vectorized scalar loop; the crossover is
/// ~128 on Zen 5 (see `docs/simd-flex-match.md`). Small header blocks are the common case, so we
/// stay scalar below the threshold.
const SIMD_MIN_HEADERS: usize = 128;

pub fn min_dist_and_count(headers: &[HeaderSeq], flex: u32) -> (u32, u32) {
    // Only AVX-512 (with VPOPCNTDQ), and only on large-enough blocks, beats the auto-vectorized
    // scalar loop (~1.5x at 512+ headers). The hand-written AVX2 path relies on a pshufb popcount
    // and benchmarks *below* the compiler's scalar auto-vectorization, so it is deliberately not on
    // the dispatch path (kept public only for benchmarking).
    #[cfg(target_arch = "x86_64")]
    {
        if headers.len() >= SIMD_MIN_HEADERS
            && is_x86_feature_detected!("avx512f")
            && is_x86_feature_detected!("avx512bw")
            && is_x86_feature_detected!("avx512vpopcntdq")
        {
            return min_dist_and_count_avx512(headers, flex);
        }
    }
    min_dist_and_count_scalar(headers, flex)
}

pub fn min_dist_and_count_scalar(headers: &[HeaderSeq], flex: u32) -> (u32, u32) {
    let mut min = u32::MAX;
    let mut count = 0u32;
    for h in headers {
        let d = coremer_bit_dist(h.0 ^ flex);
        if d < min {
            min = d;
            count = 1;
        } else if d == min {
            count += 1;
        }
    }
    (min, count)
}

#[cfg(target_arch = "x86_64")]
fn raw_u32(headers: &[HeaderSeq]) -> &[u32] {
    // HeaderSeq is #[repr(transparent)] over u32, so this reinterpret is sound.
    unsafe { slice::from_raw_parts(headers.as_ptr() as *const u32, headers.len()) }
}

/// AVX2 path. Panics if AVX2 is unavailable (guarded by [`min_dist_and_count`]); exposed for
/// benchmarking.
#[cfg(target_arch = "x86_64")]
pub fn min_dist_and_count_avx2(headers: &[HeaderSeq], flex: u32) -> (u32, u32) {
    assert!(is_x86_feature_detected!("avx2"), "AVX2 not available");
    unsafe { min_dist_and_count_avx2_impl(raw_u32(headers), flex) }
}

/// AVX-512 (VPOPCNTDQ) path. Panics if unavailable; exposed for benchmarking.
#[cfg(target_arch = "x86_64")]
pub fn min_dist_and_count_avx512(headers: &[HeaderSeq], flex: u32) -> (u32, u32) {
    assert!(
        is_x86_feature_detected!("avx512f")
            && is_x86_feature_detected!("avx512bw")
            && is_x86_feature_detected!("avx512vpopcntdq"),
        "AVX-512 VPOPCNTDQ not available"
    );
    unsafe { min_dist_and_count_avx512_impl(raw_u32(headers), flex) }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn min_dist_and_count_avx2_impl(headers: &[u32], flex: u32) -> (u32, u32) {
    use core::arch::x86_64::*;
    let n = headers.len();
    let p = headers.as_ptr();
    let vflex = _mm256_set1_epi32(flex as i32);
    let m55 = _mm256_set1_epi32(0x5555_5555u32 as i32);
    let maa = _mm256_set1_epi32(0xAAAA_AAAAu32 as i32);
    // Nibble popcount LUT (duplicated across both 128-bit lanes for pshufb).
    let lut = _mm256_setr_epi8(
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3,
        4,
    );
    let low_mask = _mm256_set1_epi8(0x0f);
    let ones8 = _mm256_set1_epi8(1);
    let ones16 = _mm256_set1_epi16(1);

    // Per-32-bit-lane flank distance for 8 headers at once.
    let dist8 = |h: __m256i| -> __m256i {
        let x = _mm256_xor_si256(h, vflex);
        let a = _mm256_and_si256(x, m55);
        let b = _mm256_srli_epi32::<1>(_mm256_and_si256(x, maa));
        let or = _mm256_or_si256(a, b);
        // Byte-wise popcount via pshufb nibble LUT ...
        let lo = _mm256_and_si256(or, low_mask);
        let hi = _mm256_and_si256(_mm256_srli_epi16::<4>(or), low_mask);
        let pc = _mm256_add_epi8(_mm256_shuffle_epi8(lut, lo), _mm256_shuffle_epi8(lut, hi));
        // ... then fold the 4 bytes of each 32-bit lane into that lane's popcount.
        let s16 = _mm256_maddubs_epi16(pc, ones8);
        _mm256_madd_epi16(s16, ones16)
    };

    // Single pass: per-lane running (min, count). Distances are 0..=16, so a sentinel of 65535
    // keeps signed compares valid, and the `i > 0` gate below stops it leaking into the result.
    let one = _mm256_set1_epi32(1);
    let mut vmin = _mm256_set1_epi32(65535);
    let mut vcount = _mm256_setzero_si256();
    let mut i = 0usize;
    while i + 8 <= n {
        let h = _mm256_loadu_si256(p.add(i) as *const __m256i);
        let d = dist8(h);
        let lt = _mm256_cmpgt_epi32(vmin, d); // d < vmin (both positive)
        let eq = _mm256_cmpeq_epi32(d, vmin);
        vmin = _mm256_min_epu32(vmin, d);
        let inc = _mm256_add_epi32(vcount, _mm256_and_si256(eq, one)); // +1 where d == old min
        vcount = _mm256_blendv_epi8(inc, one, lt); // reset to 1 where d < old min
        i += 8;
    }

    let mut mins = [0u32; 8];
    let mut counts = [0u32; 8];
    _mm256_storeu_si256(mins.as_mut_ptr() as *mut __m256i, vmin);
    _mm256_storeu_si256(counts.as_mut_ptr() as *mut __m256i, vcount);
    combine_lanes(&mins, &counts, i > 0, p, i, n, flex)
}

/// Fold per-lane (min, count) plus the scalar tail `[i, n)` into a single (min, count).
#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn combine_lanes(
    mins: &[u32],
    counts: &[u32],
    vec_valid: bool,
    p: *const u32,
    i: usize,
    n: usize,
    flex: u32,
) -> (u32, u32) {
    let mut min = u32::MAX;
    let mut count = 0u32;
    if vec_valid {
        let mut vmin = u32::MAX;
        for &m in mins {
            if m < vmin {
                vmin = m;
            }
        }
        let mut vcount = 0u32;
        for k in 0..mins.len() {
            if mins[k] == vmin {
                vcount += counts[k];
            }
        }
        min = vmin;
        count = vcount;
    }
    let mut j = i;
    while j < n {
        let d = coremer_bit_dist(*p.add(j) ^ flex);
        if d < min {
            min = d;
            count = 1;
        } else if d == min {
            count += 1;
        }
        j += 1;
    }
    (min, count)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512bw,avx512vpopcntdq")]
unsafe fn min_dist_and_count_avx512_impl(headers: &[u32], flex: u32) -> (u32, u32) {
    use core::arch::x86_64::*;
    let n = headers.len();
    let p = headers.as_ptr();
    let vflex = _mm512_set1_epi32(flex as i32);
    let m55 = _mm512_set1_epi32(0x5555_5555u32 as i32);
    let maa = _mm512_set1_epi32(0xAAAA_AAAAu32 as i32);

    let dist16 = |h: __m512i| -> __m512i {
        let x = _mm512_xor_si512(h, vflex);
        let a = _mm512_and_si512(x, m55);
        let b = _mm512_srli_epi32::<1>(_mm512_and_si512(x, maa));
        _mm512_popcnt_epi32(_mm512_or_si512(a, b))
    };

    // Single pass: per-lane running (min, count) via mask ops (VPOPCNTDQ gives the popcount).
    let one = _mm512_set1_epi32(1);
    let mut vmin = _mm512_set1_epi32(-1); // u32::MAX sentinel (unsigned masks below)
    let mut vcount = _mm512_setzero_si512();
    let mut i = 0usize;
    while i + 16 <= n {
        let h = _mm512_loadu_si512(p.add(i) as *const __m512i);
        let d = dist16(h);
        let lt = _mm512_cmplt_epu32_mask(d, vmin);
        let eq = _mm512_cmpeq_epu32_mask(d, vmin);
        vmin = _mm512_min_epu32(vmin, d);
        let inc = _mm512_mask_add_epi32(vcount, eq, vcount, one); // +1 where d == old min
        vcount = _mm512_mask_mov_epi32(inc, lt, one); // reset to 1 where d < old min
        i += 16;
    }

    let mut mins = [0u32; 16];
    let mut counts = [0u32; 16];
    _mm512_storeu_si512(mins.as_mut_ptr() as *mut __m512i, vmin);
    _mm512_storeu_si512(counts.as_mut_ptr() as *mut __m512i, vcount);
    combine_lanes(&mins, &counts, i > 0, p, i, n, flex)
}

impl HeaderSeq {
    pub fn from_raw(v: u32) -> Self {
        HeaderSeq(v)
    }

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
                // SIMD-accelerated reduction over all headers (scalar/AVX2/AVX-512 at runtime).
                let (min_dist, count) = min_dist_and_count(headers, flex.0 as u32);

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

    fn xorshift(state: &mut u64) -> u64 {
        let mut x = *state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        *state = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }

    fn headers_from_seed(seed: u64, n: usize) -> Vec<HeaderSeq> {
        let mut s = seed | 1;
        (0..n).map(|_| HeaderSeq(xorshift(&mut s) as u32)).collect()
    }

    /// The scalar reference and every available SIMD path must return identical (min, count)
    /// across sizes that exercise full vectors, tails, empties, and forced ties.
    #[test]
    fn min_dist_and_count_simd_matches_scalar() {
        // Sizes around the AVX2 (8) and AVX-512 (16) widths, plus tails and a big block.
        let sizes = [0usize, 1, 3, 7, 8, 9, 15, 16, 17, 31, 33, 100, 1000, 4096];
        for &n in &sizes {
            for seed in 0..8u64 {
                let headers = headers_from_seed(seed.wrapping_mul(0x9E37) + 1, n);
                let flex = xorshift(&mut (seed + 12345)) as u32;
                let reference = min_dist_and_count_scalar(&headers, flex);

                #[cfg(target_arch = "x86_64")]
                {
                    if is_x86_feature_detected!("avx2") {
                        assert_eq!(reference, min_dist_and_count_avx2(&headers, flex),
                            "AVX2 mismatch at n={n}, seed={seed}");
                    }
                    if is_x86_feature_detected!("avx512f")
                        && is_x86_feature_detected!("avx512bw")
                        && is_x86_feature_detected!("avx512vpopcntdq")
                    {
                        assert_eq!(reference, min_dist_and_count_avx512(&headers, flex),
                            "AVX-512 mismatch at n={n}, seed={seed}");
                    }
                }
                // The runtime-dispatched entry point must also agree.
                assert_eq!(reference, min_dist_and_count(&headers, flex),
                    "dispatch mismatch at n={n}, seed={seed}");
            }
        }
    }

    /// Force ties (all headers equal to the query) and near-ties to check the count path.
    #[test]
    fn min_dist_and_count_tie_handling() {
        for &n in &[0usize, 1, 8, 16, 17, 50] {
            let flex = 0xDEAD_BEEFu32;
            let headers: Vec<HeaderSeq> = (0..n).map(|_| HeaderSeq(flex)).collect();
            let expected = if n == 0 { (u32::MAX, 0) } else { (0u32, n as u32) };
            assert_eq!(min_dist_and_count(&headers, flex), expected, "all-equal n={n}");
            assert_eq!(min_dist_and_count_scalar(&headers, flex), expected, "scalar all-equal n={n}");
        }
    }
}
