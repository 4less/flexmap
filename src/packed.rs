//! PackedFlexmap: a Flexmap's content in less memory, with the same answers.
//!
//! A [`Flexmap`] keeps every value in a 64-bit cell (`VD`: a 28-bit reference id and a 34-bit
//! position) and, for a key of more than `HEADER_THRESHOLD` values, ceil(n/2) cells of 32-bit flank
//! cells ("header") before them, so a key of S cells holds n = S - (S + 2) / 3 values. Its key table
//! is addressed directly: a 16-bit offset for every one of the 4^C keys and a 64-bit control head
//! per block, 2.7 GB at C = 15 whatever the reference holds.
//!
//! Here:
//!  - [`PackedValues`]: a value takes the bits the reference needs, `ref_bits + pos_bits` (the
//!    largest reference id and position decide them), not 64. A key's values stay where its cells
//!    were, in units the key table already counts: cell `s` starts at bit `s * slot_bits32 / 32`,
//!    where `slot_bits32 / 32` is the least number of bits per cell that holds every key's flank
//!    cells and values (32 + W bits a value for a key with a header, which has at most 2/3 of a
//!    value per cell; W for a key without). The key table's offsets are therefore unchanged.
//!    A header key's flank cells start at the next byte (so they read as plain `u32`s), then its
//!    values; each value can be read on its own.
//!  - [`CompactKeys`]: the key table of the used keys alone (a rank structure: a presence bit per
//!    key, 448 keys to a 64-byte line whose last word counts the used keys before it, and a 24-bit
//!    start per used key). About 0.15 GB plus 3 bytes a used key, instead of 2.7 GB.
//!  - [`DirectKeys`]: the direct table as in Flexmap, for comparison and for dense indexes.
//!  - [`PackedFlexmap::for_each_prefetched`]: a lookup pipeline that fetches the key table and
//!    then the values of the lookups ahead, so that their cache misses overlap.
//!  - Files: one blob of page-aligned sections, loaded into a huge-page buffer on all threads or
//!    memory-mapped.
//!
//! The packed index is built from a Flexmap ([`PackedFlexmap::from_flexmap`]) and answers every key
//! with the same values, flank cells and order.

use std::fs::File;
use std::io::{Read, Write};
use std::os::unix::fs::FileExt;
use std::marker::PhantomData;

use kmerrs::consecutive::kmer::Kmer;

use crate::flexmap::{Flexmap, FlexmapBlob};
use crate::keys::{FMKeys, KCell};
use crate::values::VCell;
use crate::memory::{prefetch, read_parallel, zeroed_words, Backing, Section, SharedBytes};
use crate::VD;

#[cfg(not(target_endian = "little"))]
compile_error!("the packed index reads its bit fields as little-endian words");

/// Reads the `w` bits (w <= 64) at bit `bit` of `base`. At least 16 bytes from `bit / 8` on must be
/// readable (the values carry that padding).
#[inline(always)]
unsafe fn read_bits(base: *const u8, bit: u64, w: u32) -> u64 {
    let p = base.add((bit >> 3) as usize);
    let shift = (bit & 7) as u32;
    let mask = if w == 64 { u64::MAX } else { (1u64 << w) - 1 };
    if w + shift <= 64 {
        (p.cast::<u64>().read_unaligned() >> shift) & mask
    } else {
        ((p.cast::<u128>().read_unaligned() >> shift) as u64) & mask
    }
}

/// ORs the low `w` bits of `x` (the rest 0) into the zeroed bits at `bit` of `words`.
#[inline(always)]
fn put_bits(words: &mut [u64], bit: u64, x: u64, w: u32) {
    debug_assert!(w == 64 || x >> w == 0);
    let i = (bit >> 6) as usize;
    let shift = (bit & 63) as u32;
    words[i] |= x << shift;
    if shift + w > 64 {
        words[i + 1] |= x >> (64 - shift);
    }
}

/// Number of bits needed for `x` (at least 1).
fn bits_for(x: u64) -> u32 {
    (64 - x.leading_zeros()).max(1)
}

/// The widths of a packed index's values and the bits a cell of the key table's units takes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PackedLayout {
    pub ref_bits: u32,
    pub pos_bits: u32,
    pub header_threshold: u32,
    /// Bits per cell, in 1/32 bit.
    pub slot_bits32: u64,
}

impl PackedLayout {
    /// The layout for values of reference ids up to `max_ref`, positions up to `max_pos` and keys of
    /// up to `max_slots` cells.
    pub fn new(max_ref: u64, max_pos: u64, max_slots: u64, header_threshold: usize) -> Self {
        let ref_bits = bits_for(max_ref);
        let pos_bits = bits_for(max_pos);
        let w = (ref_bits + pos_bits) as u64;
        assert!(w <= 64, "a value needs {w} bits");
        let t = header_threshold as u64;
        // A key without a header: S = n values in S cells.
        let mut slot_bits32 = 32 * w;
        // A key with one: n flank cells and n values in S = n + ceil(n/2) cells, plus up to 7 bits to
        // start the flank cells at a byte. The region of S cells spans at least floor(S * bits) bits
        // wherever it starts, so the bits per cell must be at least (n (32 + W) + 7) / S.
        let mut n = t + 1;
        loop {
            let s = n + (n + 1) / 2;
            if s > max_slots {
                break;
            }
            let need = 32 * (n * (32 + w) + 7);
            slot_bits32 = slot_bits32.max((need + s - 1) / s);
            n += 1;
        }
        Self { ref_bits, pos_bits, header_threshold: header_threshold as u32, slot_bits32 }
    }

    #[inline(always)]
    pub fn entry_bits(&self) -> u32 {
        self.ref_bits + self.pos_bits
    }

    /// The first bit of cell `slot`.
    #[inline(always)]
    pub fn slot_bit(&self, slot: u64) -> u64 {
        (slot * self.slot_bits32) >> 5
    }

    /// The values a key of `slots` cells holds (as `FMValues::get_range` splits them).
    #[inline(always)]
    pub fn entries(&self, slots: u64) -> u64 {
        if slots > self.header_threshold as u64 { slots - (slots + 2) / 3 } else { slots }
    }

    /// Bits per value as stored, averaged over the cells (for reports).
    pub fn bits_per_slot(&self) -> f64 {
        self.slot_bits32 as f64 / 32.0
    }
}

/// The values of one key: its flank cells (if it has a header) and its values, read from the packed
/// bits. The counterpart of [`crate::values::VRange`].
#[derive(Clone, Copy)]
pub struct PackedRange<'a> {
    flex: &'a [u8],
    entries: *const u8,
    entry_bit: u64,
    len: usize,
    entry_bits: u32,
    pos_bits: u32,
    _bytes: PhantomData<&'a u8>,
}

// A PackedRange reads bytes that are never written while it lives (the index's sections are
// immutable once built), so it can cross threads like the slices of a VRange.
unsafe impl Send for PackedRange<'_> {}
unsafe impl Sync for PackedRange<'_> {}

/// The flank bases (2 bits each) that differ between flank cells `a` and `b`, as `HeaderSeq::dist`.
#[inline(always)]
pub fn flank_distance(a: u32, b: u32) -> u32 {
    let x = a ^ b;
    ((x & 0x5555_5555) | ((x & 0xAAAA_AAAA) >> 1)).count_ones()
}

/// Minimum flank distance from `flex` to the flank cells in `cells` (4 bytes each) and how many
/// cells have it; `(u32::MAX, 0)` for none. The same as [`crate::values::min_dist_and_count`].
#[inline]
pub fn min_dist_and_count_bytes(cells: &[u8], flex: u32) -> (u32, u32) {
    let mut min = u32::MAX;
    let mut count = 0u32;
    for c in cells.chunks_exact(4) {
        let d = flank_distance(u32::from_le_bytes([c[0], c[1], c[2], c[3]]), flex);
        if d < min {
            min = d;
            count = 1;
        } else if d == min {
            count += 1;
        }
    }
    (min, count)
}

impl<'a> PackedRange<'a> {
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline(always)]
    pub fn has_header(&self) -> bool {
        !self.flex.is_empty()
    }

    /// Value `i` as `(reference id, position)`, as `VD::get` gives a Flexmap's cell.
    #[inline(always)]
    pub fn get(&self, i: usize) -> (u64, u64) {
        assert!(i < self.len);
        let e = unsafe { read_bits(self.entries, self.entry_bit + i as u64 * self.entry_bits as u64, self.entry_bits) };
        (e >> self.pos_bits, e & ((1u64 << self.pos_bits) - 1))
    }

    /// Value `i` (< len) of a range whose values take at most 57 bits: one unaligned 8-byte read.
    #[inline(always)]
    fn get_narrow(&self, i: usize) -> (u64, u64) {
        let bit = self.entry_bit + i as u64 * self.entry_bits as u64;
        let mask = u64::MAX >> (64 - self.entry_bits);
        let e = unsafe { (self.entries.add((bit >> 3) as usize).cast::<u64>().read_unaligned() >> (bit & 7)) & mask };
        (e >> self.pos_bits, e & ((1u64 << self.pos_bits) - 1))
    }

    /// Calls `f(i, value i)` for every value, with the width test outside the loop.
    #[inline(always)]
    fn for_each_value<L: FnMut(usize, (u64, u64))>(&self, mut f: L) {
        if self.entry_bits <= 57 {
            for i in 0..self.len {
                f(i, self.get_narrow(i));
            }
        } else {
            for i in 0..self.len {
                f(i, self.get(i));
            }
        }
    }

    /// Flank cell `i`, if the key has a header.
    #[inline(always)]
    pub fn flex(&self, i: usize) -> Option<u32> {
        let c = self.flex.get(4 * i..4 * i + 4)?;
        Some(u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
    }

    /// The flank cells, 4 bytes each (empty without a header).
    pub fn flex_bytes(&self) -> &'a [u8] {
        self.flex
    }

    /// As [`crate::values::VRange::best_flex_match`]: with a header, the values whose flank cell is
    /// nearest to `flex`, with `(distance, how many)`; without one, every value.
    #[inline]
    pub fn best_flex_match<const F: usize, L>(&self, flex: &Kmer<F>, mut lambda: L)
    where
        L: FnMut(u64, u64, Option<(u32, u32)>),
    {
        if self.has_header() && self.len <= DIST_BUFFER && one_pass_seeds() {
            // Each flank's distance once, into a byte buffer (vectorized loops), then the minimum,
            // its ties and their values.
            let q = flex.0 as u32;
            let mut dist = [0u8; DIST_BUFFER];
            let dist = &mut dist[..self.len];
            for (d, c) in dist.iter_mut().zip(self.flex.chunks_exact(4)) {
                *d = flank_distance(u32::from_le_bytes([c[0], c[1], c[2], c[3]]), q) as u8;
            }
            let min_dist = *dist.iter().min().unwrap();
            let count = dist.iter().filter(|&&d| d == min_dist).count() as u32;
            let narrow = self.entry_bits <= 57;
            for (i, &d) in dist.iter().enumerate() {
                if d == min_dist {
                    let (value, rpos) = if narrow { self.get_narrow(i) } else { self.get(i) };
                    lambda(rpos, value, Some((min_dist as u32, count)));
                }
            }
        } else if self.has_header() {
            let q = flex.0 as u32;
            let (min_dist, count) = min_dist_and_count_bytes(self.flex, q);
            for (i, c) in self.flex.chunks_exact(4).enumerate() {
                let dist = flank_distance(u32::from_le_bytes([c[0], c[1], c[2], c[3]]), q);
                if dist == min_dist {
                    let (value, rpos) = if self.entry_bits <= 57 { self.get_narrow(i) } else { self.get(i) };
                    lambda(rpos, value, Some((dist, count)));
                }
            }
        } else {
            self.for_each_value(|_, (value, rpos)| lambda(rpos, value, None));
        }
    }

    /// Every value as `(position, reference id)`, as [`crate::values::VRange::all_matches`].
    #[inline]
    pub fn all_matches<L>(&self, mut lambda: L)
    where
        L: FnMut(u64, u64),
    {
        self.for_each_value(|_, (value, rpos)| lambda(rpos, value));
    }
}

/// Keys of up to this many values get the one-pass flank scan of `best_flex_match`.
const DIST_BUFFER: usize = 4096;

/// Experiment switches, read once from the environment: `FLEXMAP_ONE_PASS=0` (two-pass flank scan)
/// and `FLEXMAP_PF_LINES` (value lines prefetched per lookup, default 8).
fn one_pass_seeds() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("FLEXMAP_ONE_PASS").map(|v| v != "0").unwrap_or(true))
}

fn prefetch_lines() -> usize {
    static N: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *N.get_or_init(|| std::env::var("FLEXMAP_PF_LINES").ok().and_then(|v| v.parse().ok()).unwrap_or(8))
}

/// Bytes after the last value: a value is read 16 bytes at a time from its first byte.
const VALUE_PADDING: usize = 64;

/// The values of every key, packed ([`PackedLayout`]), addressed by the key table's cells.
pub struct PackedValues {
    pub layout: PackedLayout,
    total_slots: u64,
    data: Section,
}

impl PackedValues {
    /// The bytes the values take.
    pub fn memory_bytes(&self) -> usize {
        self.data.resident()
    }

    pub fn total_slots(&self) -> u64 {
        self.total_slots
    }

    /// The values of the key whose cells are `[start, start + slots)`, slots >= 1.
    #[inline(always)]
    pub fn range(&self, start: u64, slots: u64) -> PackedRange<'_> {
        let l = &self.layout;
        let n = l.entries(slots);
        let bit = l.slot_bit(start);
        let base = self.data.ptr();
        if slots > l.header_threshold as u64 {
            let flex_byte = ((bit + 7) >> 3) as usize;
            let flex = unsafe { std::slice::from_raw_parts(base.wrapping_add(flex_byte), 4 * n as usize) };
            PackedRange {
                flex,
                entries: unsafe { base.wrapping_add(flex_byte + 4 * n as usize) },
                entry_bit: 0,
                len: n as usize,
                entry_bits: l.entry_bits(),
                pos_bits: l.pos_bits,
                _bytes: PhantomData,
            }
        } else {
            PackedRange {
                flex: &[],
                entries: unsafe { base.wrapping_add((bit >> 3) as usize) },
                entry_bit: bit & 7,
                len: n as usize,
                entry_bits: l.entry_bits(),
                pos_bits: l.pos_bits,
                _bytes: PhantomData,
            }
        }
    }

    /// Fetches the first lines (up to `FLEXMAP_PF_LINES`, 8) of the values of cells `[start, start + slots)`.
    #[inline(always)]
    pub fn prefetch(&self, start: u64, slots: u64) {
        let max_lines = prefetch_lines();
        let first = (self.layout.slot_bit(start) >> 3) as usize & !63;
        let end = ((self.layout.slot_bit(start + slots) + 7) >> 3) as usize;
        let base = self.data.ptr();
        let mut b = first;
        let mut lines = 0;
        while b < end && lines < max_lines {
            prefetch(unsafe { base.wrapping_add(b) });
            b += 64;
            lines += 1;
        }
    }

    /// Packs a Flexmap's value cells, whose keys `each_key` lists in key order as (key, start, cells).
    fn from_cells<const T: usize>(cells: &[VCell], each_key: &dyn Fn(&mut dyn FnMut(u64, u64, u64))) -> Self {
        // Pass 1: the widths.
        let (mut max_ref, mut max_pos, mut max_slots) = (0u64, 0u64, 0u64);
        each_key(&mut |_key, start, slots| {
            max_slots = max_slots.max(slots);
            let h = if slots > T as u64 { (slots + 2) / 3 } else { 0 };
            for c in &cells[(start + h) as usize..(start + slots) as usize] {
                let (val, pos) = VD::get(c.0);
                max_ref = max_ref.max(val);
                max_pos = max_pos.max(pos);
            }
        });
        let layout = PackedLayout::new(max_ref, max_pos, max_slots.max(1), T);
        let total_slots = cells.len() as u64;
        let bytes = ((layout.slot_bit(total_slots) + 7) >> 3) as usize + VALUE_PADDING;
        let mut words = zeroed_words((bytes + 7) / 8);
        let w = layout.entry_bits();
        // Pass 2: each key's flank cells and values.
        each_key(&mut |_key, start, slots| {
            let n = layout.entries(slots);
            let bit = layout.slot_bit(start);
            let region = &cells[start as usize..(start + slots) as usize];
            let (entry_bit, values) = if slots > T as u64 {
                let h = ((slots + 2) / 3) as usize;
                let flex_byte = (bit + 7) >> 3;
                for i in 0..n as usize {
                    // Flank cell i: half i % 2 of header cell i / 2 (the low half first).
                    let flank = (region[i / 2].0 >> (32 * (i % 2))) as u32;
                    put_bits(&mut words, (flex_byte + 4 * i as u64) * 8, flank as u64, 32);
                }
                ((flex_byte + 4 * n) * 8, &region[h..h + n as usize])
            } else {
                (bit, region)
            };
            for (i, c) in values.iter().enumerate() {
                let (val, pos) = VD::get(c.0);
                put_bits(&mut words, entry_bit + i as u64 * w as u64, (val << layout.pos_bits) | pos, w);
            }
            debug_assert!(entry_bit + n * w as u64 <= layout.slot_bit(start + slots));
        });
        Self { layout, total_slots, data: Section::from_words(words, bytes) }
    }
}

impl PackedValues {
    /// The values of the keys `each_key` lists as (source first cell, cells), read from `src` and
    /// written one after the other from cell 0, in `src`'s layout. A key's cells start at a different
    /// bit than in `src`, so its flank cells and values are re-encoded, not copied.
    fn from_packed(src: &PackedValues, total_slots: u64, each_key: &dyn Fn(&mut dyn FnMut(u64, u64))) -> Self {
        let layout = src.layout;
        let bytes = ((layout.slot_bit(total_slots) + 7) >> 3) as usize + VALUE_PADDING;
        let mut words = zeroed_words((bytes + 7) / 8);
        let w = layout.entry_bits();
        let mut start = 0u64;
        each_key(&mut |src_start, slots| {
            let range = src.range(src_start, slots);
            let n = range.len() as u64;
            let bit = layout.slot_bit(start);
            let entry_bit = if range.has_header() {
                let flex_byte = (bit + 7) >> 3;
                for (i, c) in range.flex_bytes().chunks_exact(4).enumerate() {
                    let flank = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
                    put_bits(&mut words, (flex_byte + 4 * i as u64) * 8, flank as u64, 32);
                }
                (flex_byte + 4 * n) * 8
            } else {
                bit
            };
            for i in 0..n {
                let (val, pos) = range.get(i as usize);
                put_bits(&mut words, entry_bit + i * w as u64, (val << layout.pos_bits) | pos, w);
            }
            debug_assert!(entry_bit + n * w as u64 <= layout.slot_bit(start + slots));
            start += slots;
        });
        assert_eq!(start, total_slots, "cells of the listed keys");
        Self { layout, total_slots, data: Section::from_words(words, bytes) }
    }
}

/// Where a key's cells are: `locate` gives (first cell, number of cells) for a key with values.
pub trait SlotLocator: Send + Sync + Sized {
    /// Stored in the file.
    const KIND: u32;
    fn locate(&self, key: u64) -> Option<(u64, u64)>;
    /// Fetches what `locate(key)` reads first.
    fn prefetch_key(&self, key: u64);
    /// Fetches what `locate(key)` reads next, once the first fetch has arrived (nothing for a table
    /// that answers from one line).
    #[inline(always)]
    fn prefetch_key_second(&self, _key: u64) {}
    fn memory_bytes(&self) -> usize;
    fn sections(&self) -> Vec<Section>;
    fn params(&self) -> [u64; 4];
    fn from_sections(sections: &[Section], params: [u64; 4]) -> Self;
    /// Calls `f(key, first cell, cells)` for every key in `[lo, hi)` with values, in key order.
    fn for_each_in(&self, lo: u64, hi: u64, f: &mut dyn FnMut(u64, u64, u64)) {
        for key in lo..hi {
            if let Some((start, slots)) = self.locate(key) {
                f(key, start, slots);
            }
        }
    }
    /// The key table as lookups of keys `[lo, hi)` read it, and nothing else of it: `read(i, from,
    /// to)` gives bytes `[from, to)` of section `i` (of `lens[i]` bytes) as a [`Section::window`].
    /// The result must only be asked for keys in `[lo, hi)`. By default the whole table is read.
    fn read_window(
        params: [u64; 4],
        lens: &[usize],
        lo: u64,
        hi: u64,
        read: &mut dyn FnMut(usize, usize, usize) -> std::io::Result<Section>,
    ) -> std::io::Result<Self> {
        let _ = (lo, hi);
        let sections = (0..lens.len()).map(|i| read(i, 0, lens[i])).collect::<std::io::Result<Vec<_>>>()?;
        Ok(Self::from_sections(&sections, params))
    }
    /// Cells `[from, to)` holding the values of keys `[lo, hi)` -- possibly a few more -- found from
    /// what [`read_window`](Self::read_window) reads.
    fn window_cells(&self, lo: u64, hi: u64) -> (u64, u64) {
        let (mut from, mut to) = (u64::MAX, 0);
        self.for_each_in(lo, hi, &mut |_, start, slots| {
            from = from.min(start);
            to = to.max(start + slots);
        });
        if from > to { (0, 0) } else { (from, to) }
    }
    /// `n + 1` non-decreasing keys from 0 to the key space, cutting it into `n` ranges whose values
    /// take about equal cells (`total_slots / n`; a range may be empty).
    fn split(&self, n: usize, total_slots: u64) -> Vec<u64> {
        let _ = (n, total_slots);
        unimplemented!("splitting a key table of kind {}", Self::KIND)
    }
}

/// Flexmap's direct key table: a 16-bit offset for every key, a 64-bit head per CPB keys.
pub struct DirectKeys<const C: usize, const CPB: u64> {
    cells: Section,
}

impl<const C: usize, const CPB: u64> DirectKeys<C, CPB> {
    const HEAD: u64 = 4;

    pub fn from_fmkeys(keys: &FMKeys<C, CPB>) -> Self {
        let bytes = keys.data.len() * 2;
        let mut words = zeroed_words((bytes + 7) / 8);
        bytemuck::cast_slice_mut::<u64, u8>(&mut words)[..bytes].copy_from_slice(bytemuck::cast_slice(&keys.data));
        Self { cells: Section::from_words(words, bytes) }
    }

    #[inline(always)]
    fn block_index(key: u64) -> u64 {
        (key >> CPB.ilog2()) * (Self::HEAD + CPB)
    }

    #[inline(always)]
    fn head(&self, index: u64) -> u64 {
        unsafe { (self.cells.ptr() as *const u16).add(index as usize).cast::<u64>().read_unaligned() }
    }

    #[inline(always)]
    fn cell(&self, index: u64) -> u64 {
        unsafe { *(self.cells.ptr() as *const u16).add(index as usize) as u64 }
    }
}

impl<const C: usize, const CPB: u64> SlotLocator for DirectKeys<C, CPB> {
    const KIND: u32 = 1;

    #[inline(always)]
    fn locate(&self, key: u64) -> Option<(u64, u64)> {
        let mask = CPB - 1;
        let block = Self::block_index(key);
        let index = block + Self::HEAD + (key & mask);
        let base = self.head(block);
        let start = base + self.cell(index);
        let end = if key & mask != mask { base + self.cell(index + 1) } else { self.head(block + Self::HEAD + CPB) };
        if end <= start {
            return None;
        }
        Some((start, end - start))
    }

    #[inline(always)]
    fn prefetch_key(&self, key: u64) {
        prefetch(unsafe { self.cells.ptr().add(2 * Self::block_index(key) as usize) });
    }

    fn memory_bytes(&self) -> usize {
        self.cells.len()
    }

    fn sections(&self) -> Vec<Section> {
        vec![self.cells.clone()]
    }

    fn params(&self) -> [u64; 4] {
        [C as u64, CPB, 0, 0]
    }

    fn from_sections(sections: &[Section], params: [u64; 4]) -> Self {
        assert!(params[0] == C as u64 && params[1] == CPB, "direct key table of other parameters");
        assert_eq!(sections[0].len(), FMKeys::<C, CPB>::table_size() as usize * 2, "direct key table size");
        Self { cells: sections[0].clone() }
    }
}

/// The key table of the used keys alone (see the module's documentation).
pub struct CompactKeys {
    keys: u64,
    used: u64,
    /// 8 words a line: 7 words of presence bits (448 keys), then the used keys before the line.
    lines: Section,
    /// The start of every 64th rank.
    bases: Section,
    /// 3 bytes a rank (0..=used): its start less its group's base.
    offsets: Section,
}

impl CompactKeys {
    pub const KEYS_PER_LINE: u64 = 448;
    const GROUP: u64 = 64;

    /// The table of `used` keys, out of `keys` possible ones, listed by `each_key` in key order as
    /// (key, first cell, cells); the keys' cells follow one another.
    pub fn build(keys: u64, used: u64, each_key: &dyn Fn(&mut dyn FnMut(u64, u64, u64))) -> Self {
        let line_count = (keys + Self::KEYS_PER_LINE - 1) / Self::KEYS_PER_LINE;
        let mut lines = zeroed_words(line_count as usize * 8);
        let ranks = used + 1;
        let mut bases = vec![0u64; ((ranks + Self::GROUP - 1) / Self::GROUP) as usize];
        let offset_bytes = 3 * ranks as usize + 8;
        let mut offset_words = zeroed_words((offset_bytes + 7) / 8);
        let mut rank = 0u64;
        let mut next = 0u64;
        {
            let offsets = &mut bytemuck::cast_slice_mut::<u64, u8>(&mut offset_words)[..offset_bytes];
            let mut put_start = |rank: u64, start: u64| {
                let group = (rank / Self::GROUP) as usize;
                if rank % Self::GROUP == 0 {
                    bases[group] = start;
                }
                let off = start - bases[group];
                assert!(off < 1 << 24, "a group of {} keys spans {} cells (at most 2^24)", Self::GROUP, off);
                let r = 3 * rank as usize;
                offsets[r..r + 3].copy_from_slice(&(off as u32).to_le_bytes()[..3]);
            };
            each_key(&mut |key, start, slots| {
                assert!(rank == 0 || start == next, "the keys' cells do not follow one another");
                let line = (key / Self::KEYS_PER_LINE) as usize;
                let i = key % Self::KEYS_PER_LINE;
                lines[line * 8 + (i / 64) as usize] |= 1u64 << (i % 64);
                put_start(rank, start);
                next = start + slots;
                rank += 1;
            });
            assert_eq!(rank, used, "used keys");
            put_start(used, next);
        }
        let mut before = 0u64;
        for line in lines.chunks_exact_mut(8) {
            line[7] = before;
            before += line[..7].iter().map(|w| w.count_ones() as u64).sum::<u64>();
        }
        let bases_len = bases.len() * 8;
        Self {
            keys,
            used,
            lines: Section::from_words(lines, line_count as usize * 64),
            bases: Section::from_words(bases, bases_len),
            offsets: Section::from_words(offset_words, offset_bytes),
        }
    }

    pub fn used(&self) -> u64 {
        self.used
    }

    /// Calls `f(key, first cell, cells)` for every used key, in key order.
    pub fn for_each_used(&self, f: &mut dyn FnMut(u64, u64, u64)) {
        let lines = (self.keys + Self::KEYS_PER_LINE - 1) / Self::KEYS_PER_LINE;
        let mut rank = 0u64;
        for line in 0..lines {
            let p = unsafe { (self.lines.ptr() as *const u64).wrapping_add(line as usize * 8) };
            for w in 0..7u64 {
                let mut bits = unsafe { *p.add(w as usize) };
                while bits != 0 {
                    let b = bits.trailing_zeros() as u64;
                    bits &= bits - 1;
                    let start = self.start(rank);
                    f(line * Self::KEYS_PER_LINE + w * 64 + b, start, self.start(rank + 1) - start);
                    rank += 1;
                }
            }
        }
    }

    #[inline(always)]
    fn line(&self, key: u64) -> *const u64 {
        unsafe { (self.lines.ptr() as *const u64).wrapping_add((key / Self::KEYS_PER_LINE * 8) as usize) }
    }

    /// The key's index among the used keys, if it is used.
    #[inline(always)]
    pub fn rank(&self, key: u64) -> Option<u64> {
        let p = self.line(key);
        let i = key % Self::KEYS_PER_LINE;
        let w = (i / 64) as usize;
        let word = unsafe { *p.add(w) };
        if (word >> (i % 64)) & 1 == 0 {
            return None;
        }
        let mut rank = unsafe { *p.add(7) };
        for j in 0..7 {
            let c = unsafe { *p.add(j) }.count_ones() as u64;
            rank += if j < w { c } else { 0 };
        }
        Some(rank + (word & ((1u64 << (i % 64)) - 1)).count_ones() as u64)
    }

    /// The used keys before `line`, and the first rank past `line` (its own keys counted).
    fn line_ranks(&self, line: u64) -> (u64, u64) {
        let p = (self.lines.ptr() as *const u64).wrapping_add(line as usize * 8);
        let before = unsafe { *p.add(7) };
        (before, before + (0..7).map(|w| unsafe { *p.add(w) }.count_ones() as u64).sum::<u64>())
    }

    /// Ranks `[from, to]` whose starts keys `[lo, hi)` read: from the first key of `lo`'s line to the
    /// end of `hi - 1`'s line.
    fn window_ranks(&self, lo: u64, hi: u64) -> (u64, u64) {
        (self.line_ranks(lo / Self::KEYS_PER_LINE).0, self.line_ranks((hi - 1) / Self::KEYS_PER_LINE).1)
    }

    #[inline(always)]
    fn start(&self, rank: u64) -> u64 {
        let base = unsafe { *(self.bases.ptr() as *const u64).wrapping_add((rank / Self::GROUP) as usize) };
        let off = unsafe { self.offsets.ptr().wrapping_add(3 * rank as usize).cast::<u32>().read_unaligned() } & 0xFF_FFFF;
        base + off as u64
    }
}

impl SlotLocator for CompactKeys {
    const KIND: u32 = 2;

    #[inline(always)]
    fn locate(&self, key: u64) -> Option<(u64, u64)> {
        let rank = self.rank(key)?;
        let start = self.start(rank);
        Some((start, self.start(rank + 1) - start))
    }

    #[inline(always)]
    fn prefetch_key(&self, key: u64) {
        prefetch(self.line(key) as *const u8);
    }

    #[inline(always)]
    fn prefetch_key_second(&self, key: u64) {
        if let Some(rank) = self.rank(key) {
            prefetch(unsafe { self.offsets.ptr().wrapping_add(3 * rank as usize) });
            prefetch(unsafe { self.bases.ptr().wrapping_add(8 * (rank / Self::GROUP) as usize) });
        }
    }

    fn memory_bytes(&self) -> usize {
        self.lines.resident() + self.bases.resident() + self.offsets.resident()
    }

    fn sections(&self) -> Vec<Section> {
        vec![self.lines.clone(), self.bases.clone(), self.offsets.clone()]
    }

    fn params(&self) -> [u64; 4] {
        [self.keys, self.used, 0, 0]
    }

    fn from_sections(sections: &[Section], params: [u64; 4]) -> Self {
        let (keys, used) = (params[0], params[1]);
        let lines = (keys + Self::KEYS_PER_LINE - 1) / Self::KEYS_PER_LINE;
        assert_eq!(sections[0].len() as u64, lines * 64, "presence lines");
        assert_eq!(sections[1].len() as u64, (used + 1 + Self::GROUP - 1) / Self::GROUP * 8, "bases");
        assert_eq!(sections[2].len() as u64, 3 * (used + 1) + 8, "offsets");
        Self { keys, used, lines: sections[0].clone(), bases: sections[1].clone(), offsets: sections[2].clone() }
    }

    /// The lines of keys `[lo, hi)`, then the bases and offsets of their ranks.
    fn read_window(
        params: [u64; 4],
        lens: &[usize],
        lo: u64,
        hi: u64,
        read: &mut dyn FnMut(usize, usize, usize) -> std::io::Result<Section>,
    ) -> std::io::Result<Self> {
        assert!(lo < hi && hi <= params[0], "key range {lo}..{hi} of {}", params[0]);
        let (first, last) = (lo / Self::KEYS_PER_LINE, (hi - 1) / Self::KEYS_PER_LINE);
        let lines = read(0, first as usize * 64, (last as usize + 1) * 64)?;
        let empty = Section::from_words(Vec::new(), 0);
        let partial = Self { keys: params[0], used: params[1], lines: lines.clone(), bases: empty.clone(), offsets: empty };
        let (r_lo, r_hi) = partial.window_ranks(lo, hi);
        let bases = read(1, 8 * (r_lo / Self::GROUP) as usize, (8 * (r_hi / Self::GROUP + 1) as usize).min(lens[1]))?;
        // `start` reads 4 bytes at 3 * rank.
        let offsets = read(2, 3 * r_lo as usize, (3 * r_hi as usize + 8).min(lens[2]))?;
        Ok(Self::from_sections(&[lines, bases, offsets], params))
    }

    fn window_cells(&self, lo: u64, hi: u64) -> (u64, u64) {
        let (r_lo, r_hi) = self.window_ranks(lo, hi);
        (self.start(r_lo), self.start(r_hi))
    }

    /// At line boundaries: a line's first cell is the start of the rank before it, and starts grow
    /// with the line.
    fn split(&self, n: usize, total_slots: u64) -> Vec<u64> {
        let lines = (self.keys + Self::KEYS_PER_LINE - 1) / Self::KEYS_PER_LINE;
        let mut bounds = vec![0u64];
        for i in 1..n as u64 {
            let target = (total_slots as u128 * i as u128 / n as u128) as u64;
            let (mut a, mut b) = (0u64, lines);
            while a < b {
                let m = a + (b - a) / 2;
                if self.start(self.line_ranks(m).0) < target { a = m + 1 } else { b = m }
            }
            bounds.push((a * Self::KEYS_PER_LINE).min(self.keys).max(*bounds.last().unwrap()));
        }
        bounds.push(self.keys);
        bounds
    }

    fn for_each_in(&self, lo: u64, hi: u64, f: &mut dyn FnMut(u64, u64, u64)) {
        let hi = hi.min(self.keys);
        if lo >= hi {
            return;
        }
        let first = lo / Self::KEYS_PER_LINE;
        let last = (hi - 1) / Self::KEYS_PER_LINE;
        for line in first..=last {
            let p = unsafe { (self.lines.ptr() as *const u64).wrapping_add(line as usize * 8) };
            let mut rank = unsafe { *p.add(7) };
            for w in 0..7u64 {
                let mut bits = unsafe { *p.add(w as usize) };
                while bits != 0 {
                    let b = bits.trailing_zeros() as u64;
                    bits &= bits - 1;
                    let key = line * Self::KEYS_PER_LINE + w * 64 + b;
                    if key >= lo && key < hi {
                        let start = self.start(rank);
                        f(key, start, self.start(rank + 1) - start);
                    }
                    rank += 1;
                }
            }
        }
    }
}

/// The used keys in a hash table of 64-byte lines, 8 entries a line, so that a lookup reads one line
/// (the compact table reads two, one after the other). Of L lines, key k's home is k L / keys (the
/// keys are mixed, `mask_coremer`, so they spread evenly; any L, so the table does not double at a
/// power of two); its entry holds the key's low bits (the tag: a home's keys are a run of at most
/// ceil(keys / L), whose low bits differ),
/// how many lines past its home it lies (0-31), its first cell and its number of cells (never 0, so
/// an empty entry is 0). Keys go in key order, so each into the first line at or after its home with
/// room. The top bit of a line's first entry says that a key homed at or before the line lies past
/// it; a lookup reads the next line only then, so most misses read one line. At most 0.75 of the
/// entries are used: 8 / 0.75 = 10.7 bytes a used key.
pub struct HashedKeys {
    keys: u64,
    used: u64,
    homes: u64,
    /// floor(homes * 2^64 / keys): a key's home line is the high word of key * mult, so the key
    /// space need not be a power of two. For one that is, this is the old `key * homes >> key_bits`.
    mult: u64,
    tag_bits: u32,
    start_bits: u32,
    lines: Section,
}

impl HashedKeys {
    const ENTRIES: u64 = 8;
    const DISP_BITS: u32 = 5;
    const COUNT_BITS: u32 = 13;
    const MAX_LOAD: f64 = 0.75;
    /// In a line's first entry: keys homed at or before this line lie past it.
    const OVERFLOW: u64 = 1 << 63;

    /// The table of `used` keys (keys < `keys`), listed by `each_key` in key order.
    /// None if its fields do not fit 64 bits or a key would lie more than 31 lines past its home.
    pub fn build(keys: u64, used: u64, total_slots: u64, each_key: &dyn Fn(&mut dyn FnMut(u64, u64, u64))) -> Option<Self> {
        let homes = ((used as f64 / (Self::MAX_LOAD * Self::ENTRIES as f64)).ceil() as u64).max(1);
        let (mult, tag_bits) = Self::mult_and_tag_bits(keys, homes);
        let start_bits = bits_for(total_slots);
        if tag_bits + Self::DISP_BITS + start_bits + Self::COUNT_BITS > 63 {
            return None;
        }
        // One spare line past the last home, as a key may lie past it (lookups wrap to it, not around).
        let line_count = homes as usize + 32;
        let mut words = zeroed_words(line_count * 8);
        let mut cursor = 0u64;
        let mut fill = 0u64;
        let mut ok = true;
        each_key(&mut |key, start, slots| {
            let home = Self::home_of(key, mult);
            if home > cursor {
                cursor = home;
                fill = 0;
            }
            if fill == Self::ENTRIES {
                cursor += 1;
                fill = 0;
            }
            let disp = cursor - home;
            if disp > (1 << Self::DISP_BITS) - 1 || slots >= 1 << Self::COUNT_BITS {
                ok = false;
                return;
            }
            let entry = (key & ((1 << tag_bits) - 1))
                | disp << tag_bits
                | start << (tag_bits + Self::DISP_BITS)
                | slots << (tag_bits + Self::DISP_BITS + start_bits);
            words[(cursor * Self::ENTRIES + fill) as usize] = entry;
            for line in home..cursor {
                words[(line * Self::ENTRIES) as usize] |= Self::OVERFLOW;
            }
            fill += 1;
        });
        if !ok || cursor as usize >= line_count {
            return None;
        }
        Some(Self { keys, used, homes, mult, tag_bits, start_bits, lines: Section::from_words(words, line_count * 64) })
    }

    /// The home multiplier, and the tag bits that tell apart the keys sharing a home: a home holds
    /// at most ceil(2^64 / mult) consecutive keys.
    /// A key space no larger than the home lines (a slice of a key or two) is taken as one more than
    /// them, so every key still has a home line of its own.
    fn mult_and_tag_bits(keys: u64, homes: u64) -> (u64, u32) {
        let keys = keys.max(homes + 1);
        let mult = (((homes as u128) << 64) / keys as u128) as u64;
        let per_home = ((1u128 << 64) + mult as u128 - 1) / mult as u128;
        (mult, bits_for(per_home as u64))
    }

    #[inline(always)]
    fn home_of(key: u64, mult: u64) -> u64 {
        ((key as u128 * mult as u128) >> 64) as u64
    }

    #[inline(always)]
    fn line(&self, index: u64) -> *const u64 {
        unsafe { (self.lines.ptr() as *const u64).wrapping_add((index * Self::ENTRIES) as usize) }
    }
}

impl SlotLocator for HashedKeys {
    const KIND: u32 = 4; // 3: power-of-two homes (2026-10-23, never released)

    #[inline(always)]
    fn locate(&self, key: u64) -> Option<(u64, u64)> {
        let home = Self::home_of(key, self.mult);
        let tag = key & ((1 << self.tag_bits) - 1);
        let low = self.tag_bits + Self::DISP_BITS;
        let low_mask = (1u64 << low) - 1;
        for disp in 0..(1u64 << Self::DISP_BITS) {
            let p = self.line(home + disp);
            let want = tag | disp << self.tag_bits;
            for i in 0..Self::ENTRIES as usize {
                let e = unsafe { *p.add(i) } & !Self::OVERFLOW;
                if e & low_mask == want && e != 0 {
                    let start = (e >> low) & ((1u64 << self.start_bits) - 1);
                    return Some((start, (e >> (low + self.start_bits)) & ((1 << Self::COUNT_BITS) - 1)));
                }
            }
            if unsafe { *p } & Self::OVERFLOW == 0 {
                return None;
            }
        }
        None
    }

    #[inline(always)]
    fn prefetch_key(&self, key: u64) {
        prefetch(self.line(Self::home_of(key, self.mult)) as *const u8);
    }

    fn memory_bytes(&self) -> usize {
        self.lines.resident()
    }

    fn sections(&self) -> Vec<Section> {
        vec![self.lines.clone()]
    }

    fn params(&self) -> [u64; 4] {
        [self.keys, self.used, self.homes, self.start_bits as u64]
    }

    fn from_sections(sections: &[Section], params: [u64; 4]) -> Self {
        let (keys, used, homes, start_bits) = (params[0], params[1], params[2], params[3] as u32);
        assert_eq!(sections[0].len(), (homes as usize + 32) * 64, "hashed key lines");
        let (mult, tag_bits) = Self::mult_and_tag_bits(keys, homes);
        Self { keys, used, homes, mult, tag_bits, start_bits, lines: sections[0].clone() }
    }

    /// The lines a lookup of keys `[lo, hi)` reads: from `lo`'s home to 31 lines past `hi - 1`'s.
    fn read_window(
        params: [u64; 4],
        _lens: &[usize],
        lo: u64,
        hi: u64,
        read: &mut dyn FnMut(usize, usize, usize) -> std::io::Result<Section>,
    ) -> std::io::Result<Self> {
        assert!(lo < hi && hi <= params[0], "key range {lo}..{hi} of {}", params[0]);
        let (first, last) = Self::window_lines(params[0], params[2], lo, hi);
        let lines = read(0, first as usize * 64, last as usize * 64)?;
        Ok(Self::from_sections(&[lines], params))
    }

    /// The entries of the lines whose home is that of a key in `[lo, hi)`.
    fn window_cells(&self, lo: u64, hi: u64) -> (u64, u64) {
        let (first, last) = Self::window_lines(self.keys, self.homes, lo, hi);
        let (h_lo, h_hi) = (Self::home_of(lo, self.mult), Self::home_of(hi - 1, self.mult));
        let (mut from, mut to) = (u64::MAX, 0);
        for line in first..last {
            self.each_entry(line, &mut |home, start, slots| {
                if home >= h_lo && home <= h_hi {
                    from = from.min(start);
                    to = to.max(start + slots);
                }
            });
        }
        if from > to { (0, 0) } else { (from, to) }
    }

    /// At home boundaries: a key range is then the keys of a run of homes, and the first cell of a
    /// home grows with the home (keys go in key order). One pass over the lines.
    fn split(&self, n: usize, total_slots: u64) -> Vec<u64> {
        // first[h]: the least first cell of the homes from h on.
        let mut first = vec![u64::MAX; self.homes as usize + 1];
        first[self.homes as usize] = total_slots;
        for line in 0..self.homes + 32 {
            self.each_entry(line, &mut |home, start, _| {
                let f = &mut first[home as usize];
                *f = (*f).min(start);
            });
        }
        for h in (0..self.homes as usize).rev() {
            first[h] = first[h].min(first[h + 1]);
        }
        let mut bounds = vec![0u64];
        for i in 1..n as u64 {
            let target = (total_slots as u128 * i as u128 / n as u128) as u64;
            let home = first.partition_point(|&f| f < target) as u64;
            bounds.push(self.first_key_of(home).max(*bounds.last().unwrap()));
        }
        bounds.push(self.keys);
        bounds
    }
}

impl HashedKeys {
    /// Lines `[first, last)` that lookups of keys `[lo, hi)` read.
    fn window_lines(keys: u64, homes: u64, lo: u64, hi: u64) -> (u64, u64) {
        let (mult, _) = Self::mult_and_tag_bits(keys, homes);
        (Self::home_of(lo, mult), (Self::home_of(hi - 1, mult) + 32).min(homes + 32))
    }

    /// The least key homed at `home` or later (the key space past the last home).
    fn first_key_of(&self, home: u64) -> u64 {
        if home >= self.homes {
            return self.keys;
        }
        // The least k with k * mult >> 64 >= home: ceil(home * 2^64 / mult), corrected for rounding.
        let mut k = ((((home as u128) << 64) + self.mult as u128 - 1) / self.mult as u128) as u64;
        while k > 0 && Self::home_of(k - 1, self.mult) >= home {
            k -= 1;
        }
        while Self::home_of(k, self.mult) < home {
            k += 1;
        }
        k.min(self.keys)
    }

    /// Calls `f(home, first cell, cells)` for every entry of `line`.
    fn each_entry(&self, line: u64, f: &mut dyn FnMut(u64, u64, u64)) {
        let p = self.line(line);
        let low = self.tag_bits + Self::DISP_BITS;
        for i in 0..Self::ENTRIES as usize {
            let e = unsafe { *p.add(i) } & !Self::OVERFLOW;
            if e == 0 {
                continue;
            }
            let disp = (e >> self.tag_bits) & ((1 << Self::DISP_BITS) - 1);
            let start = (e >> low) & ((1u64 << self.start_bits) - 1);
            f(line - disp, start, (e >> (low + self.start_bits)) & ((1 << Self::COUNT_BITS) - 1));
        }
    }
}

const PACKED_MAGIC: [u8; 8] = *b"FMPACK01";
const PACKED_VERSION: u32 = 1;
/// Header bytes; the sections follow, each at a multiple of this (page-aligned).
const PACKED_ALIGN: usize = 4096;

/// A Flexmap in packed values with key table `K`.
pub struct PackedFlexmap<K: SlotLocator, const C: usize, const F: usize, const HEADER_THRESHOLD: usize> {
    pub keys: K,
    pub values: PackedValues,
}

pub type PackedFlexmapStd = PackedFlexmap<CompactKeys, 15, 16, 2>;

/// Lists a Flexmap's keys with values in key order, as (key, first cell, cells).
/// Lists the keys with values of a Flexmap key table (its cells, as FMKeys or a blob hold them) in key
/// order, as (key, first cell, cells).
fn table_keys<const C: usize, const CPB: u64>(cells: &[KCell]) -> impl Fn(&mut dyn FnMut(u64, u64, u64)) + '_ {
    assert_eq!(cells.len() as u64, FMKeys::<C, CPB>::table_size(), "key table of other C or CPB");
    let head = |i: u64| -> u64 { (0..4).map(|j| (cells[(i + j) as usize].0 as u64) << (16 * j)).fold(0, |a, b| a | b) };
    move |f: &mut dyn FnMut(u64, u64, u64)| {
        let block_cells = 4 + CPB;
        for block in 0..(1u64 << (2 * C)) / CPB {
            let b = block * block_cells;
            let (base, next) = (head(b), head(b + block_cells));
            for j in 0..CPB {
                let start = base + cells[(b + 4 + j) as usize].0 as u64;
                let end = if j + 1 < CPB { base + cells[(b + 5 + j) as usize].0 as u64 } else { next };
                if end > start {
                    f(block * CPB + j, start, end - start);
                }
            }
        }
    }
}

impl<const C: usize, const F: usize, const T: usize> PackedFlexmap<CompactKeys, C, F, T> {
    /// `flexmap` with packed values and a compact key table.
    pub fn from_flexmap<const CPB: u64>(flexmap: &Flexmap<C, F, CPB, T>) -> Self {
        Self::from_cells::<CPB>(&flexmap.keys.data, &flexmap.values.data)
    }

    /// A Flexmap blob (as `Flexmap::load_blob` or `FlexmapBlob::mmap_from_file` give it) with packed
    /// values and a compact key table: converts an existing index without rebuilding it.
    pub fn from_blob<const CPB: u64>(blob: &FlexmapBlob<C, F, CPB, T>) -> Self {
        Self::from_cells::<CPB>(blob.keys(), blob.values())
    }

    /// A Flexmap's key table and value cells.
    pub fn from_cells<const CPB: u64>(keys: &[KCell], values: &[VCell]) -> Self {
        let each_key = table_keys::<C, CPB>(keys);
        let mut used = 0u64;
        each_key(&mut |_, _, _| used += 1);
        let keys = CompactKeys::build(1u64 << (2 * C), used, &each_key);
        Self { keys, values: PackedValues::from_cells::<T>(values, &each_key) }
    }
}

impl<const C: usize, const F: usize, const T: usize> PackedFlexmap<HashedKeys, C, F, T> {
    /// The same index with a line-hashed key table over the same key space (a slice's is its
    /// range); shares `compact`'s values (no copy).
    pub fn from_compact(compact: &PackedFlexmap<CompactKeys, C, F, T>) -> Option<Self> {
        let each_key = |f: &mut dyn FnMut(u64, u64, u64)| compact.keys.for_each_used(f);
        let keys = HashedKeys::build(compact.keys.keys, compact.keys.used(), compact.values.total_slots, &each_key)?;
        let values = PackedValues { layout: compact.values.layout, total_slots: compact.values.total_slots, data: compact.values.data.clone() };
        Some(Self { keys, values })
    }

    /// `flexmap` with packed values and a line-hashed key table; None if the table's fields do not fit.
    pub fn from_flexmap_hashed<const CPB: u64>(flexmap: &Flexmap<C, F, CPB, T>) -> Option<Self> {
        let each_key = table_keys::<C, CPB>(&flexmap.keys.data);
        let mut used = 0u64;
        each_key(&mut |_, _, _| used += 1);
        let keys = HashedKeys::build(1u64 << (2 * C), used, flexmap.values.data.len() as u64, &each_key)?;
        Some(Self { keys, values: PackedValues::from_cells::<T>(&flexmap.values.data, &each_key) })
    }
}

impl<const C: usize, const F: usize, const T: usize, const CPB: u64> PackedFlexmap<DirectKeys<C, CPB>, C, F, T> {
    /// `flexmap` with packed values and its direct key table.
    pub fn from_flexmap_direct(flexmap: &Flexmap<C, F, CPB, T>) -> Self {
        let each_key = table_keys::<C, CPB>(&flexmap.keys.data);
        Self { keys: DirectKeys::from_fmkeys(&flexmap.keys), values: PackedValues::from_cells::<T>(&flexmap.values.data, &each_key) }
    }
}

fn put_u32(out: &mut Vec<u8>, x: u32) {
    out.extend_from_slice(&x.to_le_bytes());
}

fn put_u64(out: &mut Vec<u8>, x: u64) {
    out.extend_from_slice(&x.to_le_bytes());
}

struct HeaderReader<'a>(&'a [u8], usize);

impl HeaderReader<'_> {
    fn u32(&mut self) -> u32 {
        let x = u32::from_le_bytes(self.0[self.1..self.1 + 4].try_into().unwrap());
        self.1 += 4;
        x
    }

    fn u64(&mut self) -> u64 {
        let x = u64::from_le_bytes(self.0[self.1..self.1 + 8].try_into().unwrap());
        self.1 += 8;
        x
    }
}

fn align_up(x: usize) -> usize {
    (x + PACKED_ALIGN - 1) / PACKED_ALIGN * PACKED_ALIGN
}

impl<K: SlotLocator, const C: usize, const F: usize, const T: usize> PackedFlexmap<K, C, F, T> {
    /// The values of `key`'s cells, if it has any.
    #[inline(always)]
    pub fn get(&self, key: u64) -> Option<PackedRange<'_>> {
        let (start, slots) = self.keys.locate(key)?;
        Some(self.values.range(start, slots))
    }

    /// Calls `f(i, get(queries[i]))` for every query in order, fetching ahead: the key table of the
    /// query `3 * step` ahead, its second read (a compact table's start) `2 * step` ahead, and the
    /// values of the query `step` ahead. The same calls as a loop over `get`.
    #[inline]
    pub fn for_each_prefetched_with<'s, Fun>(&'s self, queries: &[u64], step: usize, mut f: Fun)
    where
        Fun: FnMut(usize, Option<PackedRange<'s>>),
    {
        const RING: usize = 64;
        let step = step.clamp(1, RING / 3);
        let n = queries.len();
        let mut ring = [None::<(u64, u64)>; RING];
        for i in 0..n + 3 * step {
            if i < n {
                self.keys.prefetch_key(queries[i]);
            }
            if i >= step && i - step < n {
                self.keys.prefetch_key_second(queries[i - step]);
            }
            if i >= 2 * step && i - 2 * step < n {
                let j = i - 2 * step;
                let loc = self.keys.locate(queries[j]);
                if let Some((start, slots)) = loc {
                    self.values.prefetch(start, slots);
                }
                ring[j % RING] = loc;
            }
            if i >= 3 * step && i - 3 * step < n {
                let j = i - 3 * step;
                f(j, ring[j % RING].map(|(start, slots)| self.values.range(start, slots)));
            }
        }
    }

    /// [`for_each_prefetched_with`](Self::for_each_prefetched_with) at the default distance (8).
    #[inline]
    pub fn for_each_prefetched<'s, Fun>(&'s self, queries: &[u64], f: Fun)
    where
        Fun: FnMut(usize, Option<PackedRange<'s>>),
    {
        self.for_each_prefetched_with(queries, 8, f)
    }

    /// Calls `f(key, first cell, cells)` for every key in `[lo, hi)` with values, in key order.
    pub fn for_each_key_in(&self, lo: u64, hi: u64, f: &mut dyn FnMut(u64, u64, u64)) {
        self.keys.for_each_in(lo, hi, f)
    }

    /// The keys `[lo, hi)` alone, as keys `0..hi - lo`, with a compact key table and only their
    /// values: `slice(lo, hi).get(k - lo)` answers as `get(k)` for every `k` in `[lo, hi)` (and must
    /// not be asked for keys outside `0..hi - lo`, as no table may be asked past its keys). The values
    /// keep this index's widths, so slices of one index can be read alike. Reads only the keys and
    /// values of the range; the result owns its memory.
    pub fn slice(&self, lo: u64, hi: u64) -> PackedFlexmap<CompactKeys, C, F, T> {
        assert!(lo <= hi, "key range {lo}..{hi}");
        let (mut used, mut total_slots) = (0u64, 0u64);
        self.keys.for_each_in(lo, hi, &mut |_, _, slots| {
            used += 1;
            total_slots += slots;
        });
        let each_key = |f: &mut dyn FnMut(u64, u64, u64)| {
            let mut start = 0u64;
            self.keys.for_each_in(lo, hi, &mut |key, _, slots| {
                f(key - lo, start, slots);
                start += slots;
            });
        };
        let keys = CompactKeys::build(hi - lo, used, &each_key);
        let values = PackedValues::from_packed(&self.values, total_slots, &|f: &mut dyn FnMut(u64, u64)| {
            self.keys.for_each_in(lo, hi, &mut |_, start, slots| f(start, slots));
        });
        PackedFlexmap { keys, values }
    }

    /// [`slice`](Self::slice) with a line-hashed key table over keys `0..hi - lo`; None if the
    /// table's fields do not fit.
    pub fn slice_hashed(&self, lo: u64, hi: u64) -> Option<PackedFlexmap<HashedKeys, C, F, T>> {
        PackedFlexmap::<HashedKeys, C, F, T>::from_compact(&self.slice(lo, hi))
    }

    /// Bytes of the key table and of the values.
    pub fn memory_bytes(&self) -> (usize, usize) {
        (self.keys.memory_bytes(), self.values.memory_bytes())
    }

    pub fn save(&self, path: &str) -> std::io::Result<()> {
        let mut file = std::io::BufWriter::with_capacity(8 << 20, File::create(path)?);
        self.write_to(&mut file)?;
        file.flush()
    }

    /// Writes the index as [`save`](Self::save) does, to any writer; returns the bytes written (a
    /// multiple of 4096). Inside a larger file the index must start at a multiple of 4096, so its
    /// sections keep their alignment; [`from_shared`](Self::from_shared) then reads it there.
    pub fn write_to<W: Write>(&self, file: &mut W) -> std::io::Result<u64> {
        let mut sections = self.keys.sections();
        sections.push(self.values.data.clone());
        let l = &self.values.layout;
        let mut header = Vec::new();
        header.extend_from_slice(&PACKED_MAGIC);
        for x in [PACKED_VERSION, C as u32, F as u32, T as u32, K::KIND, sections.len() as u32, l.ref_bits, l.pos_bits] {
            put_u32(&mut header, x);
        }
        put_u64(&mut header, l.slot_bits32);
        put_u64(&mut header, self.values.total_slots);
        for x in self.keys.params() {
            put_u64(&mut header, x);
        }
        let mut offset = PACKED_ALIGN;
        for s in &sections {
            put_u64(&mut header, offset as u64);
            put_u64(&mut header, s.len() as u64);
            offset = align_up(offset + s.len());
        }
        assert!(header.len() <= PACKED_ALIGN);
        header.resize(PACKED_ALIGN, 0);
        file.write_all(&header)?;
        let mut at = PACKED_ALIGN;
        for s in &sections {
            file.write_all(s.bytes())?;
            let end = align_up(at + s.len());
            file.write_all(&vec![0u8; end - at - s.len()])?;
            at = end;
        }
        Ok(at as u64)
    }

    fn from_backing(backing: Backing) -> Self {
        Self::from_bytes_at(&SharedBytes::from_arc(Section::backing(backing)), 0)
    }

    /// The index written by [`write_to`](Self::write_to) at `offset` (a multiple of 4096) inside
    /// `bytes`: a bundle file read or mapped once, of which the index is one part. Nothing is
    /// copied; the index keeps `bytes` alive.
    pub fn from_shared(bytes: &SharedBytes, offset: usize) -> Self {
        assert_eq!(offset % PACKED_ALIGN, 0, "a packed flexmap must start on a 4096-byte boundary");
        Self::from_bytes_at(bytes, offset)
    }

    fn from_bytes_at(bytes: &SharedBytes, base: usize) -> Self {
        let all = bytes.section(base, PACKED_ALIGN);
        let (layout, total_slots, params, extents) = Self::parse_header(all.bytes());
        let count = extents.len();
        let sections: Vec<Section> = extents.iter().map(|&(offset, len)| bytes.section(base + offset, len)).collect();
        let data = sections[count - 1].clone();
        assert!(data.len() as u64 >= ((layout.slot_bit(total_slots) + 7) >> 3) + VALUE_PADDING as u64, "values truncated");
        Self { keys: K::from_sections(&sections[..count - 1], params), values: PackedValues { layout, total_slots, data } }
    }

    /// The layout, cell count, key-table parameters and (offset, length) of every section, from a
    /// packed flexmap's first PACKED_ALIGN bytes.
    fn parse_header(head: &[u8]) -> (PackedLayout, u64, [u64; 4], Vec<(usize, usize)>) {
        let mut r = HeaderReader(head, 8);
        assert!(head[..8] == PACKED_MAGIC, "not a packed flexmap");
        assert_eq!(r.u32(), PACKED_VERSION, "packed flexmap version");
        assert!(
            r.u32() == C as u32 && r.u32() == F as u32 && r.u32() == T as u32,
            "packed flexmap of other C, F or HEADER_THRESHOLD"
        );
        assert_eq!(r.u32(), K::KIND, "packed flexmap of another key table");
        let count = r.u32() as usize;
        let (ref_bits, pos_bits) = (r.u32(), r.u32());
        let slot_bits32 = r.u64();
        let total_slots = r.u64();
        let params = [r.u64(), r.u64(), r.u64(), r.u64()];
        let extents = (0..count).map(|_| (r.u64() as usize, r.u64() as usize)).collect();
        (PackedLayout { ref_bits, pos_bits, header_threshold: T as u32, slot_bits32 }, total_slots, params, extents)
    }

    /// The header of the index written at `offset` of `file`, and its key table read whole.
    fn read_keys(file: &File, offset: u64, threads: usize) -> std::io::Result<(K, PackedLayout, u64, Vec<(usize, usize)>)> {
        let mut head = vec![0u8; PACKED_ALIGN];
        file.read_exact_at(&mut head, offset)?;
        let (layout, total_slots, params, extents) = Self::parse_header(&head);
        let mut sections = Vec::new();
        for &(at, len) in &extents[..extents.len() - 1] {
            sections.push(SharedBytes::read(file, offset + at as u64, len, threads)?.section(0, len));
        }
        Ok((K::from_sections(&sections, params), layout, total_slots, extents))
    }

    /// The keys `[lo, hi)` of the index written at `offset` of `file` (a multiple of 4096: a
    /// bundle's index section, or 0 for an index file of its own), read alone: of the key table
    /// the part their lookups read, of the values the cells their keys hold, each in one read on
    /// `threads` threads into huge pages. Nothing is rewritten or copied: the parts keep the whole
    /// index's addressing ([`Section::window`]), so `get(k)` answers as the whole index for every
    /// `k` in `[lo, hi)`, through the same code. It must not be asked for other keys.
    pub fn read_window(file: &File, offset: u64, lo: u64, hi: u64, threads: usize) -> std::io::Result<Self> {
        let mut head = vec![0u8; PACKED_ALIGN];
        file.read_exact_at(&mut head, offset)?;
        let (layout, total_slots, params, extents) = Self::parse_header(&head);
        let count = extents.len();
        let mut read = |i: usize, from: usize, to: usize| -> std::io::Result<Section> {
            let (at, len) = extents[i];
            let to = to.min(len);
            let from = from.min(to);
            let held = SharedBytes::read(file, offset + (at + from) as u64, to - from, threads)?;
            Ok(Section::window(&held, from, len))
        };
        let lens: Vec<usize> = extents[..count - 1].iter().map(|e| e.1).collect();
        let keys = K::read_window(params, &lens, lo, hi, &mut |i, from, to| read(i, from, to))?;
        let (first, last) = keys.window_cells(lo, hi);
        // The values of cells [first, last), as `range` reads them: from the 8-byte word holding
        // the first cell's first bit to VALUE_PADDING bytes past the last cell (unaligned reads).
        let from = (layout.slot_bit(first) >> 3) as usize & !7;
        let to = ((layout.slot_bit(last) + 7) >> 3) as usize + VALUE_PADDING;
        let data = read(count - 1, from, to)?;
        Ok(Self { keys, values: PackedValues { layout, total_slots, data } })
    }

    /// `n + 1` non-decreasing keys cutting the index at `offset` of `file` into `n` key ranges of
    /// about equal value bytes, for [`read_window`](Self::read_window). Reads only the key table.
    pub fn split_in_file(file: &File, offset: u64, n: usize, threads: usize) -> std::io::Result<Vec<u64>> {
        assert!(n > 0, "need at least one part");
        let (keys, _, total_slots, _) = Self::read_keys(file, offset, threads)?;
        Ok(keys.split(n, total_slots))
    }

    /// Reads the file into a huge-page buffer on `threads` threads.
    pub fn load(path: &str, threads: usize) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let size = file.metadata()?.len() as usize;
        let mut words = zeroed_words((size + 7) / 8);
        read_parallel(&file, 0, &mut bytemuck::cast_slice_mut::<u64, u8>(&mut words)[..size], threads)?;
        Ok(Self::from_backing(Backing::Words(words)))
    }

    /// Memory-maps the file: pages are read as lookups touch them.
    pub fn mmap(path: &str) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let map = unsafe { memmap2::Mmap::map(&file)? };
        Ok(Self::from_backing(Backing::Mapped(map)))
    }

    /// Whether `path` holds a packed flexmap these parameters read.
    pub fn is_compatible(path: &str) -> bool {
        let mut head = [0u8; 28];
        File::open(path).and_then(|mut f| f.read_exact(&mut head)).is_ok() && Self::is_compatible_header(&head)
    }

    /// Whether `head`, the first 28 bytes of a packed flexmap, is one these parameters read.
    pub fn is_compatible_header(head: &[u8]) -> bool {
        if head.len() < 28 {
            return false;
        }
        let mut r = HeaderReader(head, 8);
        head[..8] == PACKED_MAGIC
            && r.u32() == PACKED_VERSION
            && r.u32() == C as u32
            && r.u32() == F as u32
            && r.u32() == T as u32
            && r.u32() == K::KIND
    }
}

#[cfg(test)]
mod tests {
    /// Every used key of a key space that is not a power of two is found with its start and count,
    /// and every other key is absent; a power-of-two space gets the homes `key * homes >> bits` did.
    #[test]
    fn hashed_keys_any_key_space() {
        for keys in [1000u64, 4096, 999_983, 1 << 20] {
            let used: Vec<u64> = (0..keys).filter(|k| k % 7 == 3 || k % 11 == 0).collect();
            let each = |f: &mut dyn FnMut(u64, u64, u64)| {
                for (i, &k) in used.iter().enumerate() {
                    f(k, i as u64 * 3, 1 + (k % 5));
                }
            };
            let t = HashedKeys::build(keys, used.len() as u64, used.len() as u64 * 3, &each).expect("fits");
            let mut it = used.iter().peekable();
            for k in 0..keys {
                if it.peek() == Some(&&k) {
                    let i = used.iter().position(|&u| u == k).unwrap() as u64;
                    assert_eq!(t.locate(k), Some((i * 3, 1 + (k % 5))), "key {k} of {keys}");
                    it.next();
                } else {
                    assert_eq!(t.locate(k), None, "key {k} of {keys}");
                }
            }
            if keys.is_power_of_two() {
                let bits = keys.trailing_zeros();
                for k in (0..keys).step_by(97) {
                    assert_eq!(HashedKeys::home_of(k, t.mult), (k * t.homes) >> bits);
                }
            }
        }
    }

    use std::cmp::min;

    use kmerrs::consecutive::kmer::{Kmer, KmerIter};

    use super::*;
    use crate::flexmap::VRangeGetter;
    use crate::keys::FMKeys;

    /// A Flexmap of C = 6 over a random sequence repeated with mutations, so that keys have 1 to
    /// hundreds of values, with and without headers, and several reference ids.
    fn random_flexmap() -> Flexmap<6, 8, 8, 2> {
        let mut state = 0x9E3779B97F4A7C15u64;
        let mut next = move || {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            state.wrapping_mul(0x2545F4914F6CDD1D)
        };
        let base: Vec<u8> = (0..20_000).map(|_| b"ACGT"[(next() >> 62) as usize]).collect();
        let records: Vec<Vec<u8>> = (0..12)
            .map(|_| base.iter().map(|&c| if next() % 50 == 0 { b"ACGT"[(next() >> 62) as usize] } else { c }).collect())
            .collect();
        let mut keys = FMKeys::<6, 8>::new();
        for seq in &records {
            for (_, f, r) in KmerIter::<22, true>::new(seq) {
                let kmer: Kmer<22> = min(f, r);
                keys.get_kmer_cell_mut_ref(kmer.middle::<6>().0).increment();
            }
        }
        keys.build::<2>(4000);
        let mut flexmap = Flexmap::<6, 8, 8, 2>::new(keys);
        for (id, seq) in records.iter().enumerate() {
            for (pos, f, r) in KmerIter::<22, true>::new(seq) {
                let kmer: Kmer<22> = min(f, r);
                if let Some(range) = flexmap.keys.vrange(kmer.middle::<6>().0) {
                    flexmap.values.get_range_mut(range).insert(VD::set(id as u64 + 1, pos as u64), kmer.flanks::<8>().0 as u32);
                }
            }
        }
        flexmap
    }

    fn assert_same<K: SlotLocator>(flexmap: &Flexmap<6, 8, 8, 2>, packed: &PackedFlexmap<K, 6, 8, 2>) {
        let mut with_header = 0;
        for key in 0..(1u64 << 12) {
            match (flexmap.get_vrange(key), packed.get(key)) {
                (None, None) => {}
                (Some(a), Some(b)) => {
                    assert_eq!(a.positions.len(), b.len(), "key {key}: values");
                    for i in 0..b.len() {
                        assert_eq!(VD::get(a.positions[i].0), b.get(i), "key {key}: value {i}");
                    }
                    match a.header {
                        Some(h) => {
                            with_header += 1;
                            for i in 0..b.len() {
                                assert_eq!(Some(h[i].get()), b.flex(i), "key {key}: flank cell {i}");
                            }
                        }
                        None => assert!(!b.has_header(), "key {key}: header"),
                    }
                    for q in [0u64, 0xFFFF, key * 2654435761 & 0xFFFF] {
                        let flex = Kmer::<8>(q);
                        let (mut x, mut y) = (Vec::new(), Vec::new());
                        a.best_flex_match(&flex, |p, v, d| x.push((p, v, d)));
                        b.best_flex_match(&flex, |p, v, d| y.push((p, v, d)));
                        assert_eq!(x, y, "key {key}: best_flex_match");
                    }
                }
                (a, b) => panic!("key {key}: Flexmap has values {}, packed {}", a.is_some(), b.is_some()),
            }
        }
        assert!(with_header > 100, "the test needs keys with headers ({with_header})");
    }

    #[test]
    fn packed_flexmap_answers_as_flexmap() {
        let flexmap = random_flexmap();
        let compact = PackedFlexmap::<CompactKeys, 6, 8, 2>::from_flexmap(&flexmap);
        assert_same(&flexmap, &compact);
        let direct = PackedFlexmap::<DirectKeys<6, 8>, 6, 8, 2>::from_flexmap_direct(&flexmap);
        assert_same(&flexmap, &direct);
        let hashed = PackedFlexmap::<HashedKeys, 6, 8, 2>::from_flexmap_hashed(&flexmap).expect("hashed table fits");
        assert_same(&flexmap, &hashed);
        let from_compact = PackedFlexmap::<HashedKeys, 6, 8, 2>::from_compact(&compact).expect("hashed table fits");
        assert_same(&flexmap, &from_compact);
        for key in 0..(1u64 << 12) {
            assert_eq!(hashed.keys.locate(key), compact.keys.locate(key), "key {key}");
        }
        assert!(compact.values.memory_bytes() < flexmap.values.data.len() * 8);
    }

    #[test]
    fn slices_answer_as_the_whole_index() {
        let flexmap = random_flexmap();
        let compact = PackedFlexmap::<CompactKeys, 6, 8, 2>::from_flexmap(&flexmap);
        let hashed = PackedFlexmap::<HashedKeys, 6, 8, 2>::from_compact(&compact).expect("hashed table fits");
        let keys = 1u64 << 12;
        // Bounds that are not multiples of a line (448) or of anything else.
        let bounds = [0u64, 1, 447, 1000, 1001, 2900, keys];
        let mut total = 0u64;
        for w in bounds.windows(2) {
            let (lo, hi) = (w[0], w[1]);
            for slice in [compact.slice(lo, hi), hashed.slice(lo, hi)] {
                assert_eq!(slice.values.layout, compact.values.layout);
                for key in lo..hi {
                    match (compact.get(key), slice.get(key - lo)) {
                        (None, None) => {}
                        (Some(a), Some(b)) => {
                            assert_eq!((a.len(), a.has_header()), (b.len(), b.has_header()), "key {key}");
                            assert_eq!(a.flex_bytes(), b.flex_bytes(), "key {key}: flank cells");
                            for i in 0..a.len() {
                                assert_eq!(a.get(i), b.get(i), "key {key}: value {i}");
                            }
                        }
                        (a, b) => panic!("key {key}: whole {} slice {}", a.is_some(), b.is_some()),
                    }
                }
            }
            // Hashed slices: a key space of hi - lo, which is not a power of two.
            for slice in [compact.slice_hashed(lo, hi).expect("fits"), hashed.slice_hashed(lo, hi).expect("fits")] {
                for key in lo..hi {
                    match (compact.get(key), slice.get(key - lo)) {
                        (None, None) => {}
                        (Some(a), Some(b)) => {
                            assert_eq!((a.len(), a.has_header()), (b.len(), b.has_header()), "key {key}");
                            assert_eq!(a.flex_bytes(), b.flex_bytes(), "key {key}: flank cells");
                            for i in 0..a.len() {
                                assert_eq!(a.get(i), b.get(i), "key {key}: value {i}");
                            }
                        }
                        (a, b) => panic!("key {key}: whole {} hashed slice {}", a.is_some(), b.is_some()),
                    }
                }
            }
            let mut walked = Vec::new();
            hashed.for_each_key_in(lo, hi, &mut |k, s, n| walked.push((k, s, n)));
            let mut expect = Vec::new();
            compact.for_each_key_in(lo, hi, &mut |k, s, n| expect.push((k, s, n)));
            assert_eq!(walked, expect, "keys of {lo}..{hi}");
            total += compact.slice(lo, hi).values.total_slots();
        }
        assert_eq!(total, compact.values.total_slots(), "the slices partition the cells");
        let empty = compact.slice(5, 5);
        assert_eq!((empty.keys.used(), empty.values.total_slots()), (0, 0));
    }

    /// Windows read from a file (at an offset, as inside a bundle) answer as the whole index for
    /// every key of their range, for both key tables; the split covers the key space in order and
    /// the windows hold less than the whole.
    #[test]
    fn windows_answer_as_the_whole_index() {
        fn check<K: SlotLocator>(whole: &PackedFlexmap<K, 6, 8, 2>, name: &str) {
            let path = std::env::temp_dir().join(format!("flexmap-window-{}-{name}", std::process::id()));
            let offset = 3 * 4096u64;
            {
                let mut f = std::io::BufWriter::new(File::create(&path).unwrap());
                f.write_all(&vec![7u8; offset as usize]).unwrap();
                whole.write_to(&mut f).unwrap();
                f.flush().unwrap();
            }
            let file = File::open(&path).unwrap();
            let keys = 1u64 << 12;
            for n in [1usize, 2, 3, 7] {
                let bounds = PackedFlexmap::<K, 6, 8, 2>::split_in_file(&file, offset, n, 2).unwrap();
                assert_eq!((bounds.len(), bounds[0], bounds[n]), (n + 1, 0, keys), "{name} {n}: {bounds:?}");
                assert!(bounds.windows(2).all(|w| w[0] <= w[1]), "{name} {n}: {bounds:?}");
                for w in bounds.windows(2) {
                    let (lo, hi) = (w[0], w[1]);
                    if lo == hi {
                        continue;
                    }
                    let part = PackedFlexmap::<K, 6, 8, 2>::read_window(&file, offset, lo, hi, 2).unwrap();
                    if n > 1 {
                        assert!(part.values.memory_bytes() < whole.values.memory_bytes(), "{name} {n}: window holds all values");
                    }
                    for key in lo..hi {
                        match (whole.get(key), part.get(key)) {
                            (None, None) => {}
                            (Some(a), Some(b)) => {
                                assert_eq!((a.len(), a.has_header()), (b.len(), b.has_header()), "{name} key {key}");
                                assert_eq!(a.flex_bytes(), b.flex_bytes(), "{name} key {key}: flank cells");
                                for i in 0..a.len() {
                                    assert_eq!(a.get(i), b.get(i), "{name} key {key}: value {i}");
                                }
                            }
                            (a, b) => panic!("{name} key {key}: whole {} window {}", a.is_some(), b.is_some()),
                        }
                    }
                }
            }
            std::fs::remove_file(&path).unwrap();
        }
        let flexmap = random_flexmap();
        let compact = PackedFlexmap::<CompactKeys, 6, 8, 2>::from_flexmap(&flexmap);
        let hashed = PackedFlexmap::<HashedKeys, 6, 8, 2>::from_compact(&compact).expect("hashed table fits");
        check(&compact, "compact");
        check(&hashed, "hashed");
    }

    #[test]
    fn packed_flexmap_from_blob() {
        let flexmap = random_flexmap();
        let path = std::env::temp_dir().join(format!("flexmap_blob_to_packed_{}.bin", std::process::id()));
        let path = path.to_string_lossy().to_string();
        flexmap.save_blob(&path);
        let blob = crate::flexmap::FlexmapBlob::<6, 8, 8, 2>::mmap_from_file(&path);
        assert_same(&flexmap, &PackedFlexmap::<CompactKeys, 6, 8, 2>::from_blob(&blob));
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn packed_flexmap_files_round_trip() {
        let flexmap = random_flexmap();
        let packed = PackedFlexmap::<CompactKeys, 6, 8, 2>::from_flexmap(&flexmap);
        let path = std::env::temp_dir().join(format!("flexmap_packed_{}.bin", std::process::id()));
        let path = path.to_string_lossy().to_string();
        packed.save(&path).unwrap();
        assert!(PackedFlexmap::<CompactKeys, 6, 8, 2>::is_compatible(&path));
        assert!(!PackedFlexmap::<DirectKeys<6, 8>, 6, 8, 2>::is_compatible(&path));
        assert_same(&flexmap, &PackedFlexmap::<CompactKeys, 6, 8, 2>::load(&path, 3).unwrap());
        assert_same(&flexmap, &PackedFlexmap::<CompactKeys, 6, 8, 2>::mmap(&path).unwrap());
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn packed_flexmap_inside_a_larger_file() {
        // An index written after 8192 bytes of other data and followed by more: read back from the
        // shared bytes of the whole file, both read in and mapped, it answers as the original.
        let flexmap = random_flexmap();
        let packed = PackedFlexmap::<CompactKeys, 6, 8, 2>::from_flexmap(&flexmap);
        let mut file_bytes = vec![7u8; 8192];
        let written = packed.write_to(&mut file_bytes).unwrap();
        assert_eq!(written as usize, file_bytes.len() - 8192);
        assert_eq!(written % 4096, 0);
        file_bytes.extend_from_slice(b"trailing data");
        let path = std::env::temp_dir().join(format!("flexmap_shared_{}.bin", std::process::id()));
        std::fs::write(&path, &file_bytes).unwrap();
        let file = File::open(&path).unwrap();
        assert!(PackedFlexmap::<CompactKeys, 6, 8, 2>::is_compatible_header(&file_bytes[8192..]));
        for bytes in [SharedBytes::read(&file, 0, file_bytes.len(), 3).unwrap(), SharedBytes::mmap(&file).unwrap()] {
            assert_eq!(bytes.bytes(), &file_bytes[..]);
            let index = PackedFlexmap::<CompactKeys, 6, 8, 2>::from_shared(&bytes, 8192);
            drop(bytes); // the index keeps the bytes alive
            assert_same(&flexmap, &index);
        }
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn prefetched_lookups_match_get() {
        let flexmap = random_flexmap();
        let packed = PackedFlexmap::<CompactKeys, 6, 8, 2>::from_flexmap(&flexmap);
        let queries: Vec<u64> = (0..5000u64).map(|i| (i * 2654435761) & 0xFFF).collect();
        for step in [1, 3, 8, 21] {
            let mut seen = 0;
            packed.for_each_prefetched_with(&queries, step, |i, range| {
                assert_eq!(i, seen, "order");
                seen += 1;
                let expected = packed.get(queries[i]);
                assert_eq!(range.map(|r| (r.len(), r.get(0))), expected.map(|r| (r.len(), r.get(0))));
            });
            assert_eq!(seen, queries.len());
        }
    }

    #[test]
    fn layout_holds_every_key_size() {
        for (max_ref, max_pos) in [(1, 1), (1000, 5_000_000), ((1 << 28) - 1, (1 << 34) - 1), (3, (1 << 40) - 1)] {
            let l = PackedLayout::new(max_ref, max_pos, 6144, 2);
            let w = l.entry_bits() as u64;
            // The cells a key of n values takes: n, or n + ceil(n/2) with a header.
            for n in 1..4096u64 {
                let slots = if n > 2 { n + (n + 1) / 2 } else { n };
                assert_eq!(l.entries(slots), n);
                let need = if slots > 2 { n * (32 + w) + 7 } else { n * w };
                for start in [0u64, 1, 7, 12345, 999_999_937] {
                    assert!(l.slot_bit(start + slots) - l.slot_bit(start) >= need, "{max_ref} {max_pos}: {slots} cells at {start}");
                }
            }
        }
    }
}
