//! One range type for both indexes: a Flexmap's `VRange` (cells and header slices) or a packed
//! index's `PackedRange` (bit fields). Code that reads ranges through it works with either index.

use kmerrs::consecutive::kmer::Kmer;

use crate::packed::{flank_distance, PackedRange};
use crate::values::VRange;
use crate::VD;

/// A key's values from either index.
#[derive(Clone, Copy)]
pub enum AnyRange<'a> {
    Plain(VRange<'a>),
    Packed(PackedRange<'a>),
}

impl<'a> From<VRange<'a>> for AnyRange<'a> {
    fn from(r: VRange<'a>) -> Self {
        AnyRange::Plain(r)
    }
}

impl<'a> From<PackedRange<'a>> for AnyRange<'a> {
    fn from(r: PackedRange<'a>) -> Self {
        AnyRange::Packed(r)
    }
}

impl<'a> AnyRange<'a> {
    /// The key's values.
    #[inline(always)]
    pub fn len(&self) -> usize {
        match self {
            AnyRange::Plain(r) => r.positions.len(),
            AnyRange::Packed(r) => r.len(),
        }
    }

    /// Whether the key has flank cells (more than HEADER_THRESHOLD values).
    #[inline(always)]
    pub fn has_header(&self) -> bool {
        match self {
            AnyRange::Plain(r) => r.header.is_some(),
            AnyRange::Packed(r) => r.has_header(),
        }
    }

    /// Value `i` as `(reference id, position)`, as `VD::get` gives a cell.
    #[inline(always)]
    pub fn value(&self, i: usize) -> (u64, u64) {
        match self {
            AnyRange::Plain(r) => VD::get(r.positions[i].0),
            AnyRange::Packed(r) => r.get(i),
        }
    }

    /// Value `i` as a Flexmap cell (`VD::set` of it, as `VCell::set` stores it).
    #[inline(always)]
    pub fn cell(&self, i: usize) -> u64 {
        match self {
            AnyRange::Plain(r) => r.positions[i].0,
            AnyRange::Packed(r) => {
                let (value, pos) = r.get(i);
                VD::set(value, pos) & ((1 << 60) - 1)
            }
        }
    }

    /// Writes the flank distance of each value to `flex` into `out` (cleared first; nothing without
    /// a header) and returns the smallest (u32::MAX for none).
    #[inline]
    pub fn flank_distances(&self, flex: u32, out: &mut Vec<u8>) -> u32 {
        out.clear();
        let mut min = u32::MAX;
        match self {
            AnyRange::Plain(r) => {
                if let Some(h) = r.header {
                    out.reserve(h.len());
                    for cell in h {
                        let d = cell.dist(flex);
                        out.push(d as u8);
                        min = min.min(d);
                    }
                }
            }
            AnyRange::Packed(r) => {
                let cells = r.flex_bytes();
                out.reserve(cells.len() / 4);
                for c in cells.chunks_exact(4) {
                    let d = flank_distance(u32::from_le_bytes([c[0], c[1], c[2], c[3]]), flex);
                    out.push(d as u8);
                    min = min.min(d);
                }
            }
        }
        min
    }

    /// As `VRange::best_flex_match`.
    #[inline]
    pub fn best_flex_match<const F: usize, L>(&self, flex: &Kmer<F>, lambda: L)
    where
        L: FnMut(u64, u64, Option<(u32, u32)>),
    {
        match self {
            AnyRange::Plain(r) => r.best_flex_match(flex, lambda),
            AnyRange::Packed(r) => r.best_flex_match(flex, lambda),
        }
    }

    /// As `VRange::all_matches`: every value as (position, reference id).
    #[inline]
    pub fn all_matches<L>(&self, lambda: L)
    where
        L: FnMut(u64, u64),
    {
        match self {
            AnyRange::Plain(r) => r.all_matches(lambda),
            AnyRange::Packed(r) => r.all_matches(lambda),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::example::build_flexmap;
    use crate::flexmap::VRangeGetter;
    use crate::packed::{CompactKeys, PackedFlexmap};

    #[test]
    fn any_range_reads_both_indexes_alike() {
        let flexmap = build_flexmap();
        let packed = PackedFlexmap::<CompactKeys, 3, 8, 2>::from_flexmap(&flexmap);
        let (mut da, mut db) = (Vec::new(), Vec::new());
        for key in 0..64u64 {
            let (a, b) = (flexmap.get_vrange(key).map(AnyRange::from), packed.get(key).map(AnyRange::from));
            assert_eq!(a.is_some(), b.is_some());
            let (Some(a), Some(b)) = (a, b) else { continue };
            assert_eq!((a.len(), a.has_header()), (b.len(), b.has_header()));
            for i in 0..a.len() {
                assert_eq!((a.value(i), a.cell(i)), (b.value(i), b.cell(i)));
            }
            for q in [0u32, 0xFFFF, 0x1234] {
                assert_eq!(a.flank_distances(q, &mut da), b.flank_distances(q, &mut db));
                assert_eq!(da, db);
            }
        }
    }
}
