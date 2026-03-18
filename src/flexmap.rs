use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};
use std::mem;
use std::slice;

use crate::keys::{FMKeys, FMKeysHash, KCell};
use crate::values::{FMValues, HeaderSeq, VCell, VRange};

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


pub trait FlexOptions {
    
}

pub trait VRangeGetter {
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange>;
}

pub trait DBBuilder {
    fn build(options: impl FlexOptions) -> Self;
}

const FLEXMAP_BLOB_MAGIC: [u8; 8] = *b"FMBLOB01";
const FLEXMAP_BLOB_VERSION: u32 = 1;
const FLEXMAP_BLOB_ALIGN: usize = 64;

#[derive(Clone, Copy)]
#[repr(C)]
struct FlexmapBlobHeader {
    magic: [u8; 8],
    version: u32,
    c: u32,
    f: u32,
    cells_per_body: u64,
    header_threshold: u64,
    keys_len: u64,
    values_len: u64,
    keys_offset: u64,
    values_offset: u64,
}

const FLEXMAP_BLOB_HEADER_SIZE: usize = 68;

fn align_up(value: usize, align: usize) -> usize {
    let mask = align - 1;
    (value + mask) & !mask
}

fn write_u32_le(dst: &mut [u8], offset: usize, value: u32) {
    dst[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn write_u64_le(dst: &mut [u8], offset: usize, value: u64) {
    dst[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

fn read_u32_le(src: &[u8], offset: usize) -> u32 {
    let mut bytes = [0u8; 4];
    bytes.copy_from_slice(&src[offset..offset + 4]);
    u32::from_le_bytes(bytes)
}

fn read_u64_le(src: &[u8], offset: usize) -> u64 {
    let mut bytes = [0u8; 8];
    bytes.copy_from_slice(&src[offset..offset + 8]);
    u64::from_le_bytes(bytes)
}

fn encode_header(header: FlexmapBlobHeader) -> [u8; FLEXMAP_BLOB_HEADER_SIZE] {
    let mut out = [0u8; FLEXMAP_BLOB_HEADER_SIZE];
    out[0..8].copy_from_slice(&header.magic);
    write_u32_le(&mut out, 8, header.version);
    write_u32_le(&mut out, 12, header.c);
    write_u32_le(&mut out, 16, header.f);
    write_u64_le(&mut out, 20, header.cells_per_body);
    write_u64_le(&mut out, 28, header.header_threshold);
    write_u64_le(&mut out, 36, header.keys_len);
    write_u64_le(&mut out, 44, header.values_len);
    write_u64_le(&mut out, 52, header.keys_offset);
    write_u64_le(&mut out, 60, header.values_offset);
    out
}

fn decode_header(src: &[u8]) -> FlexmapBlobHeader {
    let mut magic = [0u8; 8];
    magic.copy_from_slice(&src[0..8]);
    FlexmapBlobHeader {
        magic,
        version: read_u32_le(src, 8),
        c: read_u32_le(src, 12),
        f: read_u32_le(src, 16),
        cells_per_body: read_u64_le(src, 20),
        header_threshold: read_u64_le(src, 28),
        keys_len: read_u64_le(src, 36),
        values_len: read_u64_le(src, 44),
        keys_offset: read_u64_le(src, 52),
        values_offset: read_u64_le(src, 60),
    }
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

    pub fn save(&self, keys_file: &String, values_file: &String) -> () {
        self.keys.save(keys_file);
        self.values.save(values_file);
    }

    pub fn load(keys_file: &String, values_file: &String) -> Self {
        Self {
            keys: FMKeys::load(keys_file),
            values: FMValues::load(values_file),
        }
    }

    pub fn save_blob(&self, filename: &String) -> () {
        FlexmapBlob::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::save_to_file(self, filename)
    }

    pub fn load_blob(filename: &String) -> FlexmapBlob<C, F, CELLS_PER_BODY, HEADER_THRESHOLD> {
        FlexmapBlob::<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>::load_from_file(filename)
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

    pub fn save(&self, keys_file: &String, values_file: &String) -> () {
        self.keys.save(keys_file);
        self.values.save(values_file);
    }

    pub fn load(keys_file: &String, values_file: &String) -> Self {
        Self {
            keys: FMKeysHash::load(keys_file),
            values: FMValues::load(values_file),
        }
    }
}

pub struct FlexmapBlob<
    const C: usize,
    const F: usize,
    const CELLS_PER_BODY: u64,
    const HEADER_THRESHOLD: usize,
> {
    storage: Vec<u64>,
    file_size: usize,
    keys_ptr: *const KCell,
    keys_len: usize,
    values_ptr: *const VCell,
    values_len: usize,
}

// Safety: FlexmapBlob is immutable after construction.
// - `storage` owns the backing bytes and is moved with the struct.
// - `keys_ptr`/`values_ptr` are derived from `storage` and never mutated.
// - exposed APIs only hand out shared views (`&[KCell]`, `&[VCell]`, `VRange` with shared slices).
unsafe impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    Send for FlexmapBlob<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
}

// Safety rationale is identical to `Send`: concurrent shared access is read-only.
unsafe impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    Sync for FlexmapBlob<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    FlexmapBlob<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    pub fn save_to_file(
        flexmap: &Flexmap<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>,
        filename: &String,
    ) -> () {
        let mut file = File::create(filename).expect("no file found");

        let header_size = FLEXMAP_BLOB_HEADER_SIZE;
        let keys_offset = align_up(header_size, FLEXMAP_BLOB_ALIGN);
        let keys_bytes = bytemuck::cast_slice::<KCell, u8>(&flexmap.keys.data);
        let values_offset = align_up(keys_offset + keys_bytes.len(), FLEXMAP_BLOB_ALIGN);
        let values_bytes = bytemuck::cast_slice::<VCell, u8>(&flexmap.values.data);

        let header = FlexmapBlobHeader {
            magic: FLEXMAP_BLOB_MAGIC,
            version: FLEXMAP_BLOB_VERSION,
            c: C as u32,
            f: F as u32,
            cells_per_body: CELLS_PER_BODY,
            header_threshold: HEADER_THRESHOLD as u64,
            keys_len: flexmap.keys.data.len() as u64,
            values_len: flexmap.values.data.len() as u64,
            keys_offset: keys_offset as u64,
            values_offset: values_offset as u64,
        };

        let header_bytes = encode_header(header);
        file.write_all(&header_bytes).expect("write failed");

        if keys_offset > header_bytes.len() {
            let pad = vec![0u8; keys_offset - header_bytes.len()];
            file.write_all(&pad).expect("write failed");
        }

        file.write_all(keys_bytes).expect("write failed");

        let after_keys = keys_offset + keys_bytes.len();
        if values_offset > after_keys {
            let pad = vec![0u8; values_offset - after_keys];
            file.write_all(&pad).expect("write failed");
        }

        file.write_all(values_bytes).expect("write failed");
    }

    pub fn load_from_file(filename: &String) -> Self {
        let mut file = File::open(filename).expect("no file found");
        let file_size = file.metadata().expect("metadata failed").len() as usize;
        let words = (file_size + mem::size_of::<u64>() - 1) / mem::size_of::<u64>();
        let mut storage = vec![0u64; words];

        let bytes = unsafe { slice::from_raw_parts_mut(storage.as_mut_ptr() as *mut u8, file_size) };
        file.read_exact(bytes).expect("read failed");

        let header_size = FLEXMAP_BLOB_HEADER_SIZE;
        assert!(file_size >= header_size, "file too small for header");
        let header = decode_header(&bytes[..header_size]);

        assert!(header.magic == FLEXMAP_BLOB_MAGIC, "invalid blob magic");
        assert!(header.version == FLEXMAP_BLOB_VERSION, "unsupported blob version");
        assert!(header.c as usize == C, "blob C mismatch");
        assert!(header.f as usize == F, "blob F mismatch");
        assert!(header.cells_per_body == CELLS_PER_BODY, "blob cells_per_body mismatch");
        assert!(
            header.header_threshold as usize == HEADER_THRESHOLD,
            "blob header_threshold mismatch"
        );

        let keys_offset = header.keys_offset as usize;
        let keys_len = header.keys_len as usize;
        let keys_byte_len = keys_len * mem::size_of::<KCell>();
        let keys_end = keys_offset + keys_byte_len;

        let values_offset = header.values_offset as usize;
        let values_len = header.values_len as usize;
        let values_byte_len = values_len * mem::size_of::<VCell>();
        let values_end = values_offset + values_byte_len;

        assert!(keys_end <= file_size, "keys out of file bounds");
        assert!(values_end <= file_size, "values out of file bounds");
        assert!(keys_offset >= header_size, "invalid keys offset");
        assert!(values_offset >= keys_end, "invalid values offset");

        let keys_slice: &[KCell] =
            bytemuck::try_cast_slice(&bytes[keys_offset..keys_end]).expect("invalid keys section");
        let values_slice: &[VCell] = bytemuck::try_cast_slice(&bytes[values_offset..values_end])
            .expect("invalid values section");

        Self {
            storage,
            file_size,
            keys_ptr: keys_slice.as_ptr(),
            keys_len,
            values_ptr: values_slice.as_ptr(),
            values_len,
        }
    }

    pub fn keys(&self) -> &[KCell] {
        unsafe { slice::from_raw_parts(self.keys_ptr, self.keys_len) }
    }

    pub fn values(&self) -> &[VCell] {
        unsafe { slice::from_raw_parts(self.values_ptr, self.values_len) }
    }

    pub fn file_size(&self) -> usize {
        self.file_size
    }

    pub fn backing_bytes(&self) -> &[u8] {
        unsafe { slice::from_raw_parts(self.storage.as_ptr() as *const u8, self.file_size) }
    }

    #[inline(always)]
    unsafe fn get_control_header_value(keys_ptr: *const KCell, index: usize) -> u64 {
        *keys_ptr.add(index).cast::<u64>()
    }

    #[inline(always)]
    fn kmer_to_ctrl_block_index(canonical_kmer: u64) -> usize {
        let shift = CELLS_PER_BODY.ilog2() as u64;
        let cells_per_head = 4u64;
        ((canonical_kmer >> shift) * (cells_per_head + CELLS_PER_BODY)) as usize
    }

    #[inline(always)]
    fn vrange_from_keys(&self, canonical_kmer: u64) -> Option<(usize, usize)> {
        let cells_per_head = 4u64;
        let key_block_mask = CELLS_PER_BODY - 1;

        let block_index = Self::kmer_to_ctrl_block_index(canonical_kmer);
        let key_index = block_index + cells_per_head as usize + (canonical_kmer & key_block_mask) as usize;
        let keys_ptr = self.keys_ptr;
        let ctrl_block_value = unsafe { Self::get_control_header_value(keys_ptr, block_index) };
        let key_cell = unsafe { *keys_ptr.add(key_index) };

        let value_start = ctrl_block_value as usize + key_cell.0 as usize;
        let value_end = if (canonical_kmer & key_block_mask) != key_block_mask {
            let next_cell = unsafe { *keys_ptr.add(key_index + 1) };
            ctrl_block_value as usize + next_cell.0 as usize
        } else {
            let next_block = block_index + cells_per_head as usize + CELLS_PER_BODY as usize;
            unsafe { Self::get_control_header_value(keys_ptr, next_block) as usize }
        };

        if value_end <= value_start {
            return None;
        }
        Some((value_start, value_end))
    }

    #[inline(always)]
    fn get_range_from_values(&self, range: (usize, usize)) -> VRange {
        let (start, end) = range;
        let size = end - start;
        let values_ptr = self.values_ptr;

        if size > HEADER_THRESHOLD {
            let header_size = FMValues::<F, HEADER_THRESHOLD>::get_header_size(size);
            let values_size = size - header_size;
            let header_ptr = unsafe { values_ptr.add(start) };
            let positions_ptr = unsafe { values_ptr.add(start + header_size) };
            let header = unsafe { slice::from_raw_parts(header_ptr as *const HeaderSeq, values_size) };
            let positions = unsafe { slice::from_raw_parts(positions_ptr, size - header_size) };
            VRange::new(Some(header), positions)
        } else {
            let positions_ptr = unsafe { values_ptr.add(start) };
            let positions = unsafe { slice::from_raw_parts(positions_ptr, size) };
            VRange::new(None, positions)
        }
    }
}

impl<const C: usize, const F: usize, const CELLS_PER_BODY: u64, const HEADER_THRESHOLD: usize>
    VRangeGetter for FlexmapBlob<C, F, CELLS_PER_BODY, HEADER_THRESHOLD>
{
    fn get_vrange(&self, canonical_kmer: u64) -> Option<VRange> {
        let range = self.vrange_from_keys(canonical_kmer)?;
        Some(self.get_range_from_values(range))
    }
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

#[cfg(test)]
mod tests {
    use std::env;
    use std::fs;
    use std::cmp::min;
    use std::time::{SystemTime, UNIX_EPOCH};

    use kmerrs::consecutive::kmer::KmerIter;

    use crate::example::build_flexmap;
    use crate::flexmap::{Flexmap, FlexmapHash, VRangeGetter};
    use crate::keys::FMKeysHash;
    use crate::{keys::{FMKeys, KCell}, values::VCell, VD};
    use test::{black_box, Bencher};

    fn assert_vrange_equal(left: Option<crate::values::VRange>, right: Option<crate::values::VRange>) {
        match (left, right) {
            (None, None) => {}
            (Some(l), Some(r)) => {
                assert_eq!(l.positions.len(), r.positions.len(), "positions length differs");
                for idx in 0..l.positions.len() {
                    assert_eq!(l.positions[idx].0, r.positions[idx].0, "position cell differs at {idx}");
                }

                match (l.header, r.header) {
                    (None, None) => {}
                    (Some(lh), Some(rh)) => {
                        assert_eq!(lh.len(), rh.len(), "header length differs");
                        for idx in 0..lh.len() {
                            assert_eq!(lh[idx].get(), rh[idx].get(), "header cell differs at {idx}");
                        }
                    }
                    _ => panic!("header presence differs"),
                }
            }
            _ => panic!("one range is None while the other is Some"),
        }
    }

    #[test]
    fn flexmap_blob_matches_flexmap_queries() {
        let flexmap = build_flexmap();
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time drift")
            .as_nanos();
        let path = env::temp_dir().join(format!("flexmap_blob_compare_{}_{}.bin", std::process::id(), now));
        let filename = path.to_string_lossy().to_string();

        flexmap.save_blob(&filename);
        let blob = Flexmap::<3, 8, 8, 2>::load_blob(&filename);

        assert_eq!(blob.keys().len(), flexmap.keys.data.len(), "key length mismatch");
        assert_eq!(blob.values().len(), flexmap.values.data.len(), "value length mismatch");

        let ckmer_space = 1u64 << (3 * 2);
        for ckmer in 0..ckmer_space {
            let a = flexmap.get_vrange(ckmer);
            let b = blob.get_vrange(ckmer);
            assert_vrange_equal(a, b);
        }

        fs::remove_file(&path).expect("failed to cleanup blob test file");
    }

    fn build_quiet_flexmap() -> Flexmap<3, 8, 8, 2> {
        const C: usize = 3;
        const F: usize = 8;
        const K: usize = C + F;
        let seq = "CATCGATCGTACGTGACTGCGTCGTCCTGCGTCGTCGTCGTGCTGCTGCTGTCGTCGTCGTCGTGCTGTCGTCGTA";

        let kiter = KmerIter::<K, true>::new(seq.as_bytes());
        let mut keys = FMKeys::<C, 8>::new();
        for (_, kmer_fwd, kmer_rev) in kiter.clone() {
            let kmer = min(kmer_fwd, kmer_rev);
            keys.get_kmer_cell_mut_ref(kmer.middle::<C>().0).increment();
        }
        keys.build::<2>(1000);

        let mut flexmap = Flexmap::<C, F, 8, 2>::new(keys);
        for (pos, kmer_fwd, kmer_rev) in kiter {
            let kmer = min(kmer_fwd, kmer_rev);
            let core = kmer.middle::<C>();
            let flanks = kmer.flanks::<F>();
            if let Some(range) = flexmap.keys.vrange(core.0) {
                let mut vblock = flexmap.values.get_range_mut(range);
                vblock.insert(VD::set(1, pos as u64), flanks.0 as u32);
            }
        }

        flexmap
    }

    fn build_hash_from_flexmap(flexmap: &Flexmap<3, 8, 8, 2>) -> FlexmapHash<3, 8, 2> {
        let kmer_space = 1u64 << (3 * 2);
        let mut entries = Vec::<(u32, u64, u32)>::new();
        for ckmer in 0..kmer_space {
            if let Some((start, end)) = flexmap.keys.vrange(ckmer) {
                entries.push((ckmer as u32, start as u64, (end - start) as u32));
            }
        }

        let capacity = (entries.len() * 4).next_power_of_two().max(16);
        let mut keys = FMKeysHash::with_capacity(capacity);
        for (key, start, len) in entries {
            keys.insert(key, start, len).expect("hash insert failed");
        }

        FlexmapHash::<3, 8, 2> {
            keys,
            values: flexmap.values.clone(),
        }
    }

    fn build_query_keys(iterations: usize) -> Vec<u64> {
        build_query_keys_masked(iterations, (1u64 << (3 * 2)) - 1)
    }

    fn build_query_keys_masked(iterations: usize, mask: u64) -> Vec<u64> {
        let mut x = 0x9E3779B97F4A7C15u64;
        let mut out = Vec::with_capacity(iterations);
        for _ in 0..iterations {
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            let key = x.wrapping_mul(0x2545F4914F6CDD1D) & mask;
            out.push(key);
        }
        out
    }

    fn scan_all<T: VRangeGetter>(map: &T, queries: &[u64]) -> usize {
        let mut total = 0usize;
        for &q in queries {
            if let Some(v) = map.get_vrange(q) {
                total = total.wrapping_add(v.positions.len());
            }
        }
        total
    }

    #[bench]
    fn bench_get_vrange_flexmap(b: &mut Bencher) {
        let flexmap = build_quiet_flexmap();
        let queries = build_query_keys(1 << 20);

        b.iter(|| {
            let total = scan_all(&flexmap, &queries);
            black_box(total);
        });
    }

    #[bench]
    fn bench_get_vrange_flexmap_blob(b: &mut Bencher) {
        let flexmap = build_quiet_flexmap();
        let now = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time drift")
            .as_nanos();
        let path = env::temp_dir().join(format!("flexmap_blob_bench_{}_{}.bin", std::process::id(), now));
        let filename = path.to_string_lossy().to_string();
        flexmap.save_blob(&filename);
        let blob = Flexmap::<3, 8, 8, 2>::load_blob(&filename);
        fs::remove_file(&path).expect("failed to cleanup blob bench file");

        let queries = build_query_keys(1 << 20);

        b.iter(|| {
            let total = scan_all(&blob, &queries);
            black_box(total);
        });
    }

    #[bench]
    fn bench_get_vrange_flexmap_hash(b: &mut Bencher) {
        let flexmap = build_quiet_flexmap();
        let hash = build_hash_from_flexmap(&flexmap);
        let queries = build_query_keys(1 << 20);

        b.iter(|| {
            let total = scan_all(&hash, &queries);
            black_box(total);
        });
    }

    fn ensure_bench_load_dataset() -> (String, String, String) {
        let base = std::path::PathBuf::from("result/benchload");
        std::fs::create_dir_all(&base).expect("failed to create benchload dir");
        let keys = base.join("keys.bin").to_string_lossy().to_string();
        let values = base.join("values.bin").to_string_lossy().to_string();
        let blob = base.join("flexmap.blob").to_string_lossy().to_string();

        if std::fs::metadata(&keys).is_err() || std::fs::metadata(&values).is_err() || std::fs::metadata(&blob).is_err() {
            let map = build_quiet_flexmap();
            map.keys.save(&keys);
            map.values.save(&values);
            map.save_blob(&blob);
        }
        (keys, values, blob)
    }

    #[bench]
    fn bench_repeated_load_flexmap_cached(b: &mut Bencher) {
        let (keys, values, _) = ensure_bench_load_dataset();
        // Warm page cache once before benchmark loop.
        let warm = Flexmap::<3, 8, 8, 2>::load(&keys, &values);
        black_box((warm.keys.data.len(), warm.values.data.len()));

        b.iter(|| {
            let map = Flexmap::<3, 8, 8, 2>::load(&keys, &values);
            let digest = map.keys.data.len().wrapping_add(map.values.data.len());
            black_box(digest);
        });
    }

    #[bench]
    fn bench_repeated_load_flexmap_blob_cached(b: &mut Bencher) {
        let (_, _, blob) = ensure_bench_load_dataset();
        // Warm page cache once before benchmark loop.
        let warm = Flexmap::<3, 8, 8, 2>::load_blob(&blob);
        black_box((warm.keys().len(), warm.values().len()));

        b.iter(|| {
            let map = Flexmap::<3, 8, 8, 2>::load_blob(&blob);
            let digest = map.keys().len().wrapping_add(map.values().len());
            black_box(digest);
        });
    }

    fn digest_queries<T: VRangeGetter>(map: &T, queries: &[u64]) -> u128 {
        let mut acc = 0x9E3779B97F4A7C15u128;
        for &q in queries {
            match map.get_vrange(q) {
                None => {
                    acc = acc.rotate_left(7) ^ 0xA5A5A5A5A5A5A5A5u128;
                }
                Some(vr) => {
                    acc = acc.wrapping_mul(0x100000001B3u128) ^ (vr.positions.len() as u128);
                    for pos in vr.positions {
                        acc = acc.rotate_left(13) ^ (pos.0 as u128);
                    }
                    if let Some(header) = vr.header {
                        acc ^= header.len() as u128;
                        for h in header {
                            acc = acc.rotate_left(11) ^ (h.get() as u128);
                        }
                    }
                }
            }
        }
        acc
    }

    #[test]
    fn flexmap_blob_matches_default_dataset_files() {
        let keys_file = "result/default/keys.bin";
        let values_file = "result/default/values.bin";
        let blob_file = "result/default/flexmap.blob";

        if std::fs::metadata(keys_file).is_err()
            || std::fs::metadata(values_file).is_err()
            || std::fs::metadata(blob_file).is_err()
        {
            eprintln!("Skipping default-dataset file test: run `cargo run -- build-default-dataset` first");
            return;
        }

        let regular = Flexmap::<3, 8, 8, 2>::load(&keys_file.to_string(), &values_file.to_string());
        let blob = Flexmap::<3, 8, 8, 2>::load_blob(&blob_file.to_string());

        let queries = build_query_keys_masked(200_000, (1u64 << (3 * 2)) - 1);
        let regular_digest = digest_queries(&regular, &queries);
        let blob_digest = digest_queries(&blob, &queries);

        assert_eq!(
            regular_digest, blob_digest,
            "default dataset digest mismatch between Flexmap and FlexmapBlob"
        );
    }

    #[test]
    fn flexmap_blob_matches_c15_f16_default_dataset_files() {
        let keys_file = "result/c15f16/keys.bin";
        let values_file = "result/c15f16/values.bin";
        let blob_file = "result/c15f16/flexmap.blob";

        if std::fs::metadata(keys_file).is_err()
            || std::fs::metadata(values_file).is_err()
            || std::fs::metadata(blob_file).is_err()
        {
            eprintln!("Skipping C15/F16 file test: run `cargo run --release -- build-default-dataset-c15f16 result/c15f16 <seed> <target_gb>` first");
            return;
        }

        let regular = Flexmap::<15, 16, 16, 2>::load(&keys_file.to_string(), &values_file.to_string());
        let blob = Flexmap::<15, 16, 16, 2>::load_blob(&blob_file.to_string());

        let queries = build_query_keys_masked(200_000, (1u64 << (15 * 2)) - 1);
        let regular_digest = digest_queries(&regular, &queries);
        let blob_digest = digest_queries(&blob, &queries);

        assert_eq!(
            regular_digest, blob_digest,
            "C15/F16 dataset digest mismatch between Flexmap and FlexmapBlob"
        );
    }

    #[test]
    #[ignore = "Large-data validation. Uses result/c15f16/{keys.bin,values.bin,flexmap.blob} by default; env vars optional"]
    fn flexmap_blob_matches_flexmap_c15_f16_sampled_real_data() {
        fn env_or_default_existing_path(var: &str, default_path: &str) -> Option<String> {
            let path = env::var(var).unwrap_or_else(|_| default_path.to_string());
            match std::fs::metadata(&path) {
                Ok(_) => Some(path),
                Err(e) => {
                    eprintln!(
                        "Skipping large-data test: {var} not set and default missing/unreadable.\n  expected path: {path}\n  error: {e}"
                    );
                    None
                }
            }
        }

        let Some(keys_file) = env_or_default_existing_path("FLEXMAP15_KEYS", "result/c15f16/keys.bin") else { return; };
        let Some(values_file) = env_or_default_existing_path("FLEXMAP15_VALUES", "result/c15f16/values.bin") else { return; };
        let Some(blob_file) = env_or_default_existing_path("FLEXMAP15_BLOB", "result/c15f16/flexmap.blob") else { return; };

        let query_count = env::var("FLEXMAP15_QUERY_COUNT")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(200_000);
        let key_mask = (1u64 << (15 * 2)) - 1;
        let queries = build_query_keys_masked(query_count, key_mask);

        eprintln!("Using FLEXMAP15_KEYS={keys_file}");
        eprintln!("Using FLEXMAP15_VALUES={values_file}");
        eprintln!("Using FLEXMAP15_BLOB={blob_file}");
        eprintln!("Using FLEXMAP15_QUERY_COUNT={query_count}");

        // Preflight: keys file for C=15/CPB=16 has a fixed size. Mismatches can cause UB on raw loads.
        let keys_file_size = std::fs::metadata(&keys_file).expect("keys metadata").len();
        let expected_keys_cells = FMKeys::<15, 16>::table_size();
        let expected_keys_size = expected_keys_cells
            .checked_mul(std::mem::size_of::<KCell>() as u64)
            .expect("expected keys size overflow");
        if keys_file_size != expected_keys_size {
            eprintln!(
                "Skipping large-data test: FLEXMAP15_KEYS has wrong size for C=15/CELLS_PER_BODY=16 (got {keys_file_size}, expected {expected_keys_size})"
            );
            return;
        }

        let values_file_size = std::fs::metadata(&values_file).expect("values metadata").len();
        if values_file_size % std::mem::size_of::<VCell>() as u64 != 0 {
            eprintln!(
                "Skipping large-data test: FLEXMAP15_VALUES size must be a multiple of {} bytes",
                std::mem::size_of::<VCell>()
            );
            return;
        }

        // Validate blob format/params first. This fails safely on mismatched C/F.
        let blob = Flexmap::<15, 16, 16, 2>::load_blob(&blob_file);

        // Load regular map first, compute digest, then drop to avoid doubling peak RAM.
        let regular = Flexmap::<15, 16, 16, 2>::load(&keys_file, &values_file);
        let regular_digest = digest_queries(&regular, &queries);
        drop(regular);

        let blob_digest = digest_queries(&blob, &queries);

        assert_eq!(
            regular_digest, blob_digest,
            "sampled digest mismatch between Flexmap and FlexmapBlob"
        );
    }
}
