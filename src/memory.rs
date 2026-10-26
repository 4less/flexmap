//! Memory for the large arrays of an index: transparent huge pages, and read-only byte sections
//! that share one backing (a heap buffer or a memory map), so that a loaded index is a set of views
//! into the bytes of its file.

use std::fs::File;
use std::os::unix::fs::FileExt;
use std::sync::{Arc, OnceLock};

/// Whether large allocations ask for transparent huge pages. On by default; `FLEXMAP_THP=0` turns it
/// off (to benchmark against plain 4 KB pages).
pub fn huge_pages_enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("FLEXMAP_THP").map(|v| v != "0").unwrap_or(true))
}

/// Asks the kernel for huge pages (`madvise(MADV_HUGEPAGE)`) on the 2 MB-aligned part of
/// `[ptr, ptr + len)`. A lookup reads the key table and the values at random, so with 4 KB pages
/// nearly every lookup also misses the TLB. Call it before the memory is first written, so that its
/// pages fault in as huge pages; Linux's default THP mode ("madvise") gives them only to memory that
/// asks. Does nothing for blocks under 64 MB, with `FLEXMAP_THP=0`, or off Linux.
pub fn advise_huge_pages(ptr: *const u8, len: usize) {
    if len < (64 << 20) || !huge_pages_enabled() {
        return;
    }
    #[cfg(target_os = "linux")]
    unsafe {
        const HUGE: usize = 1 << 21;
        let start = (ptr as usize + HUGE - 1) & !(HUGE - 1);
        let end = (ptr as usize + len) & !(HUGE - 1);
        if end > start {
            libc::madvise(start as *mut libc::c_void, end - start, libc::MADV_HUGEPAGE);
        }
    }
}

/// `n` zeroed words, not yet touched (large zeroed allocations are fresh pages) and advised for
/// huge pages.
pub fn zeroed_words(n: usize) -> Vec<u64> {
    let words = vec![0u64; n];
    advise_huge_pages(words.as_ptr() as *const u8, n * 8);
    words
}

/// `n` zeroed values of `T`, allocated zeroed (not written, unlike `vec![T::zeroed(); n]` for a
/// struct `T`) and advised for huge pages, so that the pages fault in as huge pages when filled.
pub fn zeroed_vec<T: bytemuck::Zeroable>(n: usize) -> Vec<T> {
    if n == 0 || std::mem::size_of::<T>() == 0 {
        return Vec::new();
    }
    let layout = std::alloc::Layout::array::<T>(n).expect("allocation too large");
    // Safety: T is Zeroable, so zeroed bytes are n valid values; the Vec frees with the same layout.
    unsafe {
        let ptr = std::alloc::alloc_zeroed(layout) as *mut T;
        if ptr.is_null() {
            std::alloc::handle_alloc_error(layout);
        }
        advise_huge_pages(ptr as *const u8, layout.size());
        Vec::from_raw_parts(ptr, n, n)
    }
}

/// Reads `buf.len()` bytes of `file` from `offset` on `threads` threads (`pread` in 32 MB pieces).
/// Each thread's pages fault in on that thread, so a large load is not bound by one core zeroing
/// and copying pages.
pub fn read_parallel(file: &File, offset: u64, buf: &mut [u8], threads: usize) -> std::io::Result<()> {
    const PIECE: usize = 32 << 20;
    let pieces = (buf.len() + PIECE - 1) / PIECE;
    let threads = threads.max(1).min(pieces.max(1));
    let next = std::sync::atomic::AtomicUsize::new(0);
    let base = buf.as_mut_ptr() as usize;
    let len = buf.len();
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                let next = &next;
                scope.spawn(move || -> std::io::Result<()> {
                    loop {
                        let piece = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if piece >= pieces {
                            return Ok(());
                        }
                        let from = piece * PIECE;
                        let to = (from + PIECE).min(len);
                        // Pieces are disjoint, so each thread writes its own bytes.
                        let part = unsafe { std::slice::from_raw_parts_mut((base + from) as *mut u8, to - from) };
                        file.read_exact_at(part, offset + from as u64)?;
                    }
                })
            })
            .collect();
        handles.into_iter().try_for_each(|h| h.join().expect("load thread panicked"))
    })
}

pub(crate) enum Backing {
    Words(Vec<u64>),
    Mapped(memmap2::Mmap),
}

impl Backing {
    fn base(&self) -> *const u8 {
        match self {
            Backing::Words(w) => w.as_ptr() as *const u8,
            Backing::Mapped(m) => m.as_ptr(),
        }
    }

    fn len(&self) -> usize {
        match self {
            Backing::Words(w) => w.len() * 8,
            Backing::Mapped(m) => m.len(),
        }
    }
}

/// Bytes that several structures share: a file (or part of one) read into a huge-page buffer, or
/// a whole file memory-mapped. [`Section`]s and a [`crate::packed::PackedFlexmap`] cut from it keep
/// it alive, so one read or one mapping can back an index and anything stored next to it in the
/// same file. Cloning shares the bytes.
#[derive(Clone)]
pub struct SharedBytes(Arc<Backing>, usize);

impl SharedBytes {
    pub(crate) fn from_arc(backing: Arc<Backing>) -> Self {
        let len = backing.len();
        Self(backing, len)
    }

    /// `len` bytes of `file` from `offset`, read on `threads` threads into huge pages
    /// ([`read_parallel`]). The buffer is 8-byte aligned and starts on a page.
    pub fn read(file: &File, offset: u64, len: usize, threads: usize) -> std::io::Result<Self> {
        let mut words = zeroed_words((len + 7) / 8);
        read_parallel(file, offset, &mut bytemuck::cast_slice_mut::<u64, u8>(&mut words)[..len], threads)?;
        Ok(Self(Arc::new(Backing::Words(words)), len))
    }

    /// The whole of `file`, memory-mapped: nothing is read until a page is touched.
    pub fn mmap(file: &File) -> std::io::Result<Self> {
        let map = unsafe { memmap2::Mmap::map(file)? };
        let len = map.len();
        Ok(Self(Arc::new(Backing::Mapped(map)), len))
    }

    /// Whether the bytes are a file mapping (pages can be dropped and re-read) rather than memory.
    pub fn is_mapped(&self) -> bool {
        matches!(*self.0, Backing::Mapped(_))
    }

    /// The mapping, when the bytes are one.
    pub fn mapping(&self) -> Option<&memmap2::Mmap> {
        match &*self.0 {
            Backing::Mapped(m) => Some(m),
            Backing::Words(_) => None,
        }
    }

    /// The bytes read or mapped (a read-in buffer is rounded up to whole words; this is not).
    pub fn len(&self) -> usize {
        self.1
    }

    pub fn bytes(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.0.base(), self.1) }
    }

    /// Bytes `[offset, offset + len)` as a section that keeps these bytes alive.
    pub fn section(&self, offset: usize, len: usize) -> Section {
        Section::view(&self.0, offset, len)
    }
}

/// A read-only run of bytes inside a shared [`Backing`]. Cloning shares the backing.
#[derive(Clone)]
pub struct Section {
    ptr: *const u8,
    len: usize,
    _backing: Arc<Backing>,
}

// The bytes are never written or moved once a Section exists, and the Arc keeps them alive.
unsafe impl Send for Section {}
unsafe impl Sync for Section {}

impl Section {
    /// The first `len` bytes of `words` (which must hold them) as a section of their own.
    pub fn from_words(words: Vec<u64>, len: usize) -> Self {
        assert!(len <= words.len() * 8);
        let backing = Arc::new(Backing::Words(words));
        Self { ptr: backing.base(), len, _backing: backing }
    }

    pub(crate) fn backing(backing: Backing) -> Arc<Backing> {
        Arc::new(backing)
    }

    /// Bytes `[offset, offset + len)` of `backing`.
    pub(crate) fn view(backing: &Arc<Backing>, offset: usize, len: usize) -> Self {
        assert!(offset + len <= backing.len(), "section out of the file's bounds");
        Self { ptr: unsafe { backing.base().add(offset) }, len, _backing: Arc::clone(backing) }
    }

    #[inline(always)]
    pub fn ptr(&self) -> *const u8 {
        self.ptr
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn bytes(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr, self.len) }
    }
}

/// Prefetches the cache line at `ptr` into all cache levels. Only a hint: results never change.
#[inline(always)]
pub fn prefetch(ptr: *const u8) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        core::arch::x86_64::_mm_prefetch(ptr as *const i8, core::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = ptr;
}
