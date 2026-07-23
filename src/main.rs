#![feature(exposed_provenance)]
use std::{cmp::min, collections::HashMap, env, fs::{self, File}, path::PathBuf, time::Instant};

use flexmap::{example::{build_flexmap, test_flexmap}, flexmap::Flexmap, keys::FMKeys};
use kmerrs::consecutive::kmer::{Kmer, KmerIter};
use bioreader::{fasta_byte_reader, fastq_byte_reader, fasta_reader, fastq_reader};
use savefile::{load, save};


fn test_simple() {
    let seq = "CATCGATCGTACGTGACTGCGTCGTCCTGCGTCGTCGTCGTGCTGCTGCTGTCGTCGTCGTCGTGCTGTCGTCGTA";
    const K: usize = 3;
    let kiter: KmerIter<3, true> = KmerIter::<K, true>::new(seq.as_bytes());
    let mut keys = FMKeys::<K, 8>::new();
    for (pos, kmer_fwd, kmer_rev) in kiter.clone() {
        let kmer = min(kmer_fwd, kmer_rev);
        keys.get_kmer_cell_mut_ref(kmer.0).increment();
    }

    let mut map = HashMap::<u64, u64>::new();
    for (pos, kmer_fwd, kmer_rev) in kiter.clone() {
        let kmer = min(kmer_fwd, kmer_rev);
        let entry = map.entry(kmer.0).or_insert_with(|| 0);
        *entry += 1;
    }

    for (pos, kmer_fwd, kmer_rev) in kiter {
        let kmer = min(kmer_fwd, kmer_rev);
        let cell = keys.get_kmer_cell(kmer.0);
        let ctrl_value = keys.get_control_head_value_from_kmer(kmer.0);
        println!("{} -> {} -> {} ({})", kmer.to_string().unwrap(), cell.0, map[&kmer.0], ctrl_value);
    }

    println!("{}", keys.get_control_head_value_from_kmer(0));
    keys.set_control_header_value_from_kmer(0, 42);
    println!("{}", keys.get_control_head_value_from_kmer(0));

    keys.build::<2>(100);

}



fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() >= 2 && args[1] == "build-test-data" {
        run_build_test_data(&args);
        return;
    }
    if args.len() >= 2 && args[1] == "build-default-dataset" {
        run_build_default_dataset(&args);
        return;
    }
    if args.len() >= 2 && args[1] == "build-default-dataset-c15f16" {
        run_build_default_dataset_c15f16(&args);
        return;
    }
    if args.len() >= 2 && args[1] == "bench-load-c15f16" {
        run_bench_load_c15f16(&args);
        return;
    }
    if args.len() >= 2 && args[1] == "bench-access-c15f16" {
        run_bench_access_c15f16(&args);
        return;
    }
    if args.len() >= 2 && args[1] == "bench-flex" {
        run_bench_flex(&args);
        return;
    }

    // println!("Main");
    // let testvec = vec![0u16, 1u16, 2u16, 3u16];
    // let ptr_to_u64: *const u64 = unsafe { transmute((&testvec[0] as *const u16).expose_provenance()) };
    // println!("{:?}", testvec);
    // let res = unsafe { *ptr_to_u64 };
    // dbg!(res);


    // let result: u64 = (testvec[0] as u64) | (testvec[1] as u64) << 16 | (testvec[2] as u64) << 32 | (testvec[3] as u64) << 48;
    // println!("res : {}", result);

    test_simple();

    
    let mut flexmap = build_flexmap();
    test_flexmap(&flexmap);

    
    const GLOBAL_VERSION: u32 = 1;
    let path = PathBuf::from("result/flexmap.save");
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };
    let _ = save(&mut file, GLOBAL_VERSION, &flexmap);

    let keys_path = "/usr/users/QIB_fr017/fritsche/ProjectsPrivate/flexalign/results/keys.bin".to_string();
    // unsafe { flexmap.keys.save_keys(&keys_path) };

    let mut file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };
    let flexmap: Flexmap<3, 8, 8, 2> = load(&mut file, GLOBAL_VERSION).expect("Loading did not work");

}

fn run_build_test_data(args: &[String]) {
    let out_dir = if args.len() >= 3 {
        PathBuf::from(&args[2])
    } else {
        PathBuf::from("result/testdata")
    };
    fs::create_dir_all(&out_dir).expect("failed to create output dir");

    let keys_path = out_dir.join("keys.bin");
    let values_path = out_dir.join("values.bin");
    let blob_path = out_dir.join("flexmap.blob");

    let flexmap = build_flexmap();
    flexmap.keys.save(&keys_path.to_string_lossy().to_string());
    flexmap.values.save(&values_path.to_string_lossy().to_string());
    flexmap.save_blob(&blob_path.to_string_lossy().to_string());

    let keys_abs = keys_path.canonicalize().expect("keys canonicalize failed");
    let values_abs = values_path.canonicalize().expect("values canonicalize failed");
    let blob_abs = blob_path.canonicalize().expect("blob canonicalize failed");

    println!("Built test data:");
    println!("  keys:   {}", keys_abs.display());
    println!("  values: {}", values_abs.display());
    println!("  blob:   {}", blob_abs.display());
    println!();
    println!("Use for large compare test:");
    println!("FLEXMAP15_KEYS={} \\", keys_abs.display());
    println!("FLEXMAP15_VALUES={} \\", values_abs.display());
    println!("FLEXMAP15_BLOB={} \\", blob_abs.display());
    println!("FLEXMAP15_QUERY_COUNT=200000 \\");
    println!("cargo test flexmap_blob_matches_flexmap_c15_f16_sampled_real_data -- --ignored --nocapture");
}

fn run_build_default_dataset(args: &[String]) {
    let out_dir = if args.len() >= 3 {
        PathBuf::from(&args[2])
    } else {
        PathBuf::from("result/default")
    };
    let seed = args.get(3).and_then(|v| v.parse::<u64>().ok()).unwrap_or(42);
    let target_gb = args.get(4).and_then(|v| v.parse::<u64>().ok()).unwrap_or(10).clamp(1, 10);
    let target_bytes = target_gb * 1024 * 1024 * 1024;
    // Empirical upper-bound estimate for this config (C=3,F=8,HT=2): ~10.67 bytes/base in values.
    // Round up so generated data tends to meet or exceed target size.
    let seq_len = ((target_bytes + 10) / 11) as usize;
    let max_range_size = args.get(5).and_then(|v| v.parse::<usize>().ok()).unwrap_or(seq_len);

    fs::create_dir_all(&out_dir).expect("failed to create output dir");
    let keys_path = out_dir.join("keys.bin");
    let values_path = out_dir.join("values.bin");
    let blob_path = out_dir.join("flexmap.blob");

    let flexmap = build_random_flexmap_c3f8(seed, seq_len, max_range_size);
    flexmap.keys.save(&keys_path.to_string_lossy().to_string());
    flexmap.values.save(&values_path.to_string_lossy().to_string());
    flexmap.save_blob(&blob_path.to_string_lossy().to_string());

    println!("Built default dataset from random DNA:");
    println!("  seed:   {seed}");
    println!("  target_gb: {target_gb}");
    println!("  target_bytes: {target_bytes}");
    println!("  length: {seq_len}");
    println!("  max_range_size: {max_range_size}");
    println!("  keys:   {}", keys_path.display());
    println!("  values: {}", values_path.display());
    println!("  blob:   {}", blob_path.display());
}

fn run_build_default_dataset_c15f16(args: &[String]) {
    let out_dir = if args.len() >= 3 {
        PathBuf::from(&args[2])
    } else {
        PathBuf::from("result/c15f16")
    };
    let seed = args.get(3).and_then(|v| v.parse::<u64>().ok()).unwrap_or(42);
    let target_gb = args.get(4).and_then(|v| v.parse::<u64>().ok()).unwrap_or(10).clamp(1, 10);
    let target_bytes = target_gb * 1024 * 1024 * 1024;
    let seq_len = ((target_bytes + 10) / 11) as usize;
    let max_range_size = args.get(5).and_then(|v| v.parse::<usize>().ok()).unwrap_or(seq_len);

    fs::create_dir_all(&out_dir).expect("failed to create output dir");
    let keys_path = out_dir.join("keys.bin");
    let values_path = out_dir.join("values.bin");
    let blob_path = out_dir.join("flexmap.blob");

    let flexmap = build_random_flexmap_c15f16(seed, seq_len, max_range_size);
    flexmap.keys.save(&keys_path.to_string_lossy().to_string());
    flexmap.values.save(&values_path.to_string_lossy().to_string());
    flexmap.save_blob(&blob_path.to_string_lossy().to_string());

    println!("Built C15/F16 dataset from random DNA:");
    println!("  seed:   {seed}");
    println!("  target_gb: {target_gb}");
    println!("  target_bytes: {target_bytes}");
    println!("  length: {seq_len}");
    println!("  max_range_size: {max_range_size}");
    println!("  keys:   {}", keys_path.display());
    println!("  values: {}", values_path.display());
    println!("  blob:   {}", blob_path.display());
}

fn run_bench_load_c15f16(args: &[String]) {
    let base = if args.len() >= 3 {
        PathBuf::from(&args[2])
    } else {
        PathBuf::from("result/c15f16")
    };
    let repeats = args.get(3).and_then(|v| v.parse::<usize>().ok()).unwrap_or(5);

    let keys = base.join("keys.bin").to_string_lossy().to_string();
    let values = base.join("values.bin").to_string_lossy().to_string();
    let blob = base.join("flexmap.blob").to_string_lossy().to_string();

    assert!(std::fs::metadata(&keys).is_ok(), "missing {}", keys);
    assert!(std::fs::metadata(&values).is_ok(), "missing {}", values);
    assert!(std::fs::metadata(&blob).is_ok(), "missing {}", blob);

    // Warm cache once.
    let warm_a = Flexmap::<15, 16, 16, 2>::load(&keys, &values);
    let warm_b = Flexmap::<15, 16, 16, 2>::load_blob(&blob);
    let _ = (warm_a.keys.data.len(), warm_b.keys().len());

    let mut regular_total = 0f64;
    let mut blob_total = 0f64;
    for _ in 0..repeats {
        let start = Instant::now();
        let map = Flexmap::<15, 16, 16, 2>::load(&keys, &values);
        let _ = map.keys.data.len();
        regular_total += start.elapsed().as_secs_f64();

        let start = Instant::now();
        let map = Flexmap::<15, 16, 16, 2>::load_blob(&blob);
        let _ = map.keys().len();
        blob_total += start.elapsed().as_secs_f64();
    }

    println!("Repeated load benchmark C15/F16 (cached):");
    println!("  repeats: {repeats}");
    println!("  regular avg: {:.6} s/load", regular_total / repeats as f64);
    println!("  blob avg:    {:.6} s/load", blob_total / repeats as f64);
}

fn run_bench_access_c15f16(args: &[String]) {
    let base = if args.len() >= 3 {
        PathBuf::from(&args[2])
    } else {
        PathBuf::from("result/c15f16")
    };
    let query_count = args.get(3).and_then(|v| v.parse::<usize>().ok()).unwrap_or(2_000_000);
    let rounds = args.get(4).and_then(|v| v.parse::<usize>().ok()).unwrap_or(3);

    let keys = base.join("keys.bin").to_string_lossy().to_string();
    let values = base.join("values.bin").to_string_lossy().to_string();
    let blob = base.join("flexmap.blob").to_string_lossy().to_string();

    assert!(std::fs::metadata(&keys).is_ok(), "missing {}", keys);
    assert!(std::fs::metadata(&values).is_ok(), "missing {}", values);
    assert!(std::fs::metadata(&blob).is_ok(), "missing {}", blob);

    let regular = Flexmap::<15, 16, 16, 2>::load(&keys, &values);
    let blob = Flexmap::<15, 16, 16, 2>::load_blob(&blob);
    let queries = build_query_keys_masked(query_count, (1u64 << (15 * 2)) - 1);

    let mut regular_total = 0f64;
    let mut blob_total = 0f64;
    let mut prefetch_total = 0f64;
    for _ in 0..rounds {
        let start = Instant::now();
        let n = scan_all(&regular, &queries);
        regular_total += start.elapsed().as_secs_f64();
        assert!(n > 0 || query_count == 0);

        let start = Instant::now();
        let n = scan_all(&blob, &queries);
        blob_total += start.elapsed().as_secs_f64();
        assert!(n > 0 || query_count == 0);

        let start = Instant::now();
        let n = scan_all_prefetched(&blob, &queries);
        prefetch_total += start.elapsed().as_secs_f64();
        assert!(n > 0 || query_count == 0);
    }

    let regular_avg = regular_total / rounds as f64;
    let blob_avg = blob_total / rounds as f64;
    let prefetch_avg = prefetch_total / rounds as f64;
    let mqs = |avg: f64| (query_count as f64 / avg) / 1e6;
    println!("Access benchmark C15/F16:");
    println!("  query_count: {query_count}");
    println!("  rounds:      {rounds}");
    println!("  regular avg:       {:.6} s/round ({:.2} Mq/s)", regular_avg, mqs(regular_avg));
    println!("  blob avg:          {:.6} s/round ({:.2} Mq/s)", blob_avg, mqs(blob_avg));
    println!("  blob+prefetch avg: {:.6} s/round ({:.2} Mq/s)", prefetch_avg, mqs(prefetch_avg));
    println!("  prefetch speedup vs blob:    {:.2}x", blob_avg / prefetch_avg);
    println!("  prefetch speedup vs regular: {:.2}x", regular_avg / prefetch_avg);
}

fn scan_all_prefetched<const C: usize, const F: usize, const CPB: u64, const HT: usize>(
    blob: &flexmap::flexmap::FlexmapBlob<C, F, CPB, HT>,
    queries: &[u64],
) -> usize {
    // FLEXMAP_PREFETCH_G overrides the rolling prefetch distance (default 24).
    // FLEXMAP_PREFETCH_BATCH chunks the query stream into calls of that many keys (0 = one call).
    let distance = env::var("FLEXMAP_PREFETCH_G").ok().and_then(|s| s.parse::<usize>().ok()).unwrap_or(24);
    let batch = env::var("FLEXMAP_PREFETCH_BATCH").ok().and_then(|s| s.parse::<usize>().ok()).unwrap_or(0);
    let mut total = 0usize;
    let mut acc = |_q: u64, vr: Option<flexmap::values::VRange>| {
        if let Some(v) = vr {
            total = total.wrapping_add(v.positions.len());
        }
    };
    if batch == 0 {
        blob.for_each_vrange_prefetched_rolling(queries, distance, &mut acc);
    } else {
        for chunk in queries.chunks(batch) {
            blob.for_each_vrange_prefetched_rolling(chunk, distance, &mut acc);
        }
    }
    total
}

fn run_bench_flex(args: &[String]) {
    use flexmap::values::{
        min_dist_and_count, min_dist_and_count_avx2, min_dist_and_count_avx512,
        min_dist_and_count_scalar, HeaderSeq,
    };

    let n = args.get(2).and_then(|v| v.parse::<usize>().ok()).unwrap_or(4096);
    let iters = args.get(3).and_then(|v| v.parse::<usize>().ok()).unwrap_or(200_000);

    // Build a header block of `n` pseudo-random packed 16-mers and a stream of query flanks.
    let mut s = 0x9E3779B97F4A7C15u64;
    let mut next = || {
        s ^= s >> 12;
        s ^= s << 25;
        s ^= s >> 27;
        s.wrapping_mul(0x2545F4914F6CDD1D)
    };
    let headers: Vec<HeaderSeq> = (0..n).map(|_| HeaderSeq::from_raw(next() as u32)).collect();
    let flexes: Vec<u32> = (0..iters).map(|_| next() as u32).collect();

    let avx2 = std::is_x86_feature_detected!("avx2");
    let avx512 = std::is_x86_feature_detected!("avx512f")
        && std::is_x86_feature_detected!("avx512bw")
        && std::is_x86_feature_detected!("avx512vpopcntdq");

    // Correctness cross-check before timing.
    for &f in flexes.iter().take(64) {
        let r = min_dist_and_count_scalar(&headers, f);
        if avx2 { assert_eq!(r, min_dist_and_count_avx2(&headers, f), "avx2 mismatch"); }
        if avx512 { assert_eq!(r, min_dist_and_count_avx512(&headers, f), "avx512 mismatch"); }
        assert_eq!(r, min_dist_and_count(&headers, f), "dispatch mismatch");
    }

    let time = |label: &str, f: &dyn Fn(&[HeaderSeq], u32) -> (u32, u32)| {
        // warm
        let mut acc = 0u64;
        for &fl in &flexes { acc = acc.wrapping_add(f(&headers, fl).0 as u64); }
        let mut best = f64::INFINITY;
        for _ in 0..5 {
            let start = Instant::now();
            let mut a = 0u64;
            for &fl in &flexes { a = a.wrapping_add(f(&headers, fl).0 as u64); }
            best = best.min(start.elapsed().as_secs_f64());
            acc = acc.wrapping_add(a);
        }
        let cmps = (n * iters) as f64; // header comparisons performed
        println!(
            "  {label:<8} {:.4} s  ({:.2} G cmp/s)   [checksum {}]",
            best, (cmps / best) / 1e9, acc & 0xffff
        );
    };

    println!("Flex-match reduction benchmark: n={n} headers, iters={iters} queries");
    println!("  CPU: avx2={avx2} avx512vpopcntdq={avx512}");
    time("scalar", &min_dist_and_count_scalar);
    #[cfg(target_arch = "x86_64")]
    {
        if avx2 { time("avx2", &min_dist_and_count_avx2); }
        if avx512 { time("avx512", &min_dist_and_count_avx512); }
    }
    time("dispatch", &min_dist_and_count);
}

fn scan_all<T: flexmap::flexmap::VRangeGetter>(map: &T, queries: &[u64]) -> usize {
    let mut total = 0usize;
    for &q in queries {
        if let Some(v) = map.get_vrange(q) {
            total = total.wrapping_add(v.positions.len());
        }
    }
    total
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

fn build_random_flexmap_c3f8(seed: u64, seq_len: usize, max_range_size: usize) -> Flexmap<3, 8, 8, 2> {
    let seq = random_dna(seed, seq_len);
    let kiter = KmerIter::<11, true>::new(&seq);
    let mut keys = FMKeys::<3, 8>::new();
    for (_, kmer_fwd, kmer_rev) in kiter.clone() {
        let kmer: Kmer<11> = if kmer_fwd < kmer_rev { kmer_fwd } else { kmer_rev };
        let cell = keys.get_kmer_cell_mut_ref(kmer.middle::<3>().0);
        assert!(
            cell.0 < u16::MAX,
            "default dataset overflow: key counter exceeded u16::MAX for C=3. \
             For very large datasets, use larger C (more key space) or smaller target size."
        );
        cell.increment();
    }
    keys.build::<2>(max_range_size);

    let mut flexmap = Flexmap::<3, 8, 8, 2>::new(keys);
    for (pos, kmer_fwd, kmer_rev) in kiter {
        let kmer: Kmer<11> = if kmer_fwd < kmer_rev { kmer_fwd } else { kmer_rev };
        let core = kmer.middle::<3>();
        let flanks = kmer.flanks::<8>();
        if let Some(range) = flexmap.keys.vrange(core.0) {
            let mut vblock = flexmap.values.get_range_mut(range);
            vblock.insert(flexmap::VD::set(1, pos as u64), flanks.0 as u32);
        }
    }
    flexmap
}

fn build_random_flexmap_c15f16(seed: u64, seq_len: usize, max_range_size: usize) -> Flexmap<15, 16, 16, 2> {
    let seq = random_dna(seed, seq_len);
    let kiter = KmerIter::<31, true>::new(&seq);
    let mut keys = FMKeys::<15, 16>::new();
    for (_, kmer_fwd, kmer_rev) in kiter.clone() {
        let kmer: Kmer<31> = if kmer_fwd < kmer_rev { kmer_fwd } else { kmer_rev };
        let cell = keys.get_kmer_cell_mut_ref(kmer.middle::<15>().0);
        assert!(
            cell.0 < u16::MAX,
            "default dataset overflow: key counter exceeded u16::MAX for C=15"
        );
        cell.increment();
    }
    keys.build::<2>(max_range_size);

    let mut flexmap = Flexmap::<15, 16, 16, 2>::new(keys);
    for (pos, kmer_fwd, kmer_rev) in kiter {
        let kmer: Kmer<31> = if kmer_fwd < kmer_rev { kmer_fwd } else { kmer_rev };
        let core = kmer.middle::<15>();
        let flanks = kmer.flanks::<16>();
        if let Some(range) = flexmap.keys.vrange(core.0) {
            let mut vblock = flexmap.values.get_range_mut(range);
            vblock.insert(flexmap::VD::set(1, pos as u64), flanks.0 as u32);
        }
    }
    flexmap
}

fn random_dna(seed: u64, len: usize) -> Vec<u8> {
    let mut x = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut seq = Vec::with_capacity(len);
    for _ in 0..len {
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        let n = (x.wrapping_mul(0x2545F4914F6CDD1D) & 0b11) as u8;
        let b = match n {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        seq.push(b);
    }
    seq
}
