# Packed flexmap against the Flexmap blob: memory and lookup speed

2026-10-23. Branch `packed-index` (from `origin/zeroable_load` at `cebfd56`, the rev flexalign
pins). Quick runs: 2 M queries per dataset and scenario, one timed round after a warm-up. The full
page with every variant and scenario is `report.html` (made by `scripts/packbench_report.py` from
`results/`).

## Result

Short reads, one lookup call per 150 bp read (as an aligner takes a read at a time), both sides
with their batch prefetch, M seed lookups per second:

| Dataset | Old index | Packed compact | Packed hashed | Old 1t | Compact 1t | Hashed 1t | Old 8t | Compact 8t | Hashed 8t |
|---|---|---|---|---|---|---|---|---|---|
| D1, 450 mixed bacterial genomes (2.05 Gbp) | 7.58 GB | 4.12 GB | 4.35 GB | 1.57 | 3.18 | 3.37 | 9.39 | 16.46 | 18.04 |
| D2, 300 *B. uniformis* genomes (1.31 Gbp) | 5.74 GB | 2.50 GB | – | 0.84 | 1.14 | – | 4.36 | 5.61 | – |
| D3, human GRCh37 chr1–8 (1.54 Gbp) | 5.76 GB | 2.39 GB | 2.61 GB | 1.88 | 3.37 | 3.51 | 11.15 | 18.59 | 18.36 |

- Every variant emits the same seeds: the hash of every (position, reference, distance, ties) per
  query agrees across all variants, and the unit tests compare every key of a test index field by
  field.
- "Old" is the blob as flexalign loads it (`load_blob`, no huge pages). Huge pages alone take the
  old blob to 2.18 / 0.96 / 2.74 M/s (1 thread, D1 / D2 / D3).
- D2 is bound by emitting ~150 seeds a lookup, which every variant pays alike; the one-pass flank
  scan, more value lines prefetched and other prefetch distances did not help there
  (`results/strategies_d2.tsv`).
- Reading every value of a key is slower packed on D2 and D3 (decoding the bit fields); random keys
  (almost all misses) are slower with the hashed key table.

## What changed

See the "What was changed, and what was tried" section of `report.html`: packed values
(`src/packed.rs`, `PackedValues`), compact and hashed key tables (`CompactKeys`, `HashedKeys`), a
prefetch pipeline (`for_each_prefetched`), huge pages and a parallel load (`src/memory.rs`), two
fixes in `FMKeys::build` and `KCell::increment`.

## Data and commands

- D1: 450 random genomes of `~/.gtdb_downloader/raw` (`results/d1_list.txt`); D2: 300 genomes of
  `.../r232/.../Bacteroides_uniformis` (`results/d2_list.txt`); D3: `NC_000001`–`NC_000008` of
  `~/Databases/genomes/GRCh37/GCF_000001405.13_GRCh37_genomic.fna`. Made with
  `scripts/prep_data.sh <dir>` (uppercase, headers `g<genome>_c<contig>`).
- Index: `packbench build <fasta> <prefix> 2048 old,compact[,hashed]` (K 31, C 15, F 16, S 4,
  16 keys per block, threshold 2); `packbench hash <compact.pk> <hashed.pk>`.
- Queries: `packbench reads <fasta> <out.q> 2000000 0.01 7 150` and
  `packbench reads <fasta> <out_long.q> 2000 0.01 11 10000`.
- Lookups: `FLEXMAP_THP=0|1 FLEXMAP_QUERIES=2000000 FLEXMAP_LONG=<long.q> packbench query
  <old|compact|hashed> heap <index> <short.q> 1,8 1` (`scripts/run_queries.sh`).
- Machine: Intel Core i9-13900HX (AVX2), 12 logical CPUs, 31 GB RAM shared with other work, Linux
  6.17, rustc nightly 1.99. Indexes on /dev/shm or SSD. Some loads ran under memory pressure, so
  load times vary; the D1 baseline, which ran while the machine swapped, was repeated.
