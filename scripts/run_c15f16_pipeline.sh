#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   scripts/run_c15f16_pipeline.sh [target_gb] [seed] [query_count] [load_repeats] [access_rounds] [log_file]
#
# Example:
#   scripts/run_c15f16_pipeline.sh 10 42 2000000 5 3 result/c15f16_pipeline.log

TARGET_GB="${1:-1}"
SEED="${2:-42}"
QUERY_COUNT="${3:-2000000}"
LOAD_REPEATS="${4:-5}"
ACCESS_ROUNDS="${5:-3}"
LOG_FILE="${6:-result/c15f16_pipeline.log}"
DATASET_DIR="result/c15f16"

mkdir -p "$(dirname "${LOG_FILE}")"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "===== C15/F16 pipeline started: $(date -Is) ====="
echo "log_file=${LOG_FILE}"
echo "target_gb=${TARGET_GB} seed=${SEED} query_count=${QUERY_COUNT} load_repeats=${LOAD_REPEATS} access_rounds=${ACCESS_ROUNDS}"

echo "[1/5] Building C15/F16 dataset in ${DATASET_DIR} (target_gb=${TARGET_GB}, seed=${SEED})"
cargo run --release -- build-default-dataset-c15f16 "${DATASET_DIR}" "${SEED}" "${TARGET_GB}"

echo "[2/5] Integrity test: regular vs blob (C15/F16)"
cargo test flexmap_blob_matches_c15_f16_default_dataset_files -- --nocapture

echo "[3/5] Access benchmark (get_vrange sampled)"
cargo run --release -- bench-access-c15f16 "${DATASET_DIR}" "${QUERY_COUNT}" "${ACCESS_ROUNDS}"

echo "[4/5] Cached load benchmark"
cargo run --release -- bench-load-c15f16 "${DATASET_DIR}" "${LOAD_REPEATS}"

echo "[5/5] Cleaning dataset"
rm -rf "${DATASET_DIR}"
echo "Done."
echo "===== C15/F16 pipeline finished: $(date -Is) ====="
