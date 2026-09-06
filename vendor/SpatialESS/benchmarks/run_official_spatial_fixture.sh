#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${PROJECT_DIR}/benchmarks/results/official_spatial_fixture"
export PATH="/home/dzf/miniforge3/envs/cellchat-acceleration/bin:${PATH}"
mkdir -p "${OUT_DIR}"

Rscript "${PROJECT_DIR}/benchmarks/prepare_official_spatial_fixture.R"
/usr/bin/time -v -o "${OUT_DIR}/official_time_verbose.txt" \
  Rscript "${PROJECT_DIR}/benchmarks/run_official_cellchat_spatial_fixture.R" \
  >"${OUT_DIR}/official_run.log" 2>&1
/usr/bin/time -v -o "${OUT_DIR}/spatialess_time_verbose.txt" \
  Rscript "${PROJECT_DIR}/benchmarks/run_spatialess_official_fixture.R" \
  >"${OUT_DIR}/spatialess_run.log" 2>&1
Rscript "${PROJECT_DIR}/benchmarks/compare_official_cellchat_spatial_fixture.R" \
  >"${OUT_DIR}/compare.log" 2>&1

cat "${OUT_DIR}/official_summary.tsv"
cat "${OUT_DIR}/spatialess_summary.tsv"
cat "${OUT_DIR}/comparison_summary.tsv"

