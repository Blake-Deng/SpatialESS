#!/usr/bin/env bash
set -u

ROOT=/home/dzf/cellchat_acceleration/SpatialESS
OUT="$ROOT/benchmarks/results/xenium_spatial_new_20260731"
ESS_R=/home/dzf/miniforge3/envs/cellchat-acceleration/bin/Rscript
V3_R=/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript
mkdir -p "$OUT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLAS_NUM_THREADS=1
export XENIUM_SPATIAL_NPERM=20 XENIUM_SPATIAL_SCALE_DISTANCE=20 XENIUM_SPATIAL_RADIUS_UM=30 XENIUM_SPATIAL_BLOCK_UM=100

run_ess() {
  local label="$1"
  export XENIUM_SPATIAL_NCELLS="$2"
  /usr/bin/time -v -o "$OUT/${label}_spatialess.time_verbose.txt" \
    "$ESS_R" "$ROOT/benchmarks/run_xenium_spatialess_new.R" \
    >"$OUT/${label}_spatialess.stdout.log" 2>"$OUT/${label}_spatialess.stderr.log"
}

run_v3() {
  local label="$1"
  export XENIUM_SPATIAL_NCELLS="$2"
  "$V3_R" "$ROOT/benchmarks/run_xenium_spatialcellchat_v3_new.R" \
    >"$OUT/${label}_spatialcellchat_v3.stdout.log" \
    2>"$OUT/${label}_spatialcellchat_v3.stderr.log"
}

echo "[$(date -Is)] ESS subset start" >> "$OUT/driver.log"
run_ess n005000 5000
echo "[$(date -Is)] SpatialCellChat v3 subset start" >> "$OUT/driver.log"
run_v3 n005000 5000
echo "[$(date -Is)] subset comparison" >> "$OUT/driver.log"
XENIUM_SPATIAL_COMPARE_LABEL=n005000 "$ESS_R" \
  "$ROOT/benchmarks/compare_xenium_spatial_new.R" \
  >"$OUT/n005000_comparison.log" 2>&1 || true

echo "[$(date -Is)] ESS full start" >> "$OUT/driver.log"
run_ess full 1156091

echo "[$(date -Is)] SpatialCellChat v3 full guarded start" >> "$OUT/driver.log"
export XENIUM_SPATIAL_NCELLS=1156091
ulimit -v $((120 * 1024 * 1024))
/usr/bin/time -v -o "$OUT/full_spatialcellchat_v3.time_verbose.txt" \
  "$V3_R" "$ROOT/benchmarks/run_xenium_spatialcellchat_v3_new.R" \
  >"$OUT/full_spatialcellchat_v3.stdout.log" \
  2>"$OUT/full_spatialcellchat_v3.stderr.log" || true

echo "[$(date -Is)] all stages finished" >> "$OUT/driver.log"
