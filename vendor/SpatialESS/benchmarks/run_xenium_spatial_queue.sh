#!/usr/bin/env bash
set -u

ROOT=/home/dzf/cellchat_acceleration/SpatialESS
OUT="$ROOT/benchmarks/results/xenium_spatial_new_20260731"
ESS_R=/home/dzf/miniforge3/envs/cellchat-acceleration/bin/Rscript
V3_R=/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript
mkdir -p "$OUT"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export XENIUM_SPATIAL_NPERM=20
export XENIUM_SPATIAL_RADIUS_UM=30
export XENIUM_SPATIAL_BLOCK_UM=100
export XENIUM_SPATIAL_SCALE_DISTANCE=20

stage() {
  printf '[%s] %s\n' "$(date -Is)" "$1" >> "$OUT/queue.log"
}

run_ess() {
  base_label="$1"
  n_cells="$2"
  panel="$3"
  export XENIUM_SPATIAL_NCELLS="$n_cells"
  export XENIUM_SPATIAL_PANEL="$panel"
  label="${base_label}_${panel}"
  stage "START ESS ${label}"
  /usr/bin/time -v -o "$OUT/${label}_spatialess.time_verbose.txt" \
    "$ESS_R" "$ROOT/benchmarks/run_xenium_spatialess_new.R" \
    >"$OUT/${label}_spatialess.stdout.log" \
    2>"$OUT/${label}_spatialess.stderr.log"
  rc=$?
  stage "END ESS ${label} rc=${rc}"
  return 0
}

run_v3() {
  base_label="$1"
  n_cells="$2"
  panel="$3"
  export XENIUM_SPATIAL_NCELLS="$n_cells"
  export XENIUM_SPATIAL_PANEL="$panel"
  label="${base_label}_${panel}"
  stage "START V3 ${label}"
  /usr/bin/time -v -o "$OUT/${label}_spatialcellchat_v3.time_verbose.txt" \
    "$V3_R" "$ROOT/benchmarks/run_xenium_spatialcellchat_v3_new.R" \
    >"$OUT/${label}_spatialcellchat_v3.stdout.log" \
    2>"$OUT/${label}_spatialcellchat_v3.stderr.log"
  rc=$?
  stage "END V3 ${label} rc=${rc}"
  return 0
}

compare() {
  label="$1"
  stage "START COMPARE ${label}"
  XENIUM_SPATIAL_COMPARE_LABEL="$label" "$ESS_R" \
    "$ROOT/benchmarks/compare_xenium_spatial_new.R" \
    >"$OUT/${label}_comparison.log" 2>&1 || true
  stage "END COMPARE ${label}"
}

stage "QUEUE START"
run_ess n000250 250 full
run_v3 n000250 250 full
compare n000250_full

run_ess n005000 5000 small
run_v3 n005000 5000 small
compare n005000_small

run_ess full 1156091 full

stage "START V3 full full guarded at 120 GiB virtual memory"
export XENIUM_SPATIAL_NCELLS=1156091
export XENIUM_SPATIAL_PANEL=full
ulimit -v $((120 * 1024 * 1024))
/usr/bin/time -v -o "$OUT/full_full_spatialcellchat_v3.time_verbose.txt" \
  "$V3_R" "$ROOT/benchmarks/run_xenium_spatialcellchat_v3_new.R" \
  >"$OUT/full_full_spatialcellchat_v3.stdout.log" \
  2>"$OUT/full_full_spatialcellchat_v3.stderr.log"
stage "END V3 full full rc=$?"
stage "QUEUE FINISHED"
