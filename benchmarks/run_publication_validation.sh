#!/bin/sh
set -eu

PROJECT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
if [ -z "$(printenv R_BIN 2>/dev/null || true)" ]; then
  R_BIN=/home/dzf/miniforge3/envs/cellchat-acceleration/bin/Rscript
else
  R_BIN=$(printenv R_BIN)
fi
if [ -z "$(printenv SPATIALESS_PUBLICATION_OUT 2>/dev/null || true)" ]; then
  OUT_ROOT=/data/dzf/SecAct_2026/results/spatialess_publication_validation
else
  OUT_ROOT=$(printenv SPATIALESS_PUBLICATION_OUT)
fi
RUNNER="$PROJECT_DIR/benchmarks/run_publication_validation.R"
mkdir -p "$OUT_ROOT"

run_one() {
  experiment=$1
  nperm=$2
  seed=$3
  contact=$4
  diffusion=$5
  replicate=$6
  tag=$(printf '%s_cells%06d_perm%05d_seed%010d_contact%03g_diffusion%03g_block%03g_rep%02d' \
    "$experiment" 5000 "$nperm" "$seed" "$contact" "$diffusion" 100 "$replicate")
  run_dir="$OUT_ROOT/$tag"
  mkdir -p "$run_dir"
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
    SPATIALESS_PUBLICATION_OUT="$OUT_ROOT" \
    SPATIALESS_EXPERIMENT="$experiment" SPATIALESS_CELLS=5000 \
    SPATIALESS_NPERM="$nperm" SPATIALESS_SEED="$seed" \
    SPATIALESS_CONTACT_UM="$contact" SPATIALESS_DIFFUSION_UM="$diffusion" \
    SPATIALESS_BLOCK_UM=100 SPATIALESS_REPLICATE="$replicate" \
    /usr/bin/time -v -o "$run_dir/time_verbose.txt" \
    "$R_BIN" "$RUNNER" > "$run_dir/stdout.log" 2> "$run_dir/stderr.log"
  printf '0\n' > "$run_dir/exit_status.txt"
}

for replicate in 1 2 3 4 5; do
  run_one reproducibility 20 20260804 15 35 "$replicate"
done
for nperm in 20 50 100 200; do
  run_one permutation "$nperm" 20260804 15 35 1
done
for seed in 20260804 20260805 20260806 20260807 20260808; do
  run_one seed_stability 100 "$seed" 15 35 1
done
for diffusion in 20 35 50 75 100; do
  run_one radius 20 20260804 15 "$diffusion" 1
done

"$R_BIN" "$PROJECT_DIR/benchmarks/summarize_publication_validation.R"
