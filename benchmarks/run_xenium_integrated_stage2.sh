#!/usr/bin/env bash
set -u

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R_BIN="${SPATIALESS_R_BIN:-/home/dzf/miniforge3/envs/cellchat-acceleration/bin/Rscript}"
OUT_ROOT="${XENIUM_INTEGRATED_OUT:-/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/spatialess_integrated_stage2}"
NPERM="${XENIUM_INTEGRATED_NPERM:-20}"
RUN_DIR="${OUT_ROOT}/cells1156091_grid50_permutation_$(printf 'perm%05d' "${NPERM}")"

export PATH="$(dirname "${R_BIN}"):/usr/bin:/bin"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p "${RUN_DIR}"
/usr/bin/time -v -o "${RUN_DIR}/time_verbose.txt" \
  "${R_BIN}" "${PROJECT_DIR}/benchmarks/run_xenium_integrated_stage2.R" \
  > "${RUN_DIR}/stdout.log" 2> "${RUN_DIR}/stderr.log"
status=$?
printf '%s\n' "${status}" > "${RUN_DIR}/exit_status.txt"
exit "${status}"
