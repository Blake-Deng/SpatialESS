#!/usr/bin/env bash
set -u

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R_BIN="${SPATIALESS_R_BIN:-/home/dzf/miniforge3/envs/cellchat-acceleration/bin/Rscript}"
OUT_ROOT="${SECACT_INTEGRATED_OUT:-/data/dzf/SecAct_2026/results/spatialess_integrated_stage2}"
CELLS="${SECACT_INTEGRATED_CELLS:-1000}"
NPERM="${SECACT_INTEGRATED_NPERM:-20}"
INFERENCE="${SECACT_INTEGRATED_INFERENCE:-permutation}"
RUN_DIR="${OUT_ROOT}/$(printf 'cells%06d_full_%s_perm%05d' "${CELLS}" "${INFERENCE}" "${NPERM}")"

export PATH="$(dirname "${R_BIN}"):/usr/bin:/bin"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p "${RUN_DIR}"
/usr/bin/time -v -o "${RUN_DIR}/time_verbose.txt" \
  "${R_BIN}" "${PROJECT_DIR}/benchmarks/run_secact_cosmx_integrated.R" \
  > "${RUN_DIR}/stdout.log" 2> "${RUN_DIR}/stderr.log"
