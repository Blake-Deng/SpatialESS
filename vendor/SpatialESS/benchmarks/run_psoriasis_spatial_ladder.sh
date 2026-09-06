#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R_BIN="${SPATIAL_V3_R_BIN:-/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript}"
OUT_ROOT="${SPATIALESS_PSORIASIS_OUT:-${PROJECT_DIR}/benchmarks/results/spatialcellchat_v3_psoriasis}"
PANEL="${SPATIAL_BENCH_PANEL:-small}"
NPERM="${SPATIAL_BENCH_NPERM:-20}"
SIZES="${SPATIAL_BENCH_SIZES:-100 250 500 905}"

export PATH="$(dirname "${R_BIN}"):/usr/bin:/bin"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export SPATIALESS_PSORIASIS_OUT="${OUT_ROOT}"
export SPATIAL_BENCH_PANEL="${PANEL}" SPATIAL_BENCH_NPERM="${NPERM}"

mkdir -p "${OUT_ROOT}"
"${R_BIN}" "${PROJECT_DIR}/benchmarks/prepare_psoriasis_spatial_v3_fixture.R"

for spots in ${SIZES}; do
  export SPATIAL_BENCH_SPOTS="${spots}"
  for method in spatialess spatialcellchat_v3; do
    export SPATIAL_BENCH_METHOD="${method}"
    run_dir="${OUT_ROOT}/$(printf 'spots%04d_%s_%s_perm%05d' "${spots}" "${PANEL}" "${method}" "${NPERM}")"
    mkdir -p "${run_dir}"
    /usr/bin/time -v -o "${run_dir}/time_verbose.txt" \
      "${R_BIN}" "${PROJECT_DIR}/benchmarks/run_psoriasis_spatial_headtohead.R" \
      > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"
  done
  "${R_BIN}" "${PROJECT_DIR}/benchmarks/compare_psoriasis_spatial_headtohead.R"
done

