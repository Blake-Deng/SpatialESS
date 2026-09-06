#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R_BIN="${SPATIAL_V3_R_BIN:-/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript}"
OUT_ROOT="${SECACT_V3_OUT:-/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead}"
CELLS="${SECACT_V3_CELLS:-250}"
PANEL="${SECACT_V3_PANEL:-small}"
NPERM="${SECACT_V3_NPERM:-20}"

export PATH="$(dirname "${R_BIN}"):/usr/bin:/bin"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export SECACT_V3_OUT="${OUT_ROOT}" SECACT_V3_CELLS="${CELLS}"
export SECACT_V3_PANEL="${PANEL}" SECACT_V3_NPERM="${NPERM}"
mkdir -p "${OUT_ROOT}"

if [[ ! -f "${OUT_ROOT}/secact_cosmx_v3_common_fixture.rds" ]]; then
  /usr/bin/time -v -o "${OUT_ROOT}/prepare_fixture.time_verbose.txt" \
    "${R_BIN}" "${PROJECT_DIR}/benchmarks/prepare_secact_cosmx_v3_fixture.R" \
    > "${OUT_ROOT}/prepare_fixture.stdout.log" \
    2> "${OUT_ROOT}/prepare_fixture.stderr.log"
fi

for method in spatialess spatialcellchat_v3; do
  export SECACT_V3_METHOD="${method}"
  run_dir="${OUT_ROOT}/$(printf 'cells%06d_%s_%s_perm%05d' "${CELLS}" "${PANEL}" "${method}" "${NPERM}")"
  mkdir -p "${run_dir}"
  /usr/bin/time -v -o "${run_dir}/time_verbose.txt" \
    "${R_BIN}" "${PROJECT_DIR}/benchmarks/run_secact_cosmx_v3_headtohead.R" \
    > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"
done

"${R_BIN}" "${PROJECT_DIR}/benchmarks/compare_secact_cosmx_v3_headtohead.R"

