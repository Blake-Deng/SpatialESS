#!/usr/bin/env bash
set -u

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R_BIN="${SPATIAL_V3_R_BIN:-/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript}"
OUT_ROOT="${SECACT_V3_OUT:-/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead}"
SIZES="${SECACT_V3_SIZES:-250 1000 5000 10000 25000 50000 100000 200000 443515}"
NPERM="${SECACT_V3_NPERM:-20}"
V3_MEMORY_GIB="${SECACT_V3_MEMORY_GIB:-120}"
V3_MEMORY_KB=$((V3_MEMORY_GIB * 1024 * 1024))

export PATH="$(dirname "${R_BIN}"):/usr/bin:/bin"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export SECACT_V3_OUT="${OUT_ROOT}" SECACT_V3_NPERM="${NPERM}"
export SECACT_V3_SIZES="${SIZES}"
mkdir -p "${OUT_ROOT}"
QUEUE_LOG="${OUT_ROOT}/ladder_queue.log"
STATUS_TSV="${OUT_ROOT}/ladder_queue_status.tsv"

if [[ ! -f "${STATUS_TSV}" ]]; then
  printf 'timestamp\tcells\tpanel\tmethod\tnperm\tevent\texit_status\tmemory_limit_gib\trun_directory\n' > "${STATUS_TSV}"
fi

stage() {
  printf '[%s] %s\n' "$(date -Is)" "$1" >> "${QUEUE_LOG}"
}

record_status() {
  cells="$1"
  method="$2"
  event="$3"
  rc="$4"
  run_dir="$5"
  printf '%s\t%s\tfull\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date -Is)" "${cells}" "${method}" "${NPERM}" "${event}" "${rc}" \
    "${V3_MEMORY_GIB}" "${run_dir}" >> "${STATUS_TSV}"
}

is_complete() {
  run_dir="$1"
  [[ -s "${run_dir}/summary.tsv" ]] &&
    [[ -s "${run_dir}/time_verbose.txt" ]] &&
    grep -q 'Exit status: 0' "${run_dir}/time_verbose.txt"
}

refresh_summary() {
  "${R_BIN}" "${PROJECT_DIR}/benchmarks/summarize_secact_cosmx_v3_ladder.R" \
    > "${OUT_ROOT}/ladder_summary_latest.log" 2>&1 || true
}

run_method() {
  cells="$1"
  method="$2"
  export SECACT_V3_CELLS="${cells}"
  export SECACT_V3_PANEL=full
  export SECACT_V3_METHOD="${method}"
  run_dir="${OUT_ROOT}/$(printf 'cells%06d_full_%s_perm%05d' "${cells}" "${method}" "${NPERM}")"
  mkdir -p "${run_dir}"
  if is_complete "${run_dir}"; then
    stage "SKIP completed cells=${cells} method=${method}"
    record_status "${cells}" "${method}" "skip_completed" 0 "${run_dir}"
    refresh_summary
    return 0
  fi

  stage "START cells=${cells} method=${method}"
  record_status "${cells}" "${method}" "start" -1 "${run_dir}"
  if [[ "${method}" == "spatialcellchat_v3" ]]; then
    (
      ulimit -v "${V3_MEMORY_KB}"
      /usr/bin/time -v -o "${run_dir}/time_verbose.txt" \
        "${R_BIN}" "${PROJECT_DIR}/benchmarks/run_secact_cosmx_v3_headtohead.R" \
        > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"
    )
    rc=$?
  else
    /usr/bin/time -v -o "${run_dir}/time_verbose.txt" \
      "${R_BIN}" "${PROJECT_DIR}/benchmarks/run_secact_cosmx_v3_headtohead.R" \
      > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"
    rc=$?
  fi
  stage "END cells=${cells} method=${method} rc=${rc}"
  record_status "${cells}" "${method}" "end" "${rc}" "${run_dir}"
  refresh_summary
  return 0
}

compare_pair() {
  cells="$1"
  ess_dir="${OUT_ROOT}/$(printf 'cells%06d_full_spatialess_perm%05d' "${cells}" "${NPERM}")"
  v3_dir="${OUT_ROOT}/$(printf 'cells%06d_full_spatialcellchat_v3_perm%05d' "${cells}" "${NPERM}")"
  if is_complete "${ess_dir}" && is_complete "${v3_dir}"; then
    export SECACT_V3_CELLS="${cells}" SECACT_V3_PANEL=full
    stage "START comparison cells=${cells}"
    "${R_BIN}" "${PROJECT_DIR}/benchmarks/compare_secact_cosmx_v3_headtohead.R" \
      > "${OUT_ROOT}/$(printf 'cells%06d_full_comparison.log' "${cells}")" 2>&1
    rc=$?
    stage "END comparison cells=${cells} rc=${rc}"
  else
    stage "SKIP comparison cells=${cells}; one method incomplete"
  fi
  refresh_summary
}

stage "LADDER QUEUE START sizes=${SIZES} nperm=${NPERM} v3_memory_gib=${V3_MEMORY_GIB}"

# Finish SpatialESS at every scale first so an expensive official run cannot
# block the million-cell-capable method from reaching the full dataset.
for cells in ${SIZES}; do
  run_method "${cells}" spatialess
done

# Then run the independent official implementation under the stated
# conventional-memory envelope and compare every successfully completed pair.
for cells in ${SIZES}; do
  run_method "${cells}" spatialcellchat_v3
  compare_pair "${cells}"
done

stage "LADDER QUEUE FINISHED"
refresh_summary

