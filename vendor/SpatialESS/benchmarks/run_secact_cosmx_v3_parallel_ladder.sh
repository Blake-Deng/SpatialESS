#!/usr/bin/env bash
set -u

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R_BIN="${SPATIAL_V3_R_BIN:-/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript}"
OUT_ROOT="${SECACT_V3_OUT:-/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead}"
SIZES="${SECACT_V3_SIZES:-5000 10000 25000 50000 100000 200000 443515}"
ESS_SIZES="${SECACT_V3_ESS_SIZES:-200000 443515}"
NPERM="${SECACT_V3_NPERM:-20}"
MAX_V3_JOBS="${SECACT_V3_MAX_V3_JOBS:-5}"
V3_MEMORY_GIB="${SECACT_V3_MEMORY_GIB:-120}"
V3_MEMORY_KB=$((V3_MEMORY_GIB * 1024 * 1024))

export PATH="$(dirname "${R_BIN}"):/usr/bin:/bin"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export SECACT_V3_OUT="${OUT_ROOT}" SECACT_V3_NPERM="${NPERM}"
mkdir -p "${OUT_ROOT}"

QUEUE_LOG="${OUT_ROOT}/parallel_ladder_queue.log"
STATUS_TSV="${OUT_ROOT}/parallel_ladder_status.tsv"
LOCK_DIR="${OUT_ROOT}/.parallel_summary.lock"

if [[ ! -f "${STATUS_TSV}" ]]; then
  printf 'timestamp\tcells\tmethod\tevent\texit_status\tmemory_limit_gib\trun_directory\n' > "${STATUS_TSV}"
fi

log_event() {
  printf '[%s] %s\n' "$(date -Is)" "$1" >> "${QUEUE_LOG}"
}

record_status() {
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date -Is)" "$1" "$2" "$3" "$4" "${V3_MEMORY_GIB}" "$5" >> "${STATUS_TSV}"
}

run_dir_for() {
  printf '%s/cells%06d_full_%s_perm%05d' "${OUT_ROOT}" "$1" "$2" "${NPERM}"
}

is_complete() {
  local run_dir="$1"
  [[ -s "${run_dir}/summary.tsv" ]] &&
    [[ -s "${run_dir}/records.rds" ]] &&
    [[ -s "${run_dir}/time_verbose.txt" ]] &&
    grep -q 'Exit status: 0' "${run_dir}/time_verbose.txt"
}

refresh_summary() {
  if mkdir "${LOCK_DIR}" 2>/dev/null; then
    SECACT_V3_SIZES="250 1000 5000 10000 25000 50000 100000 200000 443515" \
      "${R_BIN}" "${PROJECT_DIR}/benchmarks/summarize_secact_cosmx_v3_ladder.R" \
      > "${OUT_ROOT}/ladder_summary_latest.log" 2>&1 || true
    rmdir "${LOCK_DIR}" 2>/dev/null || true
  fi
}

compare_if_ready() {
  local cells="$1"
  local ess_dir v3_dir comparison_dir
  ess_dir="$(run_dir_for "${cells}" spatialess)"
  v3_dir="$(run_dir_for "${cells}" spatialcellchat_v3)"
  comparison_dir="${OUT_ROOT}/$(printf 'cells%06d_full_comparison_perm%05d' "${cells}" "${NPERM}")"
  if is_complete "${ess_dir}" && is_complete "${v3_dir}"; then
    log_event "START comparison cells=${cells}"
    SECACT_V3_CELLS="${cells}" SECACT_V3_PANEL=full \
      "${R_BIN}" "${PROJECT_DIR}/benchmarks/compare_secact_cosmx_v3_headtohead.R" \
      > "${comparison_dir}.log" 2>&1
    local rc=$?
    log_event "END comparison cells=${cells} rc=${rc}"
  fi
  refresh_summary
}

run_one() {
  local cells="$1"
  local method="$2"
  local run_dir rc
  run_dir="$(run_dir_for "${cells}" "${method}")"
  mkdir -p "${run_dir}"
  if is_complete "${run_dir}"; then
    log_event "SKIP completed cells=${cells} method=${method}"
    record_status "${cells}" "${method}" skip_completed 0 "${run_dir}"
    compare_if_ready "${cells}"
    return 0
  fi

  log_event "START cells=${cells} method=${method}"
  record_status "${cells}" "${method}" start -1 "${run_dir}"
  if [[ "${method}" == spatialcellchat_v3 ]]; then
    (
      ulimit -v "${V3_MEMORY_KB}"
      SECACT_V3_CELLS="${cells}" SECACT_V3_PANEL=full SECACT_V3_METHOD="${method}" \
        /usr/bin/time -v -o "${run_dir}/time_verbose.txt" \
        "${R_BIN}" "${PROJECT_DIR}/benchmarks/run_secact_cosmx_v3_headtohead.R" \
        > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"
    )
    rc=$?
  else
    SECACT_V3_CELLS="${cells}" SECACT_V3_PANEL=full SECACT_V3_METHOD="${method}" \
      /usr/bin/time -v -o "${run_dir}/time_verbose.txt" \
      "${R_BIN}" "${PROJECT_DIR}/benchmarks/run_secact_cosmx_v3_headtohead.R" \
      > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"
    rc=$?
  fi
  log_event "END cells=${cells} method=${method} rc=${rc}"
  record_status "${cells}" "${method}" end "${rc}" "${run_dir}"
  compare_if_ready "${cells}"
  return 0
}

export -f log_event record_status run_dir_for is_complete refresh_summary compare_if_ready run_one
export PROJECT_DIR R_BIN OUT_ROOT NPERM V3_MEMORY_GIB V3_MEMORY_KB QUEUE_LOG STATUS_TSV LOCK_DIR

log_event "PARALLEL QUEUE START ess_sizes=${ESS_SIZES} v3_sizes=${SIZES} max_v3_jobs=${MAX_V3_JOBS}"

# SpatialESS repairs are independent and light enough to run together.
for cells in ${ESS_SIZES}; do
  run_one "${cells}" spatialess &
done
wait

# Each official job remains one-threaded and independently capped at 120 GiB.
# Parallelism is across dataset sizes, preserving the per-method benchmark definition.
printf '%s\n' ${SIZES} | xargs -n 1 -P "${MAX_V3_JOBS}" bash -c \
  'run_one "$1" spatialcellchat_v3' _

# Pick up the separately running 250-cell official job and every completed pair.
for cells in 250 1000 5000 10000 25000 50000 100000 200000 443515; do
  compare_if_ready "${cells}"
done
refresh_summary
log_event "PARALLEL QUEUE FINISHED"

