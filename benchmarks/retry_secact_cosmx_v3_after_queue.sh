#!/usr/bin/env bash
set -u

QUEUE_PID="${1:-}"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_ROOT="${SECACT_V3_OUT:-/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead}"
LOG="${OUT_ROOT}/parallel_ladder_retry.log"

printf '[%s] RETRY watcher start queue_pid=%s\n' "$(date -Is)" "${QUEUE_PID}" >> "${LOG}"
if [[ -n "${QUEUE_PID}" ]]; then
  while kill -0 "${QUEUE_PID}" 2>/dev/null; do
    sleep 30
  done
fi

# Audit every requested scale after the concurrent pass. Successful directories
# are skipped; interrupted jobs and runs blocked by the old future globals
# threshold are repeated under the 120 GiB OS limit.
SECACT_V3_ESS_SIZES="200000 443515" \
SECACT_V3_SIZES="5000 10000 25000 50000 100000 200000 443515" \
SECACT_V3_MAX_V3_JOBS=1 \
  exec bash "${PROJECT_DIR}/benchmarks/run_secact_cosmx_v3_parallel_ladder.sh" \
  >> "${LOG}" 2>&1
