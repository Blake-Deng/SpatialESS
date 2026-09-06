#!/usr/bin/env bash
set -u

QUEUE_PID="${1:-}"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_ROOT="${SECACT_V3_OUT:-/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead}"
QUEUE_LOG="${OUT_ROOT}/parallel_ladder_queue.log"
WATCH_LOG="${OUT_ROOT}/parallel_ladder_watchdog.log"

printf '[%s] WATCH start queue_pid=%s\n' "$(date -Is)" "${QUEUE_PID}" >> "${WATCH_LOG}"

if [[ -n "${QUEUE_PID}" ]]; then
  while kill -0 "${QUEUE_PID}" 2>/dev/null; do
    sleep 30
  done
fi

if tail -n 1 "${QUEUE_LOG}" 2>/dev/null | grep -q 'PARALLEL QUEUE FINISHED'; then
  printf '[%s] WATCH queue completed normally\n' "$(date -Is)" >> "${WATCH_LOG}"
  exit 0
fi

printf '[%s] WATCH queue disappeared before completion; resuming\n' "$(date -Is)" >> "${WATCH_LOG}"
exec bash "${PROJECT_DIR}/benchmarks/run_secact_cosmx_v3_parallel_ladder.sh" \
  >> "${OUT_ROOT}/parallel_ladder_resume.log" 2>&1

