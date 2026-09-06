#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PATH="/home/dzf/miniforge3/envs/cellchat-acceleration/bin:${PATH}"
export SECACT_ROOT="${SECACT_ROOT:-/data/dzf/SecAct_2026}"
export SPATIALESS_NPERM="${SPATIALESS_NPERM:-100}"
export SPATIALESS_RADIUS_UM="${SPATIALESS_RADIUS_UM:-20}"
export SPATIALESS_BLOCK_UM="${SPATIALESS_BLOCK_UM:-100}"
export SPATIALESS_OUT="${SPATIALESS_OUT:-${SECACT_ROOT}/results/secact_cosmx_radius${SPATIALESS_RADIUS_UM}_perm${SPATIALESS_NPERM}}"

mkdir -p "${SPATIALESS_OUT}"
/usr/bin/time -v -o "${SPATIALESS_OUT}/time_verbose.txt" \
  Rscript "${PROJECT_DIR}/benchmarks/run_secact_cosmx_spatialess.R" \
  >"${SPATIALESS_OUT}/run.log" 2>&1

printf 'Completed. Results: %s\n' "${SPATIALESS_OUT}"

