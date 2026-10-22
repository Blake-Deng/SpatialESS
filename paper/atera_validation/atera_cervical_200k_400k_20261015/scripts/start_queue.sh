#!/usr/bin/env bash
set -uo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$task_root"
/home/dzf/miniforge3/bin/python -u scripts/queue_scale_v2.py >> logs/tmux_console.log 2>&1
queue_code=$?
printf '\nQUEUE_EXIT=%s\n' "$queue_code" >> logs/tmux_console.log
exec bash
