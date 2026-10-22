#!/usr/bin/env bash
set -euo pipefail

FIGURE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ANALYSIS_ROOT="/home/dzf/cellchat_acceleration/analysis/atera_cervical_200k_400k_20261015"
TAG="cervical_central_400000_SpatialESS_formal_b100"
EXIT_JSON="$ANALYSIS_ROOT/logs/$TAG.exit.json"
PERFORMANCE="$ANALYSIS_ROOT/tables/performance.tsv"
EXTENSION_SOURCE="$FIGURE_ROOT/source_data/current/cervical_extension_200k_400k.tsv"
LOG="$FIGURE_ROOT/logs/finalize_after_cervical_400k.log"
PYTHON="/home/dzf/miniforge3/bin/python"

mkdir -p "$(dirname "$LOG")"
printf '%s Waiting for terminal 400k SpatialESS run state.\n' "$(date --iso-8601=seconds)" >> "$LOG"
while :; do
    if [[ -f "$EXIT_JSON" && -f "$PERFORMANCE" ]]; then
        ready="$($PYTHON - "$EXIT_JSON" "$PERFORMANCE" <<'PY'
import csv
import json
import sys
exit_data = json.load(open(sys.argv[1]))
rows = list(csv.DictReader(open(sys.argv[2]), delimiter='\t'))
matches = [r for r in rows if r.get('tag') == exit_data.get('tag')]
if exit_data.get('status') != 'RUNNING' and len(matches) == 1 and matches[0].get('status') not in ('RUNNING', 'not_attempted', ''):
    print('ready')
PY
        )"
        if [[ "$ready" == "ready" ]]; then
            break
        fi
    fi
    sleep 20
done

cp "$PERFORMANCE" "$EXTENSION_SOURCE"
printf '%s Synced validated scale table.\n' "$(date --iso-8601=seconds)" >> "$LOG"
cd "$FIGURE_ROOT"
"$PYTHON" scripts/redraw.py --figure1-only >> "$LOG" 2>&1
"$PYTHON" scripts/compose_figure1.py >> "$LOG" 2>&1
"$PYTHON" scripts/finalize_figure1_reports.py >> "$LOG" 2>&1

"$PYTHON" - "$FIGURE_ROOT" <<'PY' >> "$LOG" 2>&1
import sys
from pathlib import Path
import pandas as pd
import pymupdf as fitz
root = Path(sys.argv[1])
status = pd.read_csv(root/'source_data/panels/Figure1I_completion.csv')
rows = status.query("tissue == 'cervical' and method == 'SpatialESS' and N in [200000,400000]")
assert len(rows) == 2 and not rows.status.isin(['not_attempted','running','in_progress']).any()
panel_i = fitz.open(root/'panels/Figure1I.pdf')[0].get_text()
assert 'Input size' not in panel_i
assert all(x in panel_i for x in ['200k','400k'])
assert all((root/'figures'/name).is_file() for name in
           ['Figure1.pdf','Figure1.svg','Figure1.png','Figure1_data.pdf','Figure1_data.svg','Figure1_data.png'])
print('Final figure exports verified. Cervical SpatialESS statuses:', rows[['N','status']].to_dict('records'))
PY

printf '%s Figure 1 finalized from terminal run records.\n' "$(date --iso-8601=seconds)" >> "$LOG"
