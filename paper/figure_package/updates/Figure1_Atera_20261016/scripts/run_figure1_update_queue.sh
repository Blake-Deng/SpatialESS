#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
python scripts/redraw.py --figure1-only
python scripts/compose_figure1.py
python - <<'PY'
from pathlib import Path
import pandas as pd

root = Path.cwd()
for report in ['data_checks.tsv', 'export_checks.tsv', 'source_data_checks.tsv']:
    table = pd.read_csv(root / 'reports' / report, sep='\t')
    if not table.status.eq('PASS').all():
        raise SystemExit(f'{report} contains failed checks')
status = pd.read_csv(root / 'source_data/panels/Figure1I_completion.csv')
new = status.query("record_type == 'display_only_status'")
assert len(status) == 50
assert len(new) == 4 and new.status.eq('not_attempted').all()
assert set(new.N) == {200000, 400000}
assert set(new.tissue) == {'cervical'}
assert not status.query('N == 169940').displayed_in_panel.any()
assert set(status.query('displayed_in_panel').axis_category) == {
    '1k','5k','10k','25k','50k','100k','200k','400k','450k','700k'}
print('Figure 1 update queue passed all figure, source and export checks.')
PY
