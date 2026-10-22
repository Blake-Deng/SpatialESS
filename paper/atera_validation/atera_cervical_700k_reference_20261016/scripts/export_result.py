#!/usr/bin/env python3
"""Export a compact, auditable table from the terminal reference run logs."""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'cervical_full_715512_SpatialCellChat_formal_b100'
state = json.loads((ROOT / 'logs' / f'{TAG}.exit.json').read_text())
log = (ROOT / 'logs' / f'{TAG}.log').read_text(errors='replace')
timing = (ROOT / 'logs' / f'{TAG}.time.txt').read_text(errors='replace')
requested = re.search(r'cannot allocate vector of size ([\d.]+)\s*G[bi]?', log, re.I)
failure = re.search(r'Error: ([^\n]+)', log)
peak = re.search(r'Maximum resident set size \(kbytes\): (\d+)', timing)
elapsed = re.search(r'Elapsed \(wall clock\) time .*?: ([\d:.]+)', timing)

def seconds(value):
    parts = [float(x) for x in value.split(':')]
    out = 0.0
    for part in parts:
        out = out * 60 + part
    return out

row = {
    'tag': TAG, 'tissue': 'cervical', 'fixture': 'full', 'N': state['input_cells'],
    'nominal_category': '700k', 'engine': 'SpatialCellChat', 'phase': 'formal',
    'genes': state['genes'], 'lr': state['lr_records'], 'nboot': state['nboot'],
    'seed_use': state['seed_use'], 'status': state['status'],
    'returncode': state['returncode'], 'whole_process_seconds': seconds(elapsed[1]) if elapsed else None,
    'peak_RSS_GB_decimal': int(peak[1]) * 1024 / 1e9 if peak else None,
    'address_space_limit_GiB': state['address_space_limit_GiB'],
    'requested_allocation_GiB': float(requested[1]) if requested else None,
    'failure_reason': failure[1] if failure else '',
    'reference_result_generated': (ROOT / 'results' / TAG / 'complete.json').exists(),
    'input': state['input'],
}
if elapsed:
    state['runner_elapsed_seconds'] = state.get('elapsed_seconds')
    state['elapsed_seconds'] = seconds(elapsed[1])
    state['whole_process_seconds'] = seconds(elapsed[1])
(ROOT / 'logs' / f'{TAG}.exit.json').write_text(json.dumps(state, indent=2) + '\\n')
with (ROOT / 'tables/cervical_700k_reference_attempt.tsv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=row.keys(), delimiter='\t')
    writer.writeheader()
    writer.writerow(row)
print(row)
