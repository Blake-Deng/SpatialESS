"""Reproduce figures from packaged, immutable source data; no inference."""
from pathlib import Path
import hashlib
import os
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for directory in ['figures', 'panels', 'preview', 'reports', 'logs', 'source_data/panels']:
    (ROOT / directory).mkdir(parents=True, exist_ok=True)

manifest = pd.read_csv(ROOT / 'provenance/INPUT_MANIFEST.tsv', sep='\t')
for row in manifest.itertuples():
    path = ROOT / row.staged_file
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != row.sha256:
        raise RuntimeError(f'Source checksum mismatch: {row.staged_file}')
print(f'Verified {len(manifest)} original input hashes.', flush=True)

for script, log in [('redraw.py', 'redraw.log'),
                    ('compose_and_check.py', 'export_qa.log'),
                    ('write_reports.py', 'reports.log')]:
    print(f'Running {script}', flush=True)
    with (ROOT / 'logs' / log).open('w') as handle:
        subprocess.run([sys.executable, str(ROOT / 'scripts' / script)],
                       cwd=ROOT, env=os.environ.copy(), stdout=handle,
                       stderr=subprocess.STDOUT, check=True)
print('Finished. Open reports/QA_REPORT.md and preview/*_from_pdf.png.')
