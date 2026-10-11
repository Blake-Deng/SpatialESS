"""Package checked figure masters, source data and reproducible plotting code."""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path, default=ROOT.parent / (ROOT.name + '.zip'))
args = parser.parse_args()

required = [ROOT / 'figures' / f'{name}.{ext}'
            for name in ['Figure1_data', 'Figure2', 'FigureS1', 'FigureS2', 'FigureS3', 'FigureS4']
            for ext in ['pdf', 'svg', 'png']]
required += [ROOT / name for name in ['source_data/Source_Data.xlsx',
                                     'FIGURE_PROVENANCE.md', 'FIGURE_CAPTIONS.md',
                                     'reports/QA_REPORT.md', 'README.md']]
assert all(path.is_file() and path.stat().st_size for path in required)
assert not (pd.read_csv(ROOT / 'reports/data_checks.tsv', sep='\t').status == 'FAIL').any()
export_checks = pd.read_csv(ROOT / 'reports/export_checks.tsv', sep='\t')
assert not (export_checks.loc[export_checks.figure != 'Figure1A', 'status'] != 'PASS').any()

def include(path):
    return path.is_file() and '__pycache__' not in path.parts and path.suffix != '.pyc' \
        and path.name != 'PACKAGE_SHA256.tsv'

files = sorted(path for path in ROOT.rglob('*') if include(path))
rows = [{'file': str(path.relative_to(ROOT)),
         'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
         'bytes': path.stat().st_size} for path in files]
checksums = ROOT / 'PACKAGE_SHA256.tsv'
pd.DataFrame(rows).to_csv(checksums, sep='\t', index=False)
files.append(checksums)
args.output.parent.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(args.output, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
    for path in files:
        archive.write(path, str(Path(ROOT.name) / path.relative_to(ROOT)))
with zipfile.ZipFile(args.output) as archive:
    assert archive.testzip() is None
archive_hash = hashlib.sha256(args.output.read_bytes()).hexdigest()
args.output.with_suffix('.zip.sha256').write_text(f'{archive_hash}  {args.output.name}\n')
print(json.dumps({'zip': str(args.output), 'files': len(files),
                  'bytes': args.output.stat().st_size, 'sha256': archive_hash}, indent=2))
