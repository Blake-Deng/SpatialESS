import csv, hashlib
from pathlib import Path
root = Path(__file__).resolve().parents[1]
with (root / 'SUBMISSION_SHA256.tsv').open() as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
names = [row['file'] for row in rows]
if len(names) != len(set(names)):
    raise SystemExit('FAIL duplicate manifest paths')
for row in rows:
    path = root / row['file']
    if not path.resolve().is_relative_to(root):
        raise SystemExit('FAIL path outside package: ' + row['file'])
    if not path.is_file() or path.stat().st_size != int(row['bytes']):
        raise SystemExit('FAIL missing file or byte count: ' + row['file'])
    if hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
        raise SystemExit('FAIL ' + row['file'])
print('PASS:', len(rows), 'packaged SHA256 checksums')
