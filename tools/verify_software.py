import csv,hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((root/'SUBMISSION_SHA256.tsv').open(),delimiter='\t'))
for r in rows:
 p=root/r['file']
 if not p.is_file() or p.stat().st_size!=int(r['bytes']) or hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise SystemExit('FAIL '+r['file'])
print('PASS',len(rows),'software file hashes (including source aliases)')
