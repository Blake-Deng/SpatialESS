"""Inspect or download cited public Zenodo inputs and verify provider checksums.

Default only lists files. Explicit --file selections download potentially large
raw datasets. A resulting SHA256SUMS records the actual downloaded bytes.
"""
import argparse,hashlib,json,urllib.request
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--record',default='17664023');p.add_argument('--file',action='append',default=[]);p.add_argument('--out',type=Path,default=Path('public_inputs'));a=p.parse_args()
if not a.record.isdigit():raise SystemExit('Numeric Zenodo record required')
with urllib.request.urlopen('https://zenodo.org/api/records/'+a.record) as r:meta=json.load(r)
for f in meta['files']:print(f['key'],f['size'],f['checksum'])
if a.file:
 a.out.mkdir(parents=True,exist_ok=True);(a.out/'zenodo_source_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
 for name in a.file:
  f=next((x for x in meta['files'] if x['key']==name),None)
  if f is None:raise SystemExit('Not in source record: '+name)
  if Path(name).name!=name:raise SystemExit('Explicit flat filename required')
  dest=a.out/name
  if dest.exists():raise SystemExit('Refusing to overwrite '+str(dest))
  urllib.request.urlretrieve(f['links']['self'],dest)
  algorithm,expected=f['checksum'].split(':',1);h=hashlib.new(algorithm);s=hashlib.sha256()
  with dest.open('rb') as stream:
   for block in iter(lambda:stream.read(1<<20),b''):h.update(block);s.update(block)
  if dest.stat().st_size!=f['size'] or h.hexdigest()!=expected:raise SystemExit('Source checksum failed '+name)
  with (a.out/'SHA256SUMS.txt').open('a') as stream:stream.write(s.hexdigest()+'  '+name+'\n')
