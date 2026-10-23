"""Verify and mount the versioned companion data without duplicating shared files."""
import argparse,csv,hashlib,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''):h.update(block)
 return h.hexdigest()
def records():
 with (ROOT/'docs/archive_20261023/MATERIALS.tsv').open() as f:return list(csv.DictReader(f,delimiter='\t'))
def mount(materials,selected=None):
 rows=records();chosen=rows if selected is None else [r for r in rows if r['logical_path'] in selected]
 byhash={r['canonical_sha256']:r for r in chosen}
 for h,r in byhash.items():
  p=materials/r['object_path']
  if not p.is_file() or sha(p)!=h or p.stat().st_size!=int(r['canonical_bytes']):raise RuntimeError('Missing/corrupt material: '+str(p))
 for r in chosen:
  dest=ROOT/r['logical_path'];src=(materials/r['object_path']).resolve()
  if dest.exists() or dest.is_symlink():
   if dest.resolve()!=src:raise RuntimeError('Refusing to replace existing file: '+str(dest))
  else:
   dest.parent.mkdir(parents=True,exist_ok=True);dest.symlink_to(os.path.relpath(src,dest.parent))
 print(f'PASS {len(byhash)} unique objects, {len(chosen)} logical paths',flush=True)
 return chosen
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('command',choices=['verify','mount']);p.add_argument('--materials',type=Path,default=Path(os.environ.get('SPATIALESS_MATERIALS',ROOT.parent/'SpatialESS_reproduction_20261023')));a=p.parse_args()
 if a.command=='mount':mount(a.materials)
 else:
  seen=set()
  for r in records():
   if r['canonical_sha256'] in seen:continue
   path=a.materials/r['object_path']
   if not path.is_file() or sha(path)!=r['canonical_sha256'] or path.stat().st_size!=int(r['canonical_bytes']):raise RuntimeError('Missing/corrupt '+str(path))
   seen.add(r['canonical_sha256'])
  print('PASS',len(seen),'unique companion objects')
