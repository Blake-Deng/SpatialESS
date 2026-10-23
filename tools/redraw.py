"""Portable entry points for the current main and supplementary figures."""
import argparse,os,subprocess,sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('figure',choices=['Figure1','Figure2','S1','S2','S3','S4']);p.add_argument('--materials',type=Path);args=p.parse_args()
if args.materials:
 from reproduction import mount
 mount(args.materials, json.loads((ROOT/'docs/archive_20261023/FIGURE_INPUTS.json').read_text())[args.figure])
paths={'Figure1':('paper/figure_package/updates/Figure1_Atera_20261016/scripts/redraw.py','--figure1-only'), 'Figure2':('paper/figure_package/updates/Figure2_20261012/scripts/redraw.py','--figure2-only')}
if args.figure in paths:
 script,flag=paths[args.figure];cmd=[sys.executable,str(ROOT/script),flag]
else:cmd=[sys.executable,str(ROOT/'paper/figure_package/scripts/redraw.py'),'--figure',args.figure]
# Historical environment logs are evidence; keep them outside the new run's outputs.
base=Path(cmd[1]).parents[1]
for d in ['figures','panels','preview','reports','logs','source_data/panels']:(base/d).mkdir(parents=True,exist_ok=True)
# Mounted outputs are immutable inputs. New exports replace the link, never the object.
for folder in ['figures','panels','preview','reports']:
 for dest in (base/folder).rglob('*'):
  if dest.is_symlink():dest.unlink()
for dest in list((base/'source_data/panels').glob('*'))+[base/'source_data/Figure1_Source_Data.xlsx',base/'source_data/Figure2_Source_Data.xlsx',base/'source_data/feature_mapping.csv',base/'source_data/abbreviations.csv']:
 if dest.is_symlink():dest.unlink()
# Reports and environment snapshots from this run get a separate location.
logbase=ROOT/'local_runs'/args.figure;logbase.mkdir(parents=True,exist_ok=True)
envfile=base/'logs/environment.json'
if envfile.is_symlink():envfile.unlink()
elif envfile.exists():envfile.replace(logbase/'previous_environment.json')
with (logbase/'redraw.log').open('w') as f:subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)
if args.figure=='Figure1':
 with (logbase/'composition.log').open('w') as f:subprocess.run([sys.executable,str(base/'scripts/compose_figure1.py')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)
if envfile.exists():envfile.replace(logbase/'environment.json')
print('Completed',args.figure,'from saved inputs. Run logs:',logbase)
