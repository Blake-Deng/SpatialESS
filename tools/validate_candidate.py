from pathlib import Path
import tempfile,tarfile,subprocess,json,csv,hashlib,shutil,sys,os
import argparse
p=argparse.ArgumentParser();p.add_argument('--delivery',type=Path,required=True);args=p.parse_args()
out=args.delivery.resolve();(out/'reports').mkdir(exist_ok=True);tmp=Path(tempfile.mkdtemp(prefix='SpatialESS-candidate-'));(out/'reports/VALIDATION_TEMP_PATH.txt').write_text(str(tmp)+'\n');results=[]
def call(name,cmd,cwd,env=None,expect=0):
 with (out/'reports'/f'{name}.log').open('w') as f:r=subprocess.run(cmd,cwd=cwd,stdout=f,stderr=subprocess.STDOUT,env=env)
 result=dict(check=name,exit_code=r.returncode,expected=expect,passed=r.returncode==expect);results.append(result);print(result,flush=True);(out/'reports/ACTUAL_RUNS.json').write_text(json.dumps(results,indent=2)+'\n');return r.returncode
# Extract actual candidate archives, not the development directory.
for archive in (out/'packages').glob('*.tar.gz'):
 with tarfile.open(archive) as tf:tf.extractall(tmp,filter='data')
sw=tmp/'SpatialESS_software_20261023';mat=tmp/'SpatialESS_reproduction_20261023'
call('software_manifest',[sys.executable,'tools/verify_software.py'],sw)
call('companion_manifest',[sys.executable,'tools/reproduction.py','verify','--materials',str(mat)],sw)
contracts=json.loads((sw/'docs/archive_20261023/FIGURE_INPUTS.json').read_text());rows={r['logical_path']:r for r in csv.DictReader((mat/'MATERIALS.tsv').open(),delimiter='\t')}
for target,paths in contracts.items():
 trial=tmp/target;trial.mkdir()
 with tarfile.open(out/'packages/SpatialESS_software_20261023.tar.gz') as tf:tf.extractall(trial,filter='data')
 source=trial/'SpatialESS_software_20261023';subset=trial/'materials';(subset/'objects').mkdir(parents=True)
 for p in paths:
  r=rows[p];dst=subset/r['object_path']
  if not dst.exists():shutil.copy2(mat/r['object_path'],dst)
 call('isolated_'+target,[sys.executable,'tools/redraw.py',target,'--materials',str(subset)],source)
 # No other figure-specific material exists here. Input isolation is tested by absence.
 if target=='Figure1':
  asset=source/'paper/figure_package/updates/Figure1_Atera_20261016/assets/Figure1A_user_original.pdf'
  data=asset.read_bytes();asset.unlink();asset.write_bytes(data+b'\nTAMPER_TEST\n')
  code=call('fixed_hash_rejects_tamper',[sys.executable,'paper/figure_package/updates/Figure1_Atera_20261016/scripts/compose_figure1.py'],source,expect=1)
  assert 'Frozen input SHA-256 mismatch' in (out/'reports/fixed_hash_rejects_tamper.log').read_text()
# Package build and real regression suite from the clean extracted software.
lib=tmp/'R-library';lib.mkdir()
call('R_build',['R','CMD','build','--no-build-vignettes','--no-manual',str(sw)],tmp)
call('R_install',['R','CMD','INSTALL','--library='+str(lib),str(tmp/'SpatialESS_0.1.4.tar.gz')],tmp)
env=dict(os.environ,R_LIBS_USER=str(lib))
call('R_regression_tests',['Rscript','-e','library(SpatialESS); testthat::test_dir("tests/testthat", reporter="summary", stop_on_failure=TRUE)'],sw,env)
call('ablation_source_install',['Rscript','tools/install_ablation.R',str(sw)],sw,env)
(out/'reports/ACTUAL_RUNS.json').write_text(json.dumps(results,indent=2)+'\n')
print('Candidate extraction:',tmp,flush=True)
