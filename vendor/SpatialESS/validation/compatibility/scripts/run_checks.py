"""Run local package checks with explicit toolchains and persistent statuses."""
import json
import os
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parent
rbin = Path('/home/dzf/miniforge3/envs/spatialcellchat-v3/bin')
env = dict(os.environ, PATH=f'{rbin}:/usr/local/bin:/usr/bin:/bin',
           R_LIBS_USER=str(root / 'lib_new'), PYTHON='/home/dzf/miniforge3/bin/python',
           RSCRIPT=str(rbin / 'Rscript'), OMP_NUM_THREADS='1',
           OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', _R_CHECK_FORCE_SUGGESTS_='false')
statuses = []

def run(name, args, cwd=root):
    with (root / 'logs' / f'{name}.log').open('w') as log:
        p = subprocess.run([str(x) for x in args], cwd=cwd, env=env,
                           stdout=log, stderr=subprocess.STDOUT)
    statuses.append(dict(check=name, exit_code=p.returncode))
    (root / 'logs/local_check_status.json').write_text(json.dumps(statuses, indent=2)+'\n')
    print(name, p.returncode, flush=True)
    if p.returncode:
        raise SystemExit(p.returncode)

run('sync_branch', [env['PYTHON'], root / 'sync_branch.py'])
main = root / 'work/SpatialESS'
bridge = root / 'work/SpatialESS-SPARKLE'
run('install_bridge', [rbin / 'R', 'CMD', 'INSTALL', '-l', root / 'lib_new', bridge])
run('mini_default', [rbin / 'Rscript', main / 'benchmarks/verify_lr_interface_mini.R', main])
run('main_unit_final', [rbin / 'Rscript', '-e', f'testthat::test_local("{main}", reporter="summary")'])
run('bridge_python', [env['PYTHON'], '-m', 'unittest', 'discover', '-s', 'tests/python', '-v'], bridge)
run('bridge_unit', [rbin / 'Rscript', '-e', f'testthat::test_local("{bridge}", reporter="summary")'])
run('bridge_demo', ['sh', bridge / 'examples/minimal/run_demo.sh'])
run('lmm_bundle', [rbin / 'Rscript', main / 'benchmarks/verify_gse250346_result_bundle.R',
                   main / 'results/tables'])
for path, pkg in [(main, 'SpatialESS'), (bridge, 'SpatialESSSPARKLE')]:
    run(f'{pkg}_build', [rbin / 'R', 'CMD', 'build', path])
    run(f'{pkg}_check', [rbin / 'R', 'CMD', 'check', '--no-manual', '--no-build-vignettes',
                         root / f'{pkg}_0.1.3.tar.gz'])
    text = (root / f'{pkg}.Rcheck/00check.log').read_text()
    if 'Status: OK' not in text:
        raise SystemExit(f'{pkg} check did not finish with Status: OK')
