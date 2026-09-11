"""Run independent cases concurrently, preserving exit status and GNU time logs."""
import os
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parent
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
rscript = '/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript'
jobs = []
for case in ('center_1000', 'center_5000', 'stratified_1000', 'stratified_5000'):
    log = (root / 'logs' / f'{case}.log').open('w')
    p = subprocess.Popen(['/usr/bin/time', '-v', '-o', str(root / 'logs' / f'{case}.time.txt'),
                          rscript, str(root / 'run_dual_case.R'), case],
                         env=env, stdout=log, stderr=subprocess.STDOUT)
    jobs.append((case, p, log))
    print('started', case, p.pid, flush=True)
for case, p, log in jobs:
    code = p.wait()
    log.close()
    (root / 'logs' / f'{case}.exit_code.txt').write_text(str(code)+'\n')
    print('finished', case, code, flush=True)
raise SystemExit(int(any(p.returncode != 0 for _, p, _ in jobs)))
