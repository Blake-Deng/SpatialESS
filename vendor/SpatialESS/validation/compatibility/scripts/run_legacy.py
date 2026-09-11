"""Replay representative historical inputs with frozen baseline and candidate."""
import os
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parent
script = root / 'work/SpatialESS/benchmarks/run_lr_interface_regression.R'
oldlib = '/home/dzf/cellchat_acceleration/releases/SpatialESS_LR_interface_20260910/lib_archive'
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
for n in (1000, 5000, 10000, 50000):
    for name, lib, mode in [('baseline', oldlib, 'old_manual'),
                            ('candidate', str(root / 'lib_new'), 'new_explicit')]:
        out = root / 'validation/legacy' / f'cosmx_{n}' / name
        out.mkdir(parents=True, exist_ok=True)
        with (out / 'run.log').open('w') as stream:
            cmd = ['/usr/bin/time', '-v', '-o', str(out / 'time.txt'),
                   '/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript', str(script),
                   lib, 'cosmx', str(n), mode, str(out)]
            run = subprocess.run(cmd, env=env, stdout=stream, stderr=subprocess.STDOUT)
        (out / 'exit_code.txt').write_text(str(run.returncode)+'\n')
        print(n, name, run.returncode, flush=True)
        if run.returncode:
            raise SystemExit(run.returncode)
