"""Replay frozen RAW and corrected ovarian inputs through the standalone bridge."""
import json
import os
import shutil
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parent
old = Path('/home/dzf/cellchat_acceleration/releases/SpatialESS_main_SPARKLE_20260910')
out = root / 'validation/bridge'
out.mkdir(parents=True, exist_ok=True)
shutil.copy2(old / 'ovarian_historical_database.rds', out)
cases = json.loads((old / 'inputs/cases.json').read_text())
(out / 'cases.json').write_text(json.dumps(cases, indent=2)+'\n')
env = dict(os.environ, R_LIBS_USER=str(root / 'lib_new'),
           PYTHON='/home/dzf/miniforge3/bin/python',
           OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
for case in cases:
    for mode in ['new_bridge', 'direct_main']:
        d = out / 'runs' / f"{case['condition']}_{case['cells']}" / mode
        d.mkdir(parents=True, exist_ok=True)
        with (d / 'run.log').open('w') as log:
            p = subprocess.run(['/usr/bin/time', '-v', '-o', str(d / 'time.txt'),
                '/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript',
                str(root / 'work/SpatialESS-SPARKLE/benchmark/run_bridge_regression.R'),
                str(out), case['condition'], str(case['cells']), mode, case['h5ad'], case['metadata']],
                env=env, stdout=log, stderr=subprocess.STDOUT)
        (d / 'exit_code.txt').write_text(str(p.returncode)+'\n')
        print(case['condition'], case['cells'], mode, p.returncode, flush=True)
        if p.returncode:
            raise SystemExit(p.returncode)
