import csv
import os
from pathlib import Path
import resource
import subprocess
import sys

root=Path(sys.argv[1]).resolve()
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
rscript='/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript'
subprocess.run([rscript,str(root/'scripts/prepare_official_lung.R'),str(root)],env=env,check=True)
with (root/'inputs/official_lung_manifest.tsv').open() as f: cases=list(csv.DictReader(f,delimiter='\t'))
with (root/'results/gse306130_pilot_recovery.tsv').open() as f: cases.extend(csv.DictReader(f,delimiter='\t'))
for case in cases:
    for variant in ['baseline','candidate']:
        d=root/'results'/case['case']/f'{variant}_1';d.mkdir(parents=True,exist_ok=True)
        if (d/'exit_code.txt').exists():continue
        with (d/'run.log').open('w') as log:
            p=subprocess.run(['/usr/bin/time','-v','-o',str(d/'time.txt'),rscript,str(root/'scripts/worker.R'),
                str(root),case['case'],variant,'1'],env=env,stdout=log,stderr=subprocess.STDOUT)
        (d/'exit_code.txt').write_text(str(p.returncode)+'\n')
        print(case['case'],variant,p.returncode,flush=True)
        if p.returncode:raise SystemExit(p.returncode)
