from pathlib import Path
import os,json,subprocess,datetime
root=Path(__file__).resolve().parent
script=root/'work/SpatialESS-SPARKLE/benchmark/run_bridge_regression.R'
env=os.environ.copy()
env['PATH']='/home/dzf/miniforge3/envs/spatialcellchat-v3/bin:/usr/bin:/bin'
env['PYTHON']='/home/dzf/miniforge3/bin/python'
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    env[key]='1'
queue=root/'logs/queue.tsv'
with queue.open('w',buffering=1) as log:
    for c in json.loads((root/'inputs/cases.json').read_text()):
        for mode in ['old_bridge','new_bridge','direct_main']:
            out=root/'runs'/f"{c['condition']}_{c['cells']}"/mode; out.mkdir(parents=True)
            env['R_LIBS_USER']=(str(root/'lib_old_bridge')+':/home/dzf/cellchat_acceleration/releases/SpatialESS_LR_interface_20260910/lib_old') if mode=='old_bridge' else str(root/'lib')
            log.write(f"{datetime.datetime.now().isoformat()}\t{c['condition']}\t{c['cells']}\t{mode}\tstarted\n")
            cmd=['/usr/bin/time','-v','-o',str(out/'time.txt'),'Rscript',str(script),str(root),c['condition'],str(c['cells']),mode,c['h5ad'],c['metadata']]
            with (out/'stdout.log').open('w') as stdout,(out/'stderr.log').open('w') as stderr:
                run=subprocess.run(cmd,env=env,stdout=stdout,stderr=stderr)
            log.write(f"{datetime.datetime.now().isoformat()}\t{c['condition']}\t{c['cells']}\t{mode}\tfinished\t{run.returncode}\n")
            if run.returncode: raise SystemExit(run.returncode)
    log.write('ALL_COMPLETED\n')
