"""Resumable serial Atera queue; finite resources, pinned package, audited failures."""
import csv
import fcntl
import hashlib
import json
import os
import resource
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT/'config/design.json').read_text())
ENV = dict(os.environ, **CFG['threads'], PYTHONUNBUFFERED='1')
ENV['PATH'] = str(Path(CFG['rscript']).parent)+':/usr/bin:/bin'
CPU = min(os.sched_getaffinity(0))


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%S%z')


def write_json(path, value):
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2,default=str)+'\n')
    temp.replace(path)


def log(*values):
    line = stamp()+' '+' '.join(str(v) for v in values)
    print(line,flush=True)
    with (ROOT/'logs/queue.log').open('a') as handle:
        handle.write(line+'\n')


def limits():
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    limit = int(CFG['address_space_limit_GiB']*1024**3)
    resource.setrlimit(resource.RLIMIT_AS,(limit,limit))
    os.sched_setaffinity(0,{CPU})


def run(tag, command, timeout, meta=None):
    meta = meta or {}
    marker = ROOT/'logs'/f'{tag}.exit.json'
    if marker.exists():
        previous = json.loads(marker.read_text())
        # Failures are retained and not repeatedly retried without intervention.
        log('RESUME',tag,previous['status'])
        return previous['returncode']==0
    attempt = dict(tag=tag,status='RUNNING',started=stamp(),runner_pid=os.getpid(),
                   command=command,timeout_seconds=timeout,**meta)
    write_json(ROOT/'reports/progress.json',attempt)
    log('START',tag)
    begun = time.time()
    log_path = ROOT/'logs'/f'{tag}.log'
    time_path = ROOT/'logs'/f'{tag}.time.txt'
    with log_path.open('w') as handle:
        proc = subprocess.Popen(['/usr/bin/time','-v','-o',str(time_path),*command],
          stdout=handle,stderr=subprocess.STDOUT,env=ENV,preexec_fn=limits,start_new_session=True)
        attempt.update(child_pid=proc.pid,log=str(log_path))
        write_json(ROOT/'reports/progress.json',attempt)
        try:
            code = proc.wait(timeout=timeout)
            status = 'completed' if code==0 else f'failed_exit_{code}'
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGKILL); proc.wait()
            code,status=124,'timeout'
    attempt.update(returncode=code,status=status,finished=stamp(),
                   elapsed_seconds=time.time()-begun)
    write_json(marker,attempt)
    log('END',tag,status,f"{attempt['elapsed_seconds']:.2f}s")
    return code==0


def install():
    target = ROOT/'source/SpatialESS'
    marker = ROOT/'reports/source_snapshot.json'
    if not marker.exists():
        source = Path(CFG['source_release'])
        for name in ['R','src','man','tests']:
            shutil.copytree(source/name,target/name,dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns('*.o','*.so','__pycache__'))
        for name in ['DESCRIPTION','NAMESPACE','.Rbuildignore','README.md']:
            shutil.copy2(source/name,target/name)
        hashes = {}
        for path in sorted(target.rglob('*')):
            if path.is_file():
                hashes[str(path.relative_to(target))] = hashlib.sha256(path.read_bytes()).hexdigest()
        write_json(marker,dict(source=str(source),version=CFG['expected_version'],
                               snapshot_at=stamp(),sha256=hashes))
    if not (ROOT/'lib/SpatialESS/DESCRIPTION').exists():
        r = str(Path(CFG['rscript']).with_name('R'))
        if not run('install_spatialess',[r,'CMD','INSTALL','--preclean','--no-multiarch',
                    '-l',str(ROOT/'lib'),str(target)],1800):
            raise RuntimeError('Isolated SpatialESS install failed; see install log')
    check = f'.libPaths(c("{ROOT / "lib"}",.libPaths())); stopifnot(as.character(packageVersion("SpatialESS"))=="{CFG["expected_version"]}"); print(sessionInfo()); print(packageDescription("SpatialCellChat")); print(body(SpatialCellChat::computeCommunProb)); print(body(SpatialCellChat::computeAvgCommunProb))'
    if not run('package_environment',[CFG['rscript'],'-e',check],120):
        raise RuntimeError('Package version preflight failed')
    audit = ROOT/'logs/package_environment.log'
    write_json(ROOT/'reports/environment.json',dict(python=sys.version,runner_pid=os.getpid(),
      cpu=CPU,threads=CFG['threads'],rscript=CFG['rscript'],spatialess_version=CFG['expected_version'],
      reference_implementation_log_sha256=hashlib.sha256(audit.read_bytes()).hexdigest(),
      address_space_limit_GiB=CFG['address_space_limit_GiB']))


def task(tissue,fixture,engine,phase,boot):
    tag=f'{tissue}_{fixture}_{engine}_{phase}_b{boot}'
    timeout = CFG['reference_timeout_seconds'] if engine=='SpatialCellChat' else CFG['ess_timeout_seconds']
    command = [CFG['rscript'],str(ROOT/'scripts/worker.R'),str(ROOT),tissue,fixture,engine,phase,str(boot),tag]
    ok=run(tag,command,timeout,dict(tissue=tissue,fixture=fixture,engine=engine,phase=phase,nboot=boot))
    summarize()
    return tag,ok


def compare(ess,reference):
    return run('compare_'+ess,[CFG['python'],str(ROOT/'scripts/compare.py'),str(ROOT),ess,reference],180)


def summarize():
    with (ROOT/'logs/summarize.log').open('a') as handle:
        code=subprocess.run([CFG['python'],str(ROOT/'scripts/summarize.py')],
              env=ENV,stdout=handle,stderr=subprocess.STDOUT).returncode
    if code:
        log('WARNING summary export failed; raw runs preserved; inspect summarize.log')


def save_plan():
    rows=[]
    def add(tissue,fixture,engine,phase,boot):
        rows.append(dict(position=len(rows)+1,tissue=tissue,fixture=fixture,engine=engine,
                         phase=phase,nboot=boot,tag=f'{tissue}_{fixture}_{engine}_{phase}_b{boot}'))
    for tissue in CFG['dataset_order']:
        for phase,boot in [('smoke',5),('pilot',20)]:
            for engine in ['SpatialESS','SpatialCellChat']:
                add(tissue,'central_1000',engine,phase,boot)
    for tissue in CFG['dataset_order']:
        for fixture in [f'central_{n}' for n in CFG['central_sizes']]+['full']:
            add(tissue,fixture,'SpatialESS','formal',CFG['nboot'])
            if fixture=='central_1000':
                add(tissue,fixture,'SpatialCellChat','formal',CFG['nboot'])
    for tissue in CFG['dataset_order']:
        for size in CFG['sampling_sizes']:
            for rep in range(1,CFG['sampling_repetitions']+1):
                for design in ['stratified','block']:
                    add(tissue,f'{design}_{size}_rep{rep}','SpatialESS','formal',CFG['nboot'])
        add(tissue,'central_5000','SpatialCellChat','formal',CFG['nboot'])
    with (ROOT/'config/run_plan.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=rows[0].keys(),delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def main():
    lock=(ROOT/'logs/queue.lock').open('w')
    try:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit('Another queue already holds the lock')
    log('QUEUE START',ROOT,'single CPU',CPU)
    install()
    available=[]
    for tissue in CFG['dataset_order']:
        if not run('prepare_matrix_'+tissue,[CFG['python'],str(ROOT/'scripts/prepare_matrix.py'),tissue],7200):
            log('EXCLUDE tissue: matrix preparation failed',tissue); continue
        if not run('prepare_fixtures_'+tissue,[CFG['rscript'],str(ROOT/'scripts/prepare_fixtures.R'),str(ROOT),tissue],7200):
            log('EXCLUDE tissue: fixture preparation failed',tissue); continue
        ess,ok=task(tissue,'central_1000','SpatialESS','smoke',5)
        if not ok:
            log('EXCLUDE tissue: ESS smoke failed',tissue); continue
        reference,refok=task(tissue,'central_1000','SpatialCellChat','smoke',5)
        if refok:
            strict=compare(ess,reference)
            row=next(csv.DictReader((ROOT/'results'/ess/'agreement.tsv').open(),delimiter='\t'))
            safe=(row['same_input_contract']=='True' and row['finite']=='True' and
                  row['active_jaccard']=='1.0' and float(row['max_absolute_probability_error'])<=1e-12)
            write_json(ROOT/'reports'/f'{tissue}_preflight.json',dict(
              tissue=tissue,strict_pass=strict,probability_key_gate=safe,
              note='Strict failures remain explicit; p-value exceptions do not imply equivalence.'))
            if not safe:
                log('HALT tissue: smoke probability/key mismatch',tissue); continue
        else:
            log('WARNING reference smoke failed; ESS-only measurements will remain unpaired',tissue)
        ess,ok=task(tissue,'central_1000','SpatialESS','pilot',20)
        if not ok:
            log('EXCLUDE tissue: full-LR pilot failed',tissue); continue
        reference,refok=task(tissue,'central_1000','SpatialCellChat','pilot',20)
        if refok: compare(ess,reference)
        available.append(tissue)
    # Both tissues reach preflight before the large, potentially long jobs.
    for tissue in available:
        for fixture in [f'central_{n}' for n in CFG['central_sizes']]+['full']:
            ess,ok=task(tissue,fixture,'SpatialESS','formal',CFG['nboot'])
            if fixture=='central_1000':
                ref,refok=task(tissue,fixture,'SpatialCellChat','formal',CFG['nboot'])
                if ok and refok: compare(ess,ref)
    for tissue in available:
        for size in CFG['sampling_sizes']:
            for rep in range(1,CFG['sampling_repetitions']+1):
                for design in ['stratified','block']:
                    task(tissue,f'{design}_{size}_rep{rep}','SpatialESS','formal',CFG['nboot'])
        ess=f'{tissue}_central_5000_SpatialESS_formal_b{CFG["nboot"]}'
        ref,refok=task(tissue,'central_5000','SpatialCellChat','formal',CFG['nboot'])
        if refok and (ROOT/'results'/ess/'complete.json').exists(): compare(ess,ref)
    summarize()
    failures=[json.loads(p.read_text()) for p in (ROOT/'logs').glob('*.exit.json')
              if json.loads(p.read_text())['returncode']!=0]
    status='FINISHED_WITH_ISSUES' if failures or len(available)!=len(CFG['dataset_order']) else 'COMPLETE'
    write_json(ROOT/'reports/progress.json',dict(status=status,finished=stamp(),
      included_tissues=available,failed_tasks=[r['tag'] for r in failures]))
    log('QUEUE FINISHED',status)


if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        write_json(ROOT/'reports/queue_error.json',dict(status='FAILED',error=repr(exc),time=stamp()))
        write_json(ROOT/'reports/progress.json',dict(status='FAILED',error=repr(exc),time=stamp()))
        log('QUEUE FAILED',repr(exc))
        raise
