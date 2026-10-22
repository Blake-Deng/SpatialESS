"""Serial, resumable, tmux-safe ablation runner with automatic table exports."""
import csv
import fcntl
import hashlib
import json
import os
import random
import resource
import signal
import subprocess
import sys
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT/'config/design.json').read_text())
ENV = dict(os.environ, **CONFIG['threads'])
ENV['PATH'] = str(Path(CONFIG['rscript']).parent) + ':/usr/bin:/bin'
CPU = min(os.sched_getaffinity(0))


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%S%z')


def limits():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    limit = int(CONFIG['address_space_limit_GiB']*1024**3)
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    os.sched_setaffinity(0, {CPU})


def log(*values):
    line = stamp() + ' ' + ' '.join(str(v) for v in values)
    print(line, flush=True)
    with (ROOT/'logs/queue.log').open('a') as handle:
        handle.write(line+'\n')


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def run_worker(case, variant, boot, rep, phase):
    tag = f'{case}_{variant}_b{boot}_rep{rep}'
    out = ROOT/'results'/phase/tag
    if (out/'summary.tsv').exists() and (ROOT/'logs'/f'{tag}.exit.json').exists():
        if json.loads((ROOT/'logs'/f'{tag}.exit.json').read_text())['returncode'] == 0:
            log('RESUME skip', phase, tag)
            return True
    log('START', phase, tag)
    begun = time.time()
    with (ROOT/'logs'/f'{phase}_{tag}.log').open('w') as handle:
        command = ['/usr/bin/time', '-v', '-o', str(ROOT/'logs'/f'{tag}.time.txt'),
                   CONFIG['rscript'], str(ROOT/'scripts/worker.R'), str(ROOT),
                   case, variant, str(boot), str(rep), phase]
        try:
            process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT,
                                       env=ENV, preexec_fn=limits, start_new_session=True)
            code = process.wait(timeout=CONFIG['timeout_seconds'])
            status = 'completed' if code == 0 else f'failed_exit_{code}'
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            code = 124; status = 'timeout'
    exit_data = dict(tag=tag, phase=phase, returncode=code, status=status,
                     elapsed_seconds=time.time()-begun, finished=stamp())
    (ROOT/'logs'/f'{tag}.exit.json').write_text(json.dumps(exit_data, indent=2)+'\n')
    log('END', phase, tag, status, f"{exit_data['elapsed_seconds']:.2f}s")
    return code == 0


def compare(case, variant, boot, rep, phase, reference='full'):
    command = [CONFIG['rscript'], str(ROOT/'scripts/compare.R'), str(ROOT), phase,
               case, str(boot), variant, str(rep), reference]
    path = ROOT/'logs'/f'compare_{phase}_{case}_{variant}_b{boot}_rep{rep}.log'
    with path.open('w') as handle:
        result = subprocess.run(command, env=ENV, preexec_fn=limits,
                                stdout=handle, stderr=subprocess.STDOUT)
    return result.returncode == 0


def summary():
    with (ROOT/'logs/summarize.log').open('a') as handle:
        subprocess.run([sys.executable, str(ROOT/'scripts/summarize.py')],
                       stdout=handle, stderr=subprocess.STDOUT, check=True)


def update_validation(plan):
    for row in plan:
        out = ROOT/'results/formal'/row['tag']
        if not (out/'summary.tsv').exists() or (out/'validation.tsv').exists():
            continue
        ref = ROOT/'results/formal'/f"{row['case']}_full_b{row['nboot']}_rep{row['repetition']}"
        if not (ref/'summary.tsv').exists():
            continue
        ok = compare(row['case'], row['variant'], row['nboot'], row['repetition'], 'formal')
        log('VALIDATION', row['tag'], 'PASS' if ok else 'FAIL')


def build_plan():
    manifest = {row['case']: row for row in read_tsv(ROOT/'reports/input_manifest.tsv')}
    rng = random.Random(CONFIG['seed_for_run_order'])
    plan = []
    for block in ['core', 'permutation_diagnostic', 'extension']:
        for rep in range(1, CONFIG['repetitions']+1):
            batch = []
            if block in ['core', 'extension']:
                cases = CONFIG['core_cases'] if block == 'core' else CONFIG['extension_cases']
                settings = [(case, int(manifest[case]['nboot']), v)
                            for case in cases for v in CONFIG['variants']]
            else:
                d = CONFIG['permutation_diagnostic']
                settings = [(d['case'], boot, v) for boot in d['nboot'] for v in d['variants']]
            for case, boot, variant in settings:
                batch.append(dict(block=block, case=case, cells=int(manifest[case]['cells']),
                                  variant=variant, nboot=boot, repetition=rep,
                                  tag=f'{case}_{variant}_b{boot}_rep{rep}'))
            rng.shuffle(batch); plan.extend(batch)
    with (ROOT/'config/run_plan.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=plan[0].keys(), delimiter='\t')
        writer.writeheader(); writer.writerows(plan)
    return plan


def archive():
    output = ROOT/'SpatialESS_ablation_results.zip'
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED) as bundle:
        for folder in ['tables', 'reports', 'config', 'scripts', 'logs']:
            for path in sorted((ROOT/folder).rglob('*')):
                if path.is_file() and '__pycache__' not in path.parts:
                    bundle.write(path, path.relative_to(ROOT))
        for folder in ['source/vanilla', 'source/instrumented']:
            for path in sorted((ROOT/folder).rglob('*')):
                if path.is_file() and path.suffix not in {'.o', '.so'}:
                    bundle.write(path, path.relative_to(ROOT))
        for path in sorted((ROOT/'results').rglob('*')):
            if path.is_file() and path.name in {'summary.tsv', 'validation.tsv', 'contract.rds', 'sessionInfo.txt'}:
                bundle.write(path, path.relative_to(ROOT))
        bundle.write(ROOT/'README.md', 'README.md')


def main():
    lock = (ROOT/'logs/queue.lock').open('w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit('A queue is already running in this directory.')
    os.sched_setaffinity(0, {CPU})
    info = dict(started=stamp(), queue_pid=os.getpid(), cpu=CPU, thread_controls=CONFIG['threads'],
                rscript=CONFIG['rscript'], platform=sys.platform,
                address_space_limit_GiB=CONFIG['address_space_limit_GiB'])
    (ROOT/'reports/runner_environment.json').write_text(json.dumps(info, indent=2)+'\n')
    log('QUEUE START', 'CPU', CPU, 'single-thread', 'root', ROOT)
    for folder in ['source/vanilla', 'source/instrumented']:
        entries = []
        for path in sorted((ROOT/folder).rglob('*')):
            if path.is_file() and path.suffix not in {'.o', '.so'}:
                entries.append(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+str(path.relative_to(ROOT)))
        (ROOT/'reports'/('source_'+Path(folder).name+'_SHA256SUMS.txt')).write_text('\n'.join(entries)+'\n')
    manifest = {row['case']: row for row in read_tsv(ROOT/'reports/input_manifest.tsv')}
    plan = build_plan(); summary()
    if not (ROOT/'reports/preflight_pass.json').exists():
        # First verify a tiny input including the nonzero group-fraction filter.
        for case, phase, boot in [('cosmx_5000', 'tiny', 5),
                                  ('cosmx_5000', 'preflight', 20),
                                  ('cervical_center_5000_default', 'preflight', 100),
                                  ('cervical_stratified_5000_default', 'preflight', 100)]:
            for variant in ['vanilla', *CONFIG['variants']]:
                if not run_worker(case, variant, boot, 0, phase):
                    raise RuntimeError(f'Preflight worker failed: {case} {variant}; inspect logs')
                reference = 'vanilla' if variant == 'full' else 'full'
                if variant != 'vanilla' and not compare(case, variant, boot, 0, phase, reference):
                    raise RuntimeError(f'Preflight numerical validation failed: {case} {variant}')
        checks = []
        for phase in ['tiny', 'preflight']:
            for path in sorted((ROOT/'results'/phase).glob('*/validation.tsv')):
                checks.extend(read_tsv(path))
        with (ROOT/'reports/preflight_validation.tsv').open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=checks[0].keys(), delimiter='\t')
            writer.writeheader(); writer.writerows(checks)
        (ROOT/'reports/preflight_pass.json').write_text(json.dumps(dict(
            status='PASS', finished=stamp(), checks='Vanilla vs Full; all variants vs Full; exact graph edges/distances',
            cases=['tiny filtered CosMx', 'CosMx 5k', 'Xenium central 5k', 'Xenium stratified 5k']), indent=2)+'\n')
    log('PREFLIGHT PASS; formal jobs', len(plan))
    for i, row in enumerate(plan, 1):
        (ROOT/'reports/progress.json').write_text(json.dumps(dict(
            status='RUNNING', position=i, planned=len(plan), current=row, updated=stamp()), indent=2)+'\n')
        run_worker(row['case'], row['variant'], row['nboot'], row['repetition'], 'formal')
        update_validation(plan); summary()
    update_validation(plan); summary()
    status = json.loads((ROOT/'reports/table_status.json').read_text())
    status['status'] = 'COMPLETE' if status['validated'] == len(plan) else 'FINISHED_WITH_ISSUES'
    status['finished'] = stamp()
    (ROOT/'reports/progress.json').write_text(json.dumps(status, indent=2)+'\n')
    archive()
    log('QUEUE FINISHED', json.dumps(status), 'tables/SpatialESS_ablation_table.xlsx')


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        (ROOT/'reports/queue_error.json').write_text(json.dumps(dict(
            status='FAILED', error=repr(exc), time=stamp()), indent=2)+'\n')
        log('QUEUE FAILED', repr(exc))
        raise
