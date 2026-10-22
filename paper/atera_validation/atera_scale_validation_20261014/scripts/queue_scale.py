"""Bounded sequential ladder with immutable reuse and recoverable running jobs."""
import fcntl
import json
import os
import platform
import resource
import signal
import subprocess
import sys
import time
from pathlib import Path

from scale_common import (BASE, CFG, ROOT, SIZES, archive_attempt, failure_class, fixture,
                          process_alive, sha256, stamp, tag, time_metrics,
                          validate_success, write_json)

ENV = dict(os.environ, **CFG['threads'], PYTHONUNBUFFERED='1')
ENV['PATH'] = str(Path(CFG['rscript']).parent) + ':/usr/bin:/bin'
CPU = int(CFG['cpu'])


def log(*values):
    line = stamp() + ' ' + ' '.join(str(x) for x in values)
    print(line, flush=True)
    with (ROOT / 'logs/queue.log').open('a') as handle:
        handle.write(line + '\n')


def child_limits():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    limit = int(CFG['address_space_limit_GiB'] * 1024 ** 3)
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    os.sched_setaffinity(0, {CPU})


def summarize():
    with (ROOT / 'logs/summarize.log').open('a') as handle:
        result = subprocess.run([CFG['python'], str(ROOT / 'scripts/summarize_scale.py')],
                                stdout=handle, stderr=subprocess.STDOUT, env=ENV)
    if result.returncode:
        raise RuntimeError('Summary export failed; see logs/summarize.log')


def sample_memory(attempt):
    pid = attempt.get('r_pid')
    if not pid:
        children = Path(f'/proc/{attempt["child_pid"]}/task/{attempt["child_pid"]}/children')
        if children.exists() and children.read_text().split():
            attempt['r_pid'] = int(children.read_text().split()[0])
            pid = attempt['r_pid']
    if pid:
        status = Path(f'/proc/{pid}/status')
        if status.exists():
            for line in status.read_text().splitlines():
                if line.startswith('VmHWM:'):
                    attempt['observed_R_high_water_RSS_KiB'] = max(
                        attempt.get('observed_R_high_water_RSS_KiB', 0), int(line.split()[1]))
    attempt['heartbeat'] = stamp()
    write_json(ROOT / 'logs' / (attempt['tag'] + '.running.json'), attempt)
    write_json(ROOT / 'reports/progress.json', attempt)


def finish(attempt, code, timed_out=False):
    name = attempt['tag']
    logfile = ROOT / 'logs' / (name + '.log')
    measures = time_metrics(ROOT / 'logs' / (name + '.time.txt'))
    if code == 0 and attempt.get('engine') and not validate_success(name):
        code = 3
    state = failure_class(code, logfile, timed_out)
    attempt.update(returncode=code, status=state, finished=stamp(),
                   elapsed_seconds=time.time() - attempt['started_epoch'], measurements=measures)
    write_json(ROOT / 'logs' / (name + '.exit.json'), attempt)
    running = ROOT / 'logs' / (name + '.running.json')
    if running.exists():
        running.unlink()
    log('END', name, state, 'seconds', round(attempt['elapsed_seconds'], 2))
    summarize()
    return code == 0


def existing_run(name):
    running = ROOT / 'logs' / (name + '.running.json')
    if not running.exists():
        return None
    attempt = json.loads(running.read_text())
    pid = attempt['child_pid']
    if process_alive(pid):
        cmdline = Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\0', b' ').decode(errors='replace')
        if name not in cmdline or str(ROOT) not in cmdline:
            raise RuntimeError('Recovered PID identity mismatch: ' + name)
    attempt.update(runner_pid=os.getpid(), adopted=True)
    log('ADOPT', name, 'existing time process', pid)
    remaining = max(0, attempt['timeout_seconds'] - (time.time() - attempt['started_epoch']))
    deadline = time.time() + remaining
    timed_out = False
    while process_alive(pid):
        sample_memory(attempt)
        if time.time() >= deadline:
            try:
                os.killpg(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            timed_out = True
            break
        time.sleep(min(15, max(.1, deadline - time.time())))
    measures = time_metrics(ROOT / 'logs' / (name + '.time.txt'))
    code = 124 if timed_out else measures.get('reported_exit_status')
    if not timed_out and measures.get('termination_signal'):
        code = 128 + measures['termination_signal']
    return finish(attempt, code, timed_out)


def run(name, command, timeout, meta):
    marker = ROOT / 'logs' / (name + '.exit.json')
    if marker.exists():
        previous = json.loads(marker.read_text())
        if previous.get('returncode') == 0:
            if meta.get('engine') and not validate_success(name):
                raise RuntimeError('Successful marker without complete validated result: ' + name)
            log('RESUME COMPLETED', name)
            return True
        log('RETAIN FAILED', name, previous['status'])
        return False
    adopted = existing_run(name)
    if adopted is not None:
        return adopted
    log('START', name, 'timeout', timeout)
    attempt = dict(tag=name, status='RUNNING', started=stamp(), started_epoch=time.time(),
                   command=command, timeout_seconds=timeout, runner_pid=os.getpid(), **meta)
    logfile = ROOT / 'logs' / (name + '.log')
    timefile = ROOT / 'logs' / (name + '.time.txt')
    with logfile.open('w') as handle:
        proc = subprocess.Popen(['/usr/bin/time', '-v', '-o', str(timefile)] + command,
                                stdout=handle, stderr=subprocess.STDOUT, env=ENV,
                                preexec_fn=child_limits, start_new_session=True)
        attempt['child_pid'] = proc.pid
        deadline = time.time() + timeout
        timed_out = False
        while True:
            sample_memory(attempt)
            try:
                code = proc.wait(timeout=min(15, max(.1, deadline - time.time())))
                break
            except subprocess.TimeoutExpired:
                if time.time() >= deadline:
                    try:
                        os.killpg(proc.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    proc.wait()
                    code, timed_out = 124, True
                    break
    return finish(attempt, code, timed_out)


def setup():
    for folder in ['logs', 'reports', 'inputs', 'results', 'comparisons', 'tables']:
        (ROOT / folder).mkdir(exist_ok=True)
    if not (ROOT / 'lib').exists():
        (ROOT / 'lib').symlink_to(BASE / 'lib')
    for tissue in CFG['tissues']:
        target = ROOT / 'inputs' / tissue
        target.mkdir(exist_ok=True)
        for path in (BASE / 'inputs' / tissue).iterdir():
            if not (target / path.name).exists():
                (target / path.name).symlink_to(path)
    if not (ROOT / 'inputs/database.rds').exists():
        (ROOT / 'inputs/database.rds').symlink_to(BASE / 'inputs/database.rds')
    for tissue in CFG['tissues']:
        for size in SIZES:
            name = tag(tissue, size, 'SpatialESS')
            source = BASE / 'results' / name
            target = ROOT / 'results' / name
            if source.exists() and not target.exists():
                target.symlink_to(source)
    # Archive only known setup failures, without automatic retries of inference failures.
    for marker in sorted((ROOT / 'logs').glob('*10000*.exit.json')):
        row = json.loads(marker.read_text())
        logfile = marker.with_name(marker.name.replace('.exit.json', '.log'))
        if (row.get('returncode') == 1 and logfile.exists() and
                'packageVersion("SpatialESS")' in logfile.read_text()):
            archive_attempt(row['tag'], 'Startup-only failure before inference; isolated 0.1.4 library path corrected.')
    write_json(ROOT / 'reports/runner_environment.json', dict(
        python=sys.version, cpu=CPU, platform=platform.platform(), threads=CFG['threads'],
        rscript=CFG['rscript'], address_space_limit_GiB=CFG['address_space_limit_GiB'],
        reference_timeout_seconds=CFG['reference_timeout_seconds'],
        source_root=CFG['source_root'], base_analysis=str(BASE),
        scripts_sha256={p.name: sha256(p) for p in sorted((ROOT / 'scripts').glob('*')) if p.is_file()}))


def ess_task(tissue, size):
    name = tag(tissue, size, 'SpatialESS')
    if (ROOT / 'results' / name / 'complete.json').exists():
        if not validate_success(name):
            raise RuntimeError('Invalid reusable ESS completion: ' + name)
        log('REUSE ESS', name, 'original results read only')
        return True
    return task(tissue, size, 'SpatialESS')


def task(tissue, size, engine):
    name = tag(tissue, size, engine)
    timeout = CFG['reference_timeout_seconds'] if engine == 'SpatialCellChat' else CFG['ess_timeout_seconds']
    command = [CFG['rscript'], str(ROOT / 'scripts/worker.R'), str(ROOT), tissue,
               fixture(size), engine, 'formal', '100', name]
    return run(name, command, timeout, dict(tissue=tissue, fixture=fixture(size),
                                           engine=engine, phase='formal', nboot=100))


def compare(tissue, size):
    ess = tag(tissue, size, 'SpatialESS')
    reference = tag(tissue, size, 'SpatialCellChat')
    return run('compare_' + ess, [CFG['python'], str(ROOT / 'scripts/compare_scale.py'), ess, reference],
               600, dict(tissue=tissue, fixture=fixture(size), operation='comparison'))


def main():
    lock = (ROOT / 'logs/queue.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit('Another scale queue holds the lock')
    setup()
    summarize()
    log('QUEUE START', 'CPU', CPU, '12h reference limit', 'sequential stages, per-tissue stopping')
    if not run('verify_extension', [CFG['rscript'], str(ROOT / 'scripts/prepare_extension.R'), str(ROOT), str(BASE)], 7200, {}):
        raise RuntimeError('Input preflight failed')
    active = set(CFG['tissues'])
    states = {}
    for size in SIZES:
        outcomes = []
        for tissue in CFG['tissues']:
            if tissue not in active:
                continue
            if not ess_task(tissue, size):
                active.remove(tissue)
                states[tissue] = dict(first_uncompleted_size=size, reason='ESS failure')
                outcomes.append(dict(tissue=tissue, size=size, status='ESS_FAILED'))
                continue
            reference_ok = task(tissue, size, 'SpatialCellChat')
            agreement_ok = reference_ok and compare(tissue, size)
            outcomes.append(dict(tissue=tissue, size=size, reference_completed=reference_ok,
                                 strict_agreement=agreement_ok))
            if not reference_ok or not agreement_ok:
                active.remove(tissue)
                states[tissue] = dict(first_stopped_size=size,
                                     reason='reference failure/budget' if not reference_ok else 'agreement discrepancy')
        write_json(ROOT / 'reports' / f'stage_{size}.json', dict(size=size, outcomes=outcomes))
        summarize()
        if not active:
            break
    status = 'COMPLETE_TESTED_PLAN' if len(active) == len(CFG['tissues']) else 'FINISHED_WITH_BUDGET_OR_VALIDATION_STOP'
    write_json(ROOT / 'reports/progress.json', dict(status=status, finished=stamp(),
                                                   stopped_tissues=states, reached_end=list(sorted(active))))
    summarize()
    log('QUEUE FINISHED', status)


if __name__ == '__main__':
    import runpy
    runpy.run_path(str(ROOT / 'scripts/queue_scale_v2.py'), run_name='__main__')
