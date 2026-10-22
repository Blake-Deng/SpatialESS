#!/usr/bin/env python3
"""Run the full 715,512-cell fixed-LR SpatialCellChat reference attempt."""
import json
import os
import re
import resource
import signal
import subprocess
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / 'config/design.json').read_text())
(ROOT / 'tables').mkdir(exist_ok=True)
TAG = 'cervical_full_715512_SpatialCellChat_formal_b100'
EXIT = ROOT / 'logs' / f'{TAG}.exit.json'
LOG = ROOT / 'logs' / f'{TAG}.log'
TIME = ROOT / 'logs' / f'{TAG}.time.txt'
R = '/home/dzf/miniforge3/envs/spatialcellchat-v3/bin/Rscript'
WORKER = ROOT / 'scripts/worker.R'


def stamp():
    return datetime.now().astimezone().isoformat(timespec='seconds')


def write_status(value):
    temp = EXIT.with_suffix('.partial')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(EXIT)


def child_limits():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    limit = int(CFG['address_space_limit_GiB'] * 1024**3)
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    os.sched_setaffinity(0, {CFG['cpu_affinity']})


def main():
    if EXIT.exists():
        raise SystemExit(f'Refusing to overwrite an existing terminal record: {EXIT}')
    env = dict(os.environ, **CFG['thread_environment'], PYTHONUNBUFFERED='1')
    env['PATH'] = str(Path(R).parent) + ':/usr/bin:/bin'
    command = ['/usr/bin/time', '-v', '-o', str(TIME), R, str(WORKER),
               str(ROOT), 'cervical', 'full', 'SpatialCellChat', 'formal', '100', TAG]
    attempt = {
        'tag': TAG, 'status': 'RUNNING', 'started': stamp(),
        'input': str(ROOT / 'inputs/cervical/full.rds'),
        'input_cells': CFG['input_cells'], 'nominal_figure_category': '700k',
        'genes': CFG['genes'], 'lr_records': CFG['lr_records'], 'nboot': CFG['nboot'],
        'seed_use': CFG['seed_use'], 'cpu_affinity': CFG['cpu_affinity'],
        'threads': CFG['thread_environment'],
        'address_space_limit_GiB': CFG['address_space_limit_GiB'],
        'timeout_seconds': CFG['timeout_seconds'], 'command': command,
    }
    with LOG.open('w') as log:
        log.write(f'{stamp()} START {TAG}\n')
        log.flush()
        proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                env=env, preexec_fn=child_limits, start_new_session=True)
        attempt['pid'] = proc.pid
        write_status(attempt)
        started = time.monotonic()
        deadline = started + CFG['timeout_seconds']
        timed_out = False
        while proc.poll() is None:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                timed_out = True
                os.killpg(proc.pid, signal.SIGKILL)
                break
            time.sleep(min(30, remaining))
        code = proc.wait()
    log_text = LOG.read_text(errors='replace')
    if timed_out:
        status = 'timeout'
    elif code == 0 and (ROOT / 'results' / TAG / 'complete.json').exists():
        status = 'completed'
    elif re.search(r'cannot allocate|std::bad_alloc|vector memory exhausted|Cannot allocate memory',
                   log_text, re.I):
        status = 'allocation_failure_under_96_GiB_virtual_memory_budget'
    elif code < 0:
        status = f'terminated_by_signal_{-code}'
    else:
        status = 'failed'
    attempt.update(status=status, returncode=code, finished=stamp(),
                   elapsed_seconds=round(time.monotonic() - started, 3))
    if TIME.exists():
        time_text = TIME.read_text(errors='replace')
        match = re.search(r'Maximum resident set size \(kbytes\): (\d+)', time_text)
        if match:
            attempt['peak_RSS_GB_decimal'] = int(match.group(1)) * 1024 / 1e9
    write_status(attempt)
    with LOG.open('a') as log:
        log.write(f'{stamp()} END {TAG} status={status} returncode={code}\n')
    subprocess.run(['/home/dzf/miniforge3/bin/python',
                    str(ROOT / 'scripts/export_result.py')], check=True)
    print(json.dumps(attempt, indent=2), flush=True)
    return 0 if status == 'completed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
