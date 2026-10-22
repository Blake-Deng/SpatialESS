"""Transfer the active measurement to the fixed queue without killing the worker."""
import json
import os
import signal
import sys
import time
from datetime import datetime
from pathlib import Path

from scale_common import ROOT, process_alive, stamp, write_json

progress = json.loads((ROOT / 'reports/progress.json').read_text())
if progress.get('status') != 'RUNNING' or not progress.get('runner_pid'):
    raise SystemExit('No live task eligible for transfer')
runner, child = int(progress['runner_pid']), int(progress['child_pid'])
runner_command = Path(f'/proc/{runner}/cmdline').read_bytes().decode(errors='replace')
child_command = Path(f'/proc/{child}/cmdline').read_bytes().decode(errors='replace')
if 'queue_scale.py' not in runner_command or progress['tag'] not in child_command:
    raise SystemExit('PID identity does not match the planned measurement')
children = Path(f'/proc/{child}/task/{child}/children').read_text().split()
if len(children) != 1 or not process_alive(int(children[0])):
    raise SystemExit('Expected active R process not found')
r_pid = int(children[0])
progress.update(started_epoch=datetime.strptime(progress['started'], '%Y-%m-%dT%H:%M:%S%z').timestamp(),
                r_pid=r_pid, status='RUNNING', handoff_timestamp=stamp())
write_json(ROOT / 'logs' / (progress['tag'] + '.running.json'), progress)
write_json(ROOT / 'reports/queue_handoff.json', dict(old_runner=runner, active_time_pid=child,
                                                    active_r_pid=r_pid, tag=progress['tag'],
                                                    started=progress['started'], time=stamp()))
os.kill(runner, signal.SIGTERM)
for _ in range(100):
    if not process_alive(runner):
        break
    time.sleep(.1)
else:
    raise SystemExit('Previous queue did not release control')
if not process_alive(child) or not process_alive(r_pid):
    raise SystemExit('R measurement exited during transfer; checkpoint remains recoverable')
os.execv(sys.executable, [sys.executable, '-u', str(ROOT / 'scripts/queue_scale_v2.py')])
