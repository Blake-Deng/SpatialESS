"""Live task PID, CPU time, memory and budget; tmux alone is not proof of work."""
import json
import subprocess
from scale_common import ROOT, process_alive

path = ROOT / 'reports/progress.json'
if path.exists():
    row = json.loads(path.read_text())
    print(json.dumps(row, indent=2))
    pids = [row[field] for field in ['runner_pid', 'child_pid', 'r_pid'] if row.get(field)]
    print('Live PIDs:', {pid: process_alive(pid) for pid in pids})
    if pids:
        subprocess.run(['ps', '-p', ','.join(str(p) for p in pids), '-o', 'pid,ppid,etime,time,pcpu,rss,stat,comm'])
    if row.get('tag'):
        log = ROOT / 'logs' / (row['tag'] + '.log')
        if log.exists():
            print('\nLog tail:\n' + '\n'.join(log.read_text(errors='replace').splitlines()[-10:]))
if (ROOT / 'reports/table_status.json').exists():
    print('\nTables:\n' + (ROOT / 'reports/table_status.json').read_text())
