"""Read persisted queue status; report live processes rather than a stale tmux name."""
import json
import os
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'reports/progress.json'
if path.exists():
    row=json.loads(path.read_text())
    print(json.dumps(row,indent=2))
    for field in ['runner_pid','child_pid']:
        if field in row:
            try:
                os.kill(int(row[field]),0)
                print(field,'LIVE')
            except ProcessLookupError:
                print(field,'NOT LIVE')
    if row.get('log') and Path(row['log']).exists():
        print('\nCurrent task log tail:\n'+'\n'.join(Path(row['log']).read_text(errors='replace').splitlines()[-15:]))
else:
    print('No progress checkpoint yet')
if (ROOT/'reports/table_status.json').exists():
    print('\nTable status:\n'+(ROOT/'reports/table_status.json').read_text())
