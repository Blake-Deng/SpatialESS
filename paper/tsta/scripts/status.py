import json
import os
from pathlib import Path

root = Path(__file__).resolve().parents[1]
path = root / "queue_status.json"
if not path.exists():
    print("No queue status yet. Inspect logs/setup_console.log and logs/queue_console.log.")
else:
    obj = json.loads(path.read_text())
    pid = obj.get("supervisor_pid")
    cmd = Path("/proc") / str(pid) / "cmdline"
    obj["supervisor_alive_now"] = cmd.exists() and str(root).encode() in cmd.read_bytes()
    print(json.dumps(obj, indent=2, ensure_ascii=False))
    for sid in obj.get("active", []):
        p = root / "results" / sid / "progress.json"
        if p.exists():
            print(sid, json.loads(p.read_text()))
