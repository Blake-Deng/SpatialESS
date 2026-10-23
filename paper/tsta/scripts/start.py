"""Run one real sample and graph oracle before starting the full parallel queue."""
import json
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
cfg = json.loads((root / "config/design.json").read_text())


def run(args, name, allow_sample_failures=False):
    print("START", name, flush=True)
    with (root / "logs" / name).open("a") as log:
        result = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, cwd=root)
    print("END", name, "returncode", result.returncode, flush=True)
    if result.returncode and not (allow_sample_failures and result.returncode == 2):
        raise SystemExit(result.returncode)
    return result.returncode


if not (root / "inputs/database.rds").exists():
    run([cfg["rscript"], "--vanilla", str(root / "scripts/freeze_database.R"), str(root)], "database_freeze.log")
if not (root / "reports/pilot_validation.json").exists():
    pilot = [sys.executable, str(root / "scripts/queue.py"), "--limit", "1", "--workers", "1"]
    if "--retry-pilot" in sys.argv:
        pilot.append("--retry-failed")
    run(pilot, "pilot_queue.log")
    run([cfg["rscript"], "--vanilla", str(root / "scripts/validate_pilot.R"), str(root)], "pilot_validation.log")
queue_returncode = run([sys.executable, str(root / "scripts/queue.py")], "queue_console.log",
                      allow_sample_failures=True)
run([sys.executable, str(root / "scripts/finalize.py")], "finalize.log")
raise SystemExit(queue_returncode)
