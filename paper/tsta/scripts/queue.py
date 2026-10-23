"""Persistent, bounded parallel queue. Resume validated successes; retain failures."""
import argparse
import csv
import fcntl
import hashlib
import json
import os
import platform
import signal
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "config/design.json").read_text())
STOP = False


def timestamp():
    return datetime.now().astimezone().isoformat(timespec="seconds")


def atomic_json(path, obj):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")
    tmp.replace(path)


def event(name, **kwargs):
    row = dict(time=timestamp(), event=name, **kwargs)
    with (ROOT / "logs/queue.jsonl").open("a") as h:
        h.write(json.dumps(row, ensure_ascii=False) + "\n")
    print(json.dumps(row, ensure_ascii=False), flush=True)


def read_json(path, default=None):
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, ValueError):
        return default


def time_metrics(path):
    vals = {}
    if not path.exists():
        return vals
    fields = {"Elapsed (wall clock) time": "whole_process_seconds",
              "Maximum resident set size (kbytes)": "peak_rss_kib"}
    for line in path.read_text().splitlines():
        for prefix, key in fields.items():
            if line.lstrip().startswith(prefix):
                value = line.rsplit(": ", 1)[-1].strip()
                if key == "whole_process_seconds":
                    parts = value.split(":")
                    vals[key] = sum(float(v)*60**i for i, v in enumerate(reversed(parts)))
                else:
                    vals[key] = int(value)
    if "peak_rss_kib" in vals:
        vals["peak_rss_GB_decimal"] = vals["peak_rss_kib"]*1024/1e9
    return vals


def write_tsv(path, rows, required):
    fields = list(required)
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as h:
        writer = csv.DictWriter(h, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def report(manifest, running, final=False):
    states = []; perf = []
    for item in manifest:
        sid = item["sample_id"]
        runfile = ROOT / "run_state" / (sid + ".json")
        rec = read_json(runfile, {})
        progress = read_json(ROOT / "results" / sid / "progress.json", {})
        row = {key: item[key] for key in ("sample_id", "Study", "Cancer", "Sample_Name", "input_spots")}
        row.update(status=rec.get("status", "PENDING"), stage=progress.get("stage", ""),
                   returncode=rec.get("returncode", ""), detail=rec.get("detail", ""),
                   pid=rec.get("pid", ""), **rec.get("time_metrics", {}))
        if sid in running:
            row["status"] = "RUNNING"
        states.append(row)
        summary = ROOT / "results" / sid / "summary.tsv"
        if row["status"] == "COMPLETE" and summary.exists():
            with summary.open() as h:
                r = next(csv.DictReader(h, delimiter="\t"))
            r.update(rec.get("time_metrics", {}))
            worker_limit = rec.get("worker_limit", 1 if sid == "TSTA_0001" else CFG["workers"])
            r["concurrency_limit"] = worker_limit
            r["measurement_scope"] = f"{worker_limit} concurrent slots; single thread per sample"
            perf.append(r)
    counts = dict(Counter(r["status"] for r in states))
    snapshot = dict(updated=timestamp(), supervisor_pid=os.getpid(),
                    status="RUNNING" if not final else "FINISHED_WITH_FAILURES" if any(
                        k not in ("COMPLETE",) and n for k, n in counts.items()) else "COMPLETE",
                    total=len(manifest), counts=counts, workers_configured=CFG["workers"],
                    active=list(running), completed_spots=sum(int(r["N"]) for r in perf),
                    walltime_scope="concurrent throughput, not isolated reference-paired benchmark")
    atomic_json(ROOT / "queue_status.json", snapshot)
    write_tsv(ROOT / "tables/sample_status.tsv", states,
              ["sample_id", "Study", "Cancer", "Sample_Name", "input_spots", "status", "stage"])
    write_tsv(ROOT / "tables/tsta_performance.tsv", perf,
              ["sample_id", "Study", "Cancer", "Sample_Name", "N", "inference_seconds",
               "whole_process_seconds", "peak_rss_GB_decimal"])
    with (ROOT / "reports/RUN_SUMMARY.md").open("w") as h:
        h.write("# TSTA SpatialESS run\n\nUpdated: " + timestamp() + "\n\n")
        h.write(f"All {len(manifest)} samples are included in the queue. Status: {counts}.\n\n")
        h.write(f"Completed retained spots: {snapshot['completed_spots']:,}.\n\n")
        h.write("Saved array IDs are converted with the archived SecAct calWeights() formula. "
                "These are conventional Visium spots, not segmented single cells. "
                "Groups are per-sample expression graph clusters, not cell-type annotations. "
                "All complete-coverage human database LR records are tested with 100 global "
                "label permutations and a 200 um spot radius; cell-contact restriction is disabled.\n\n")
        h.write("Each worker is CPU-pinned and single-threaded. GNU time records complete process "
                "wall time and peak RSS; inference time is separately measured inside R. "
                "Sixteen concurrent sample workers measure batch throughput under contention, "
                "not isolated speedup against SpatialCellChat. No reference-equivalence or new "
                "cell-type-specific biological finding is established by this feasibility run.\n")
    return snapshot


def terminate(signum, frame):
    global STOP
    STOP = True
    event("STOP_REQUESTED", signal=signum)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--retry-failed", action="store_true")
    args = parser.parse_args()
    lock = (ROOT / "run_state/queue.lock").open("a")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    with (ROOT / "tables/input_manifest.tsv").open() as h:
        manifest = list(csv.DictReader(h, delimiter="\t"))
    if args.limit:
        manifest = manifest[:args.limit]
    assert (ROOT / "inputs/database.rds").exists()
    scripts_digest = {}
    for path in sorted((ROOT / "scripts").glob("*")):
        if path.is_file():
            scripts_digest[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    atomic_json(ROOT / "reports" / f"runner_sources_{int(time.time())}.json", scripts_digest)
    for signum in (signal.SIGTERM, signal.SIGINT):
        signal.signal(signum, terminate)
    affinity = set(os.sched_getaffinity(0))
    cpus = [c for c in CFG["cpus"] if c in affinity]
    nworkers = args.workers or CFG["workers"]
    assert len(cpus) >= nworkers
    env = dict(os.environ)
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[k] = "1"
    env["R_LIBS_USER"] = str(ROOT / "lib")
    pending = []
    for row in manifest:
        sid = row["sample_id"]
        old = read_json(ROOT / "run_state" / (sid + ".json"), {})
        complete = ROOT / "results" / sid / "complete.json"
        if old.get("status") == "COMPLETE" and complete.exists():
            continue
        if old.get("status") in ("FAILED", "TIMEOUT", "INTERRUPTED", "EXCLUDED") and not args.retry_failed:
            continue
        # A stale RUNNING worker must not be launched a second time.
        if old.get("status") == "RUNNING":
            cmdline = Path("/proc") / str(old.get("pid", -1)) / "cmdline"
            if cmdline.exists() and str(ROOT).encode() in cmdline.read_bytes():
                raise RuntimeError(f"Worker {sid} is still alive; cannot resume into duplicate inference")
        if complete.exists() and old.get("status") != "COMPLETE":
            raise RuntimeError(f"{sid} has a result without verified process metrics; inspect before resuming")
        pending.append(row)
    running = {}; free_cpus = cpus[:nworkers]
    event("QUEUE_START", pending=len(pending), workers=nworkers, cpus=free_cpus,
          package=CFG["package_version"], nboot=CFG["nboot"])
    last_report = 0
    while pending or running:
        while pending and free_cpus and not STOP:
            row = pending.pop(0); sid = row["sample_id"]; cpu = free_cpus.pop(0)
            previous = read_json(ROOT / "run_state" / (sid + ".json"), {})
            if previous:
                previous_attempt = previous.get("started", timestamp()).replace(":", "")
                atomic_json(ROOT / "logs" / f"{sid}.{previous_attempt}.previous_state.json", previous)
            if row["coordinate_status"] != "VALID_SAVED_ARRAY_GRID":
                atomic_json(ROOT / "run_state" / (sid + ".json"),
                            dict(status="EXCLUDED", detail="Unvalidated saved spatial array identifiers"))
                free_cpus.append(cpu)
                continue
            attempt = int(time.time())
            logpath = ROOT / "logs" / f"{sid}.{attempt}.log"
            timepath = ROOT / "logs" / f"{sid}.{attempt}.time.txt"
            cmd = ["/usr/bin/time", "-v", "-o", str(timepath),
                   "/usr/bin/timeout", "--signal=TERM", "--kill-after=30s",
                   str(CFG["worker_timeout_seconds"]) + "s", "/usr/bin/prlimit",
                   "--as=" + str(CFG["worker_address_space_GiB"]*1024**3),
                   "--core=0", "/usr/bin/taskset", "-c", str(cpu),
                   CFG["rscript"], "--vanilla", str(ROOT / "scripts/worker.R"), str(ROOT), sid]
            log = logpath.open("w")
            proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                                    env=env, start_new_session=True)
            state = dict(status="RUNNING", sample_id=sid, pid=proc.pid, cpu=cpu, worker_limit=nworkers,
                         started=timestamp(), command=cmd, log=str(logpath),
                         time_log=str(timepath))
            atomic_json(ROOT / "run_state" / (sid + ".json"), state)
            running[sid] = dict(proc=proc, log=log, cpu=cpu, state=state,
                                timepath=timepath, started=time.monotonic())
            event("WORKER_START", sample_id=sid, pid=proc.pid, cpu=cpu,
                  study=row["Study"], sample=row["Sample_Name"])
        for sid, run in list(running.items()):
            rc = run["proc"].poll()
            if rc is None:
                continue
            run["log"].close()
            out = ROOT / "results" / sid
            success = rc == 0 and all((out / n).exists() for n in (
                "complete.json", "summary.tsv", "records.tsv.gz", "SpatialESS_result.rds"))
            state = run["state"]
            state.update(status="COMPLETE" if success else "TIMEOUT" if rc == 124 else "FAILED",
                         returncode=rc, finished=timestamp(),
                         supervisor_elapsed_seconds=time.monotonic()-run["started"],
                         time_metrics=time_metrics(run["timepath"]))
            if rc == 0 and not success:
                state["detail"] = "Worker returned zero but required output files are missing"
            if not success:
                lines = Path(state["log"]).read_text(errors="replace").splitlines()
                state["failure_tail"] = lines[-15:]
            atomic_json(ROOT / "run_state" / (sid + ".json"), state)
            free_cpus.append(run["cpu"])
            del running[sid]
            event("WORKER_END", sample_id=sid, status=state["status"], returncode=rc,
                  **state["time_metrics"])
        if time.monotonic()-last_report > 15:
            report(manifest, running)
            last_report = time.monotonic()
        if STOP:
            # Finish already-started samples, then exit without starting new ones.
            pending = []
        if pending or running:
            time.sleep(2)
    snapshot = report(manifest, running, final=True)
    event("QUEUE_END", counts=snapshot["counts"])
    return 0 if snapshot["counts"].get("FAILED", 0) + snapshot["counts"].get("TIMEOUT", 0) == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
