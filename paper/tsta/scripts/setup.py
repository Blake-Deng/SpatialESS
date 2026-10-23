"""Freeze input manifest, source snapshots and executable/environment provenance."""
import csv
import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "config/design.json").read_text())


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    for name in ("logs", "tables", "reports", "results", "lib", "inputs", "vendor", "run_state"):
        (ROOT / name).mkdir(exist_ok=True)
    assert not (ROOT / "queue_status.json").exists(), "Setup already used; resume queue.py instead."
    inp = Path(CFG["input_root"])
    source = Path(CFG["source_root"])
    destination = ROOT / "vendor/SpatialESS"
    if not destination.exists():
        shutil.copytree(source, destination, ignore=shutil.ignore_patterns(
            ".git", "paper", "benchmarks", "*.o", "*.so", "*.zip", "*.tar.gz", "tests"))
    records = []
    for p in sorted(destination.rglob("*")):
        if p.is_file():
            records.append({"path": str(p.relative_to(destination)), "sha256": sha256(p)})
    with (ROOT / "reports/package_source_sha256.tsv").open("w") as h:
        writer = csv.DictWriter(h, records[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    secact = Path(CFG["secact_root"])
    for name in ("Secretome_s0_ifun.R", "Secretome_s2_coexpr_1_vst.R"):
        shutil.copy2(secact / "code" / name, ROOT / "vendor" / name)
    with (inp / "stat.csv").open() as h:
        studies = {r[""]: r for r in csv.DictReader(h)}
    rows = []
    with (inp / "meta.csv").open() as h:
        meta = list(csv.DictReader(h))
    assert len(meta) == 1258
    for idx, r in enumerate(meta, start=1):
        p = inp / "data" / r["Study"] / r["Sample_Name"] / "counts.tsv.gz"
        assert p.is_file()
        with gzip.open(p, "rt") as h:
            header = h.readline().rstrip("\r\n").split("\t")
        if header[0] in ("", '""'):
            header = header[1:]
        header = [v.strip('"') for v in header]
        valid = all(re.fullmatch(r"\d+x\d+", v) for v in header)
        if valid:
            valid = all(sum(map(int, v.split("x"))) % 2 == 0 for v in header)
        valid = valid and len(set(header)) == len(header)
        rows.append(dict(sample_id=f"TSTA_{idx:04d}", Study=r["Study"],
                         Cancer=r["Cancer"], Sample_Name=r["Sample_Name"],
                         Normal=r["Normal"], source_publication=studies[r["Study"]]["PMID"],
                         counts_path=str(p), compressed_bytes=p.stat().st_size,
                         input_spots=len(header),
                         coordinate_status="VALID_SAVED_ARRAY_GRID" if valid else "NEEDS_REVIEW"))
    assert len({(r['Study'], r['Sample_Name']) for r in rows}) == len(rows)
    with (ROOT / "tables/input_manifest.tsv").open("w") as h:
        writer = csv.DictWriter(h, rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    for name in ("meta.csv", "stat.csv"):
        shutil.copy2(inp / name, ROOT / "inputs" / name)
    env = dict(os.environ)
    env["PATH"] = str(Path(CFG["rscript"]).parent) + os.pathsep + env.get("PATH", "")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[name] = "1"
    env["MAKEFLAGS"] = "-j2"
    cmd = [str(Path(CFG["rscript"]).with_name("R")), "CMD", "INSTALL", "--preclean",
           "--no-multiarch", "--library=" + str(ROOT / "lib"), str(destination)]
    with (ROOT / "logs/package_install.log").open("w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    source_digest = hashlib.sha256(json.dumps(records, sort_keys=True).encode()).hexdigest()
    env_summary = dict(checked_at=datetime.now().astimezone().isoformat(),
                       python=sys.version, platform=platform.platform(),
                       cpu_affinity=sorted(os.sched_getaffinity(0)),
                       config_sha256=sha256(ROOT / "config/design.json"),
                       source_manifest_sha256=source_digest,
                       secact_coordinate_source_sha256=sha256(ROOT / "vendor/Secretome_s0_ifun.R"),
                       samples=len(rows), studies=len(studies),
                       cancer_types=len({r["Cancer"] for r in rows}),
                       input_spots=sum(r["input_spots"] for r in rows),
                       valid_saved_coordinate_grids=sum(r["coordinate_status"] == "VALID_SAVED_ARRAY_GRID" for r in rows))
    (ROOT / "reports/setup.json").write_text(json.dumps(env_summary, indent=2) + "\n")
    print(json.dumps(env_summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
