"""Export only observed outcomes; queued or missing runs have no measurements."""
import csv
import json
from scale_common import CFG, ROOT, SIZES, fixture, read_row, stamp, tag, write_json


def summarize():
    rows = []
    agreements = []
    for size in SIZES:
        for engine in ["SpatialESS", "SpatialCellChat"]:
            name = tag("cervical", size, engine)
            exit_path = ROOT / "logs" / (name + ".exit.json")
            live_path = ROOT / "logs" / (name + ".running.json")
            row = dict(tissue="cervical", fixture=fixture(size), N=size, engine=engine,
                       status="not_attempted", genes=18028, lr=3081, nboot=100, seed_use=1,
                       whole_process_seconds="", inference_seconds="", peak_RSS_GB_decimal="",
                       returncode="", active_records="", input_md5="", tag=name)
            if exit_path.exists():
                attempt = json.loads(exit_path.read_text())
                row.update(status=attempt["status"], returncode=attempt["returncode"])
                row["whole_process_seconds"] = attempt["measurements"].get("whole_process_seconds", "")
                row["peak_RSS_GB_decimal"] = attempt["measurements"].get("peak_RSS_GB_decimal", "")
                if row["status"] == "completed":
                    saved = read_row(ROOT / "results" / name / "summary.tsv")
                    for key in ["inference_seconds", "active_records", "input_md5"]:
                        row[key] = saved[key]
            elif live_path.exists():
                row["status"] = "running"
            rows.append(row)
        comparison = ROOT / "comparisons" / tag("cervical", size, "SpatialESS") / "agreement.tsv"
        if comparison.exists():
            agreements.append(read_row(comparison))
    for filename, records in [("performance.tsv", rows), ("agreement.tsv", agreements)]:
        path = ROOT / "tables" / filename
        tmp = path.with_suffix(".tsv.tmp")
        fields = list(records[0]) if records else ["N", "strict_pass"]
        with tmp.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader(); writer.writerows(records)
        tmp.replace(path)
    write_json(ROOT / "reports/summary.json", dict(updated=stamp(), outcomes=rows,
        comparisons=agreements, memory_budget=CFG["memory_budget_definition"],
        figure_snapshot_policy=CFG["figure_policy"]))
    lines = ["# Cervical 200k / 400k supplemental runs", "", "Updated: " + stamp(), "",
             "| N | Engine | Status | Whole-process seconds | Peak RSS (decimal GB) |",
             "|---|---|---|---|---|"]
    for row in rows:
        lines.append("| {N} | {engine} | {status} | {whole_process_seconds} | {peak_RSS_GB_decimal} |".format(**row))
    lines += ["", "Inputs: 18,028 genes, 3,081 fixed LR records, 100 permutations, seed 1.",
              "Both methods use the identical saved input and one CPU/thread, sequentially.",
              "Per-worker budget: 96 GiB virtual address space; SpatialCellChat 12 h, SpatialESS 24 h.",
              "Resource measurements for unsuccessful attempts are failure observations, not successful runtime results.",
              "Compare numerical agreement only if both methods complete. No overexpression-based LR selection is run.",
              "The delivered figure remains a pre-run snapshot with 200k/400k marked not attempted.", ""]
    (ROOT / "reports/results_summary.md").write_text("\n".join(lines))


if __name__ == "__main__":
    summarize()
