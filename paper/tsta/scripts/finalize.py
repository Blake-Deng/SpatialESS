"""Summarize actual completed runs; preserve the full denominator and failures."""
import csv
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def summarize_groups(performance, status, key):
    inputs = defaultdict(list)
    completed = defaultdict(list)
    for row in status:
        inputs[row[key]].append(row)
    for row in performance:
        completed[row[key]].append(row)
    output = []
    for name, items in sorted(inputs.items()):
        runs = completed[name]
        row = {key: name, "input_samples": len(items), "completed_samples": len(runs),
               "failed_samples": sum(r['status'] in ('FAILED', 'TIMEOUT') for r in items),
               "excluded_samples": sum(r['status'] == 'EXCLUDED' for r in items),
               "pending_samples": sum(r['status'] in ('PENDING', 'RUNNING') for r in items)}
        for field in ("N", "genes", "lr", "inference_seconds", "whole_process_seconds", "peak_rss_GB_decimal"):
            vals = [float(r[field]) for r in runs if r.get(field, "") != ""]
            if vals:
                row.update({field+"_min": min(vals), field+"_median": median(vals), field+"_max": max(vals)})
        output.append(row)
    return output


def write_tsv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w") as h:
        writer = csv.DictWriter(h, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    performance = read_tsv(ROOT / "tables/tsta_performance.tsv")
    status = read_tsv(ROOT / "tables/sample_status.tsv")
    manifest = read_tsv(ROOT / "tables/input_manifest.tsv")
    for row in performance:
        folder = ROOT / "results" / row["sample_id"]
        complete = json.loads((folder / "complete.json").read_text())
        assert int(row["N"]) == complete["N"] and int(row["E"]) == complete["E"]
        assert row["package_version"] == "0.1.4" and int(row["nboot"]) == 100
        assert (folder / "graph_validation.txt").exists()
    assert len(manifest) == len(status) == 1258
    by_cancer = summarize_groups(performance, status, "Cancer")
    by_study = summarize_groups(performance, status, "Study")
    write_tsv(ROOT / "tables/cancer_summary.tsv", by_cancer)
    write_tsv(ROOT / "tables/study_summary.tsv", by_study)
    check = dict(checked_at=datetime.now().astimezone().isoformat(),
                 denominator_samples=len(manifest), completed=len(performance),
                 all_completed_sample_metadata_checks_pass=True,
                 incomplete=[r for r in status if r["status"] != "COMPLETE"],
                 inference_comparison="not performed; SpatialESS-only feasibility run",
                 spatial_unit="multicellular Visium spots; per-sample expression clusters")
    (ROOT / "reports/final_validation.json").write_text(json.dumps(check, indent=2) + "\n")
    import pandas as pd
    with pd.ExcelWriter(ROOT / "tables/TSTA_SpatialESS_results.xlsx", engine="openpyxl") as book:
        for name, rows in (("input_manifest", manifest), ("sample_status", status),
                           ("performance", performance), ("by_cancer", by_cancer),
                           ("by_study", by_study)):
            pd.DataFrame(rows).to_excel(book, sheet_name=name, index=False)
            sheet = book.sheets[name]
            sheet.freeze_panes = "A2"
            sheet.auto_filter.ref = sheet.dimensions
            for col in sheet.columns:
                width = max(len(str(c.value or "")) for c in col[:min(len(col), 50)])
                sheet.column_dimensions[col[0].column_letter].width = min(48, max(12, width+2))
    print(json.dumps({"samples": len(manifest), "completed": len(performance),
                      "incomplete": len(check['incomplete']),
                      "workbook": str(ROOT / "tables/TSTA_SpatialESS_results.xlsx")}, indent=2))


if __name__ == "__main__":
    main()
