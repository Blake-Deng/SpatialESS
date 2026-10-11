#!/usr/bin/env python3
"""Create matched SpatialESS inputs using one shared RAW-derived RCTD label set."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path("/home/dzf/cellchat_acceleration/SpatialESS_RNA_leakage_20260902")
BASE = ROOT / "inputs" / "spatialess_matched"
OUT = ROOT / "inputs" / "spatialess_rctd"
RCTD = (
    ROOT / "software" / "SPARKLE" / "evaluation" / "reports"
    / "rctd_ovarian" / "first_type" / "RAW_first_type.csv"
)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(BASE / "cell_metadata.tsv", sep="\t")
    labels = pd.read_csv(RCTD)
    required = {"cell_id", "first_type"}
    if not required.issubset(labels.columns):
        raise ValueError(f"RCTD file lacks columns: {sorted(required - set(labels.columns))}")
    if labels["cell_id"].duplicated().any():
        raise ValueError("RCTD output contains duplicated cell_id values")

    label_map = labels.set_index("cell_id")["first_type"]
    mapped = metadata["cell_id"].map(label_map)
    metadata["group"] = mapped.fillna("Unknown").astype(str)
    metadata.to_csv(OUT / "cell_metadata.tsv", sep="\t", index=False)

    for name in (
        "genes.tsv",
        "cells.tsv",
        "raw_log1p_cp10k.mtx.gz",
        "sparkle_log1p_cp10k.mtx.gz",
    ):
        destination = OUT / name
        if not destination.exists() and not destination.is_symlink():
            destination.symlink_to(BASE / name)

    counts = metadata["group"].value_counts().sort_index()
    summary = {
        "cells": int(len(metadata)),
        "mapped_cells": int(mapped.notna().sum()),
        "unknown_cells": int(mapped.isna().sum()),
        "mapped_fraction": float(mapped.notna().mean()),
        "cell_types": int((counts.index != "Unknown").sum()),
        "label_source": "RAW RCTD first_type; shared unchanged by RAW and SPARKLE",
        "cell_type_counts": {str(k): int(v) for k, v in counts.items()},
    }
    (OUT / "preparation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

