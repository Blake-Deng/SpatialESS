#!/usr/bin/env python3
"""Prepare matched RAW and SPARKLE inputs for SpatialESS."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


ROOT = Path("/home/dzf/cellchat_acceleration/SpatialESS_RNA_leakage_20260902")
SPARKLE_ROOT = ROOT / "software" / "SPARKLE"
H5 = SPARKLE_ROOT / "evaluation" / "data" / "ovarian" / (
    "Visium_HD_Human_Ovarian_Cancer_FF_feature_slice.h5"
)
H5AD_ROOT = SPARKLE_ROOT / "evaluation" / "reports" / "h5ad"
TAG = "ovarian_x1000-1800_y300-1100"
RAW_H5AD = H5AD_ROOT / f"{TAG}_raw.h5ad"
SPARKLE_H5AD = H5AD_ROOT / f"{TAG}_SPARKLE.h5ad"
GENE_LIST = ROOT / "inputs" / "cellchatdb_human_genes.txt"
OUT = ROOT / "inputs" / "spatialess_matched"
X_RANGE = (1000, 1800)
Y_RANGE = (300, 1100)
CHUNK_SIZE = 256


def library_sizes(x) -> np.ndarray:
    totals = np.zeros(x.shape[0], dtype=np.float64)
    for start in range(0, x.shape[0], CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, x.shape[0])
        totals[start:stop] = np.asarray(x[start:stop, :], dtype=np.float64).sum(axis=1)
    return totals


def normalized_selected(x, indices: np.ndarray, totals: np.ndarray) -> sparse.csc_matrix:
    selected = np.asarray(x[:, indices], dtype=np.float64)
    selected *= (1e4 / np.maximum(totals, 1.0))[:, None]
    np.log1p(selected, out=selected)
    selected[selected == 0] = 0
    return sparse.csc_matrix(selected.T)


def write_matrix(path: Path, matrix: sparse.spmatrix) -> None:
    with gzip.open(path, "wb", compresslevel=5) as handle:
        mmwrite(handle, matrix, precision=16)


def derive_spatial_metadata(cell_ids: np.ndarray) -> tuple[pd.DataFrame, dict]:
    with h5py.File(H5, "r") as handle:
        seg = handle["segmentations/cell_segmentation_mask"]
        rows = seg["row"][:].astype(np.int32, copy=False)
        cols = seg["col"][:].astype(np.int32, copy=False)
        labels = seg["data"][:].astype(np.int64, copy=False)
        keep = (
            (cols >= X_RANGE[0]) & (cols <= X_RANGE[1])
            & (rows >= Y_RANGE[0]) & (rows <= Y_RANGE[1])
        )
        rows = rows[keep]
        cols = cols[keep]
        labels = labels[keep]

        cluster = handle[
            "secondary_analysis/clustering/"
            "square_008um_gene_expression_graphclust"
        ]
        cluster_rows = cluster["row"][:].astype(np.int32, copy=False)
        cluster_cols = cluster["col"][:].astype(np.int32, copy=False)
        cluster_values = cluster["data"][:].astype(np.int16, copy=False)

    max_cluster = int(cluster_values.max())
    cluster_map = np.zeros(
        (int(cluster_rows.max()) + 1, int(cluster_cols.max()) + 1),
        dtype=np.int16,
    )
    cluster_map[cluster_rows, cluster_cols] = cluster_values
    pixel_clusters = cluster_map[rows // 4, cols // 4]

    id_to_position = {int(value): i for i, value in enumerate(cell_ids)}
    positions = np.fromiter(
        (id_to_position.get(int(value), -1) for value in labels),
        dtype=np.int32,
        count=len(labels),
    )
    matched = positions >= 0
    positions = positions[matched]
    rows = rows[matched]
    cols = cols[matched]
    pixel_clusters = pixel_clusters[matched]

    counts = np.bincount(positions, minlength=len(cell_ids)).astype(np.int64)
    x_sum = np.bincount(
        positions, weights=cols.astype(np.float64) * 2.0 + 1.0,
        minlength=len(cell_ids),
    )
    y_sum = np.bincount(
        positions, weights=rows.astype(np.float64) * 2.0 + 1.0,
        minlength=len(cell_ids),
    )
    if np.any(counts == 0):
        missing = cell_ids[counts == 0]
        raise ValueError(f"{len(missing)} h5ad cells lack segmentation pixels")

    votes = np.zeros((len(cell_ids), max_cluster + 1), dtype=np.int32)
    informative = pixel_clusters > 0
    np.add.at(votes, (positions[informative], pixel_clusters[informative]), 1)
    groups = votes[:, 1:].argmax(axis=1) + 1
    vote_max = votes[:, 1:].max(axis=1)
    groups[vote_max == 0] = 0

    frame = pd.DataFrame({
        "cell": [f"Cell_{int(value)}" for value in cell_ids],
        "cell_id": cell_ids.astype(np.int64),
        "x_um": x_sum / counts,
        "y_um": y_sum / counts,
        "group": np.where(groups > 0, [f"RAW_cluster_{x:02d}" for x in groups], "Unknown"),
        "segmentation_pixels": counts,
        "cluster_vote_pixels": vote_max,
        "cluster_vote_fraction": vote_max / counts,
    })
    summary = {
        "n_cells": int(len(cell_ids)),
        "n_raw_graph_clusters": max_cluster,
        "n_unknown": int(np.sum(groups == 0)),
        "median_cluster_vote_fraction": float(np.median(frame["cluster_vote_fraction"])),
        "min_cluster_vote_fraction": float(frame["cluster_vote_fraction"].min()),
    }
    return frame, summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    requested = {line.strip() for line in GENE_LIST.read_text().splitlines() if line.strip()}

    raw = ad.read_h5ad(RAW_H5AD, backed="r")
    corrected = ad.read_h5ad(SPARKLE_H5AD, backed="r")
    try:
        if raw.shape != corrected.shape:
            raise ValueError("RAW and SPARKLE h5ad shapes differ")
        if not raw.obs_names.equals(corrected.obs_names):
            raise ValueError("RAW and SPARKLE cell order differs")
        if not raw.var_names.equals(corrected.var_names):
            raise ValueError("RAW and SPARKLE gene order differs")

        indices = np.array(
            [i for i, gene in enumerate(raw.var_names) if gene in requested],
            dtype=np.int64,
        )
        genes = raw.var_names[indices].astype(str).to_numpy()
        cell_ids = raw.obs["cell_id"].to_numpy(dtype=np.int64)
        raw_totals = library_sizes(raw.X)
        corrected_totals = library_sizes(corrected.X)
        raw_matrix = normalized_selected(raw.X, indices, raw_totals)
        corrected_matrix = normalized_selected(corrected.X, indices, corrected_totals)
    finally:
        raw.file.close()
        corrected.file.close()

    metadata, spatial_summary = derive_spatial_metadata(cell_ids)
    cells = metadata["cell"].to_numpy()
    pd.Series(genes).to_csv(OUT / "genes.tsv", index=False, header=False)
    pd.Series(cells).to_csv(OUT / "cells.tsv", index=False, header=False)
    metadata.to_csv(OUT / "cell_metadata.tsv", sep="\t", index=False)
    write_matrix(OUT / "raw_log1p_cp10k.mtx.gz", raw_matrix)
    write_matrix(OUT / "sparkle_log1p_cp10k.mtx.gz", corrected_matrix)

    summary = {
        **spatial_summary,
        "n_cellchat_genes_requested": len(requested),
        "n_cellchat_genes_present": int(len(genes)),
        "raw_total_counts": float(raw_totals.sum()),
        "sparkle_total_counts": float(corrected_totals.sum()),
        "count_retained_fraction": float(corrected_totals.sum() / raw_totals.sum()),
        "normalization": "per-cell log1p(CP10K), denominator uses all 18,786 genes",
        "group_definition": "majority 10x RAW square_008um graph cluster over each segmented cell",
    }
    (OUT / "preparation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
