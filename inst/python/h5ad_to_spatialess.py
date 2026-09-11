#!/usr/bin/env python3
"""Export a SPARKLE H5AD matrix as a validated SpatialESS bundle."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert H5AD expression and spatial metadata for SpatialESS."
    )
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--group-col", required=True)
    parser.add_argument("--coordinate-cols", required=True)
    parser.add_argument("--gene-list", type=Path)
    parser.add_argument("--metadata-tsv", type=Path)
    parser.add_argument("--cell-col", default="cell")
    parser.add_argument("--layer")
    parser.add_argument(
        "--normalization",
        choices=("log1p_cp10k", "none"),
        default="log1p_cp10k",
    )
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def aligned_metadata(
    adata: ad.AnnData,
    metadata_tsv: Path | None,
    cell_col: str,
) -> pd.DataFrame:
    cells = adata.obs_names.astype(str)
    if metadata_tsv is None:
        metadata = adata.obs.copy()
        metadata.index = cells
        return metadata

    metadata = pd.read_csv(metadata_tsv, sep="\t", dtype={cell_col: str})
    if cell_col not in metadata:
        raise ValueError(f"metadata TSV lacks cell column: {cell_col}")
    if metadata[cell_col].duplicated().any():
        raise ValueError("metadata TSV contains duplicate cell identifiers")
    metadata = metadata.set_index(cell_col, drop=False)
    missing = cells.difference(metadata.index)
    if len(missing):
        raise ValueError(
            f"metadata TSV lacks {len(missing)} H5AD cells; first={missing[0]}"
        )
    return metadata.loc[cells].copy()


def select_genes(var_names: pd.Index, gene_list: Path | None) -> tuple[np.ndarray, list[str], list[str]]:
    if var_names.has_duplicates:
        raise ValueError("H5AD var_names must be unique")
    if gene_list is None:
        indices = np.arange(len(var_names), dtype=np.int64)
        return indices, var_names.astype(str).tolist(), []

    requested = [
        line.strip() for line in gene_list.read_text().splitlines() if line.strip()
    ]
    requested_set = set(requested)
    indices = np.flatnonzero(var_names.astype(str).isin(requested_set))
    selected = var_names[indices].astype(str).tolist()
    missing = sorted(requested_set.difference(selected))
    if not selected:
        raise ValueError("none of the requested genes are present in the H5AD")
    return indices.astype(np.int64), selected, missing


def export_matrix(
    matrix,
    selected_indices: np.ndarray,
    normalization: str,
    chunk_size: int,
) -> tuple[sparse.csc_matrix, float, float]:
    chunks = []
    input_total = 0.0
    output_total = 0.0
    for start in range(0, matrix.shape[0], chunk_size):
        stop = min(start + chunk_size, matrix.shape[0])
        block = matrix[start:stop, :]
        if sparse.issparse(block):
            if block.data.size and (
                not np.isfinite(block.data).all() or np.min(block.data) < 0
            ):
                raise ValueError("H5AD expression must be finite and non-negative")
            totals = np.asarray(block.sum(axis=1), dtype=np.float64).ravel()
            selected = block[:, selected_indices].astype(np.float64).toarray()
        else:
            block = np.asarray(block, dtype=np.float64)
            if not np.isfinite(block).all() or np.min(block) < 0:
                raise ValueError("H5AD expression must be finite and non-negative")
            totals = block.sum(axis=1)
            selected = block[:, selected_indices]

        input_total += float(totals.sum())
        if normalization == "log1p_cp10k":
            selected *= (10000.0 / np.maximum(totals, 1.0))[:, None]
            np.log1p(selected, out=selected)
        output_total += float(selected.sum())
        selected[selected == 0] = 0
        chunks.append(sparse.csr_matrix(selected))

    cell_by_gene = sparse.vstack(chunks, format="csr")
    gene_by_cell = cell_by_gene.transpose().tocsc()
    gene_by_cell.eliminate_zeros()
    return gene_by_cell, input_total, output_total


def main() -> None:
    args = parse_args()
    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be positive")
    coordinate_cols = [value.strip() for value in args.coordinate_cols.split(",")]
    if len(coordinate_cols) not in (2, 3) or any(not x for x in coordinate_cols):
        raise ValueError("--coordinate-cols requires two or three names")

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        if not args.overwrite:
            raise FileExistsError(
                f"output directory is not empty: {args.output_dir}"
            )
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad, backed="r")
    try:
        if adata.obs_names.has_duplicates:
            raise ValueError("H5AD obs_names must be unique")
        matrix = adata.X if args.layer is None else adata.layers[args.layer]
        input_genes = int(matrix.shape[1])
        metadata = aligned_metadata(adata, args.metadata_tsv, args.cell_col)
        required = [args.group_col, *coordinate_cols]
        missing_columns = [name for name in required if name not in metadata]
        if missing_columns:
            raise ValueError(f"metadata lacks columns: {missing_columns}")
        if metadata[args.group_col].isna().any():
            raise ValueError("group labels contain missing values")
        coordinates = metadata[coordinate_cols].apply(
            pd.to_numeric, errors="raise"
        )
        if not np.isfinite(coordinates.to_numpy()).all():
            raise ValueError("coordinates contain non-finite values")

        indices, genes, missing_genes = select_genes(
            adata.var_names, args.gene_list
        )
        expression, input_total, output_total = export_matrix(
            matrix, indices, args.normalization, args.chunk_size
        )
        cells = adata.obs_names.astype(str).tolist()
    finally:
        adata.file.close()

    matrix_path = args.output_dir / "expression.mtx.gz"
    with gzip.open(matrix_path, "wb", compresslevel=5) as handle:
        mmwrite(handle, expression, precision=17, symmetry="general")
    (args.output_dir / "genes.tsv").write_text("\n".join(genes) + "\n")
    (args.output_dir / "cells.tsv").write_text("\n".join(cells) + "\n")

    exported_metadata = pd.DataFrame({
        "cell": cells,
        "group": metadata[args.group_col].astype(str).to_numpy(),
    })
    for index, name in enumerate(coordinate_cols, start=1):
        exported_metadata[f"coordinate_{index}"] = coordinates[name].to_numpy()
    exported_metadata.to_csv(
        args.output_dir / "metadata.tsv", sep="\t", index=False
    )

    outputs = ["expression.mtx.gz", "genes.tsv", "cells.tsv", "metadata.tsv"]
    manifest = {
        "schema_version": 1,
        "source_h5ad": str(args.h5ad.resolve()),
        "layer": args.layer or "X",
        "normalization": args.normalization,
        "normalization_denominator": (
            "all H5AD genes per cell"
            if args.normalization == "log1p_cp10k"
            else "not applicable"
        ),
        "cells": len(cells),
        "input_genes": input_genes,
        "exported_genes": len(genes),
        "missing_requested_genes": missing_genes,
        "groups": int(exported_metadata["group"].nunique()),
        "coordinate_columns": coordinate_cols,
        "input_expression_sum": input_total,
        "exported_expression_sum": output_total,
        "nnz": int(expression.nnz),
        "sha256": {
            name: sha256(args.output_dir / name) for name in outputs
        },
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
