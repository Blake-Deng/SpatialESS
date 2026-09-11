#!/usr/bin/env python3
"""Create a tiny SPARKLE-like corrected H5AD for the integration demo."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


HERE = Path(__file__).resolve().parent
rng = np.random.default_rng(20260904)
cells = [f"cell_{index:02d}" for index in range(24)]
genes = ["L", "R", "BG1", "BG2"]
group = np.array(["Sender"] * 12 + ["Receiver"] * 12)
counts = rng.poisson(0.2, size=(24, 4)).astype(np.float64)
counts[:12, 0] += rng.poisson(5, size=12)
counts[12:, 1] += rng.poisson(5, size=12)
counts[:, 2:] += rng.poisson(1, size=(24, 2))
counts *= 0.92

obs = pd.DataFrame(
    {
        "annotation": group,
        "x_um": np.tile(np.arange(6) * 10.0, 4),
        "y_um": np.repeat(np.arange(4) * 10.0, 6),
    },
    index=cells,
)
var = pd.DataFrame(index=genes)
adata = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
adata.uns["source"] = "synthetic SPARKLE-like corrected counts"
adata.write_h5ad(HERE / "demo_sparkle.h5ad")
print(f"Wrote {HERE / 'demo_sparkle.h5ad'}: {adata.shape}")
