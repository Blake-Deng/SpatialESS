import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[2]
BRIDGE = ROOT / "inst" / "python" / "h5ad_to_spatialess.py"


class H5ADBridgeTest(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        expression = sparse.csr_matrix(
            np.array(
                [
                    [10.0, 0.0, 0.0, 90.0],
                    [0.0, 20.0, 80.0, 0.0],
                    [5.0, 5.0, 0.0, 0.0],
                ]
            )
        )
        obs = pd.DataFrame(
            {
                "annotation": ["A", "B", "A"],
                "x_um": [0.0, 10.0, 20.0],
                "y_um": [1.0, 11.0, 21.0],
            },
            index=["c1", "c2", "c3"],
        )
        var = pd.DataFrame(index=["L", "R", "G3", "G4"])
        self.h5ad = self.root / "input.h5ad"
        ad.AnnData(X=expression, obs=obs, var=var).write_h5ad(self.h5ad)

    def tearDown(self):
        self.tempdir.cleanup()

    def run_bridge(self, output, *extra):
        command = [
            sys.executable,
            str(BRIDGE),
            str(self.h5ad),
            str(output),
            "--group-col",
            "annotation",
            "--coordinate-cols",
            "x_um,y_um",
            *extra,
        ]
        subprocess.run(command, check=True, capture_output=True, text=True)

    def test_all_gene_denominator_and_gene_selection(self):
        genes = self.root / "genes.txt"
        genes.write_text("L\nR\nMISSING\n")
        output = self.root / "bundle"
        self.run_bridge(
            output,
            "--gene-list",
            str(genes),
            "--normalization",
            "log1p_cp10k",
            "--chunk-size",
            "2",
        )

        matrix = mmread(output / "expression.mtx.gz").toarray()
        expected = np.array(
            [
                [np.log1p(1000.0), 0.0, np.log1p(5000.0)],
                [0.0, np.log1p(2000.0), np.log1p(5000.0)],
            ]
        )
        np.testing.assert_allclose(matrix, expected, rtol=0, atol=1e-12)
        self.assertEqual((output / "cells.tsv").read_text().splitlines(), ["c1", "c2", "c3"])
        self.assertEqual((output / "genes.tsv").read_text().splitlines(), ["L", "R"])

        manifest = json.loads((output / "manifest.json").read_text())
        self.assertEqual(manifest["normalization_denominator"], "all H5AD genes per cell")
        self.assertEqual(manifest["missing_requested_genes"], ["MISSING"])
        self.assertEqual(manifest["cells"], 3)
        self.assertEqual(manifest["input_genes"], 4)
        self.assertEqual(manifest["exported_genes"], 2)

    def test_external_metadata_is_realigned_to_h5ad_order(self):
        metadata = pd.DataFrame(
            {
                "cell": ["c3", "c1", "c2"],
                "annotation": ["C", "A", "B"],
                "x_um": [30.0, 10.0, 20.0],
                "y_um": [31.0, 11.0, 21.0],
            }
        )
        metadata_path = self.root / "metadata.tsv"
        metadata.to_csv(metadata_path, sep="\t", index=False)
        output = self.root / "external_bundle"
        self.run_bridge(
            output,
            "--metadata-tsv",
            str(metadata_path),
            "--cell-col",
            "cell",
            "--normalization",
            "none",
        )

        exported = pd.read_csv(output / "metadata.tsv", sep="\t")
        self.assertEqual(exported["cell"].tolist(), ["c1", "c2", "c3"])
        self.assertEqual(exported["group"].tolist(), ["A", "B", "C"])
        self.assertEqual(exported["coordinate_1"].tolist(), [10.0, 20.0, 30.0])


if __name__ == "__main__":
    unittest.main()
