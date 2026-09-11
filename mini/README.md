# SpatialESS mini

This small repository is a reviewer-facing reproduction of the main
SpatialESS engine. It uses an eight-cell synthetic expression matrix, a
complex/cofactor LR example, a CSR radius graph and the same streamed
probability code used by the full package. It does not require any large
public dataset.

## Run

From this mini directory, install the parent package first:

```bash
R CMD INSTALL ..
./run.sh
```

The script requires the official SpatialCellChat package. It checks sparse
triMean and group molecular scoring against its corresponding helper functions
(not a fresh full official SpatialCellChat inference run).
It exits non-zero if either maximum difference exceeds
1e-14 and writes `spatialess_mini_records.tsv` on success.

## What is demonstrated

- sparse `dgCMatrix` input is accepted without dense conversion;
- only LR-referenced genes are scored;
- spatial support is represented by a CSR radius graph;
- complex, cofactor, agonist and antagonist terms use the CellChat formulas;
- the resulting group-pair probabilities match the direct reference at
  floating-point tolerance.

The full patient-level workflow is available in the main package:

```r
prepared <- prepare_multisample_communication(
  sample_results, sample_metadata, unit = "sample"
)
lmm_result <- fit_multisample_communication_lmm(
  prepared, design = ~ condition,
  coefficient = "conditiondisease", patient_col = "patient_id"
)
```

For the real-data LR interface regression, run from the repository root:

```bash
Rscript benchmarks/verify_lr_interface_mini.R .
```

This uses the bundled cervical 1k fixture and full-precision expected RDS
from dimension-patched official V3 (percentage filtering unmodified), not
an ESS-generated expected result. It covers the 8/82 percentage boundary.
This real-data mini does not require the SpatialCellChat package: its database
is embedded in the fixture.
