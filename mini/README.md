# SpatialESS mini

This small repository is a reviewer-facing reproduction of the main
SpatialESS engine. It uses an eight-cell synthetic expression matrix, a
complex/cofactor LR example, a CSR radius graph and the same streamed
probability code used by the full package. It does not require any large
public dataset.

## Run

Install the main package from the sibling repository first:

```bash
R CMD INSTALL ../SpatialESS
./run.sh
```

The script checks the sparse triMean and LR probability against a direct
CellChat reference. It exits non-zero if either maximum difference exceeds
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
  sample_results, sample_metadata, unit = "patient"
)
glm_result <- fit_multisample_communication_glm(
  prepared, design = ~ condition, coefficient = "conditiondisease"
)
```
