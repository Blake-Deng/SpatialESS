# SpatialESS

SpatialESS is an independent, high-performance implementation of the tested
SpatialCellChat V3 spatial communication workflow. It keeps the published
probability, distance, Hill-function, triMean and permutation definitions,
while replacing dense intermediate objects with sparse, streamed operations.
The package also provides downstream multi-patient communication models, including a patient-random-intercept linear mixed model.

The validated single-sample inference entry point is:

```r
SpatialESS::spatialess()
```

SpatialESS is the validated single-sample inference engine.
The package exposes one validated single-sample engine and a downstream multi-patient analysis layer.

Agreement with official SpatialCellChat is reported as numerical fidelity, not
as biological ground truth. Biological support is assessed separately at the
patient and pathway levels.

## Install

```r
install.packages("remotes")
remotes::install_github("Blake-Deng/SpatialESS")
```

## Minimal workflow

```r
library(Matrix)
library(SpatialESS)

# Load the database from the package used in the V3 comparison.
# If a database is already loaded, use that SAME object for all three tables.
db_env <- new.env()
data("CellChatDB.human", package = "SpatialCellChat", envir = db_env)
db <- db_env$CellChatDB.human

# expression: normalized signaling-gene x cell dgCMatrix, with gene symbols.
# Keep the complete signaling matrix; do not subset it to selected LR genes.
coverage <- filter_cellchat_lr(
  lr = db$interaction, genes = rownames(expression),
  complex = db$complex, cofactor = db$cofactor
)
write.csv(coverage$audit, "lr_coverage.csv", row.names = FALSE)
write.csv(coverage$missing_genes, "lr_missing_genes.csv", row.names = FALSE)

result <- spatialess(
  expression = expression,
  coordinates = coordinates,
  group = cell_type,
  lr = coverage$lr,
  complex = db$complex,
  cofactor = db$cofactor,
  tol = 5, interaction.range = 30,
  contact.range = 10,
  nboot = 100,
  seed.use = 1
)

head(result$records)
```

Targeted Xenium/CosMx panels do not measure every gene in CellChatDB. Missing
ligand/receptor genes or any required complex subunit must be reported before
inference. The filter preserves LR order and checks gene presence; it does not
test overexpression, spatial variability, or significance. Cofactors retain
the existing available-gene handling and are reported separately when missing.
Zero-expression gene rows remain measured genes. Do not fill unmeasured genes
with invented expression values. A runnable, database-free example is in
[`examples/lr_panel_example.R`](examples/lr_panel_example.R).

Alternatively, pass the full interaction table and explicitly set
`lr_missing = "filter"` in `spatialess()`. Exclusions produce a warning; the
report is saved in `result$lr_filter$audit`, `$excluded`, and `$missing_genes`.
The default `lr_missing = "error"` retains strict behavior, with an actionable
diagnostic. Calling the filter alone does not modify the original database:
you must pass `coverage$lr` into inference.

For official numerical fidelity comparisons, use the same prepared signaling
expression matrix, selected LR, complex/cofactor tables, coordinates, group
order, seed and parameters in both engines. Official V3 can additionally
select spatially variable LR; its selected `chat@LR$LRsig` can be used directly.
Coverage-filtered candidates are not equivalent to this additional selection.
V3 refers to the workflow generation; the benchmark SpatialCellChat package
reports version 0.1.0. This does not certify every CellChat 2.x database version.

Version 0.1.3 retains the panel-coverage interface and fixes the official V3
percentage-filter boundary in both observed scores and permutations. It uses
the official `format(..., digits = 1)` comparison on the small group-level
fraction table, rather than comparing unrounded fractions. This can intentionally
change results from 0.1.2 when `min.percent > 0`. See
[`RELEASE_STATUS.md`](RELEASE_STATUS.md) for local checks and regression evidence.

`spatialess()` does not run ALRA or MERINGUE internally. Official
`preProcessing()` performs ALRA and overwrites `chat@data.signaling`; if used,
pass that same resulting matrix and selected `chat@LR$LRsig` to SpatialESS.
Official `computeCommunProb()` produces cell-level results; run
`computeAvgCommunProb()` before comparing group-level probabilities and p-values.
Missing output stages must not be counted as valid zero results. See
[`docs/V3_COMPATIBILITY.md`](docs/V3_COMPATIBILITY.md).

## Main validation evidence

The compact result tables are in results/tables/. The full interpretation and
provenance are in docs/MAIN_VALIDATION.md.

### Numerical fidelity

- GSE250346: 45 samples from 35 patients. SpatialESS completed all 45. The
  unmodified official workflow completed 41/45; TILD117LA completed after an
  official sparse-matrix dimension repair; three larger samples reached a
  Matrix::crossprod() signal 11 at FN1-CD44 because the complete sparse outer
  product exceeded the 32-bit dgCMatrix index limit.
- GSE306130: four 1,000-cell samples, 33 LR pairs and 20 permutations. This is
  an independent four-sample fidelity pilot, not the main patient-level
  disease analysis.
- GSE313006: an independent 1,000-cell fixture with 133 LR pairs, tested with
  20 and 100 permutations.
- CosMx Cancerous Liver: paired official/SpatialESS comparisons from 250
  through 50,000 cells, with SpatialESS-only continuation to larger scales.

Across the historical completed paired comparisons at their recorded settings,
active-record Jaccard was 1, probability
differences were at floating-point scale, p-values agreed exactly where
reported, and no significance threshold calls differed. These historical
results do not certify every parameter combination: the subsequent cervical
external test exposed a default percentage-filter mismatch in 0.1.2. Its
failure evidence is retained separately from the 0.1.3 acceptance tests.

### Multi-patient linear mixed model

The primary biological dataset is GSE250346: 45 spatial slices from 35 patients (9 controls and 26 pulmonary fibrosis patients). SpatialESS first infers communication independently within each slice. The primary cohort model then retains all slices and accounts for repeated slices from the same patient:

```text
log1p(score) ~ condition + (1 | patient)
```

```r
prepared <- prepare_multisample_communication(
  results = sample_results,
  sample_metadata = sample_metadata,
  unit = "sample"
)

lmm_result <- fit_multisample_communication_lmm(
  prepared = prepared,
  design = ~ condition,
  coefficient = "conditiondisease",
  patient_col = "patient_id",
  min_units = 6,
  min_nonzero = 2,
  REML = FALSE
)
```

`fit_multisample_communication_lmm()` uses `lme4::lmer()`. It checks fixed-effect design rank, records zero-score fraction and response variation, standardizes each feature response for numerical conditioning, and returns estimates on the original log1p-score scale. It tries `bobyqa`, `nloptwrap` and `Nelder_Mead`, preferring a warning-free fit, and records optimizer, random-effect variance, residual variance, maximum gradient and diagnostic messages. Wald normal-approximation p-values are used and BH correction is applied only to fits with `status = "ok"`.

The simple slice-level model is retained as a comparator:

```r
slice_glm <- fit_multisample_communication_glm(
  prepared = prepared,
  design = ~ condition,
  coefficient = "conditiondisease"
)
```

The returned LMM status fields are interpreted as follows:

- `ok`, `singular = FALSE`: converged mixed model with nonzero patient variance;
- `ok`, `singular = TRUE`: converged boundary model with patient variance at or near zero;
- `convergence_warning`: coefficients are retained for diagnosis but excluded from FDR;
- `insufficient_zero_support`: zero-score fraction exceeded the configured support threshold;
- `insufficient_variation`: response variation was below the configured threshold;
- `numerical_failure`: all configured optimizers failed with a recognized numerical error;
- `fit_failed`: all configured optimizers failed with an unrecognized error;
- `insufficient_nonzero`: too few slices express the communication feature.

A singular fit is a boundary diagnostic, not a software failure. The completed 2026-09-07 robust GSE250346 rerun includes 45 slices from 35 patients (9 Control, 26 PF): 39,570 valid LMM fits (22,262 non-singular, 17,308 singular), 1,984 convergence warnings, 292 numerical failures and 8,194 features below the nonzero-slice threshold. The LMM identified 5,147 BH-significant features versus 4,172 for the slice-level GLM; 4,108 overlapped (Jaccard 0.7883), with effect-estimate Spearman rho 0.9966 across valid comparisons. These are Wald-normal screening results, not experimentally established disease mechanisms. All 39,242 previously valid effect estimates are unchanged. Full tables, audits, source and input checksums are indexed in [`docs/LMM_RESULTS_20260907.md`](docs/LMM_RESULTS_20260907.md); the old summaries remain in `results/tables/gse250346_lmm_baseline/`. See [`docs/REVIEWER_GUIDE.md`](docs/REVIEWER_GUIDE.md) for model diagnostics and inference limitations.

The existing patient-aggregated GLM, leave-one-patient-out analysis, aggregation choices and minimum-cell thresholds remain sensitivity analyses. These results assess cohort-level robustness; they do not turn inferred communication into experimental truth.

### Scale and resource evidence

The CosMx ladder provides paired runtime and peak-RSS measurements. The
million-cell Xenium run is retained as a stress test for the current SpatialESS
engine. Large source objects and raw matrices are deliberately excluded from
Git; their paths, checksums and regeneration commands belong in a local
provenance manifest.

## Validation summary

The primary fidelity benchmark used nested cells from the CosMx CancerousLiver
sample, all 373 computable LR interactions, 20 permutations, the same seed,
and one thread. From 250 through 50,000 cells:

- active-record Jaccard was 1;
- probability Pearson correlation was 1;
- maximum absolute probability difference was at most `7.81e-18`;
- p-values were identical; and
- no result crossed the 0.05 threshold differently.

At 50,000 cells, SpatialESS completed in 10.48 s with 0.802 GiB peak RSS.
SpatialCellChat V3 required 23,919 s and 60.78 GiB. The official comparator did
not complete at 100,000 cells or above; SpatialESS completed 443,515 cells in
101.75 s with 1.97 GiB peak RSS.

Large Xenium and CosMx runs are treated as scalability stress tests. Their
permutation counts, minimum attainable p-values and resource limits are
reported in the accompanying provenance tables rather than being presented as
independent biological truth.

See [the main validation report](docs/MAIN_VALIDATION.md) and
[the reproduction notes](docs/REPRODUCE_MAIN.md)
for exact parameters, interpretation boundaries, and evidence paths.

## Repository layout

R/                  R interfaces, SpatialESS engine and multi-patient LMM
src/                C++17 kernels
tests/              unit and numerical regression tests
benchmarks/         reproducible benchmark drivers
results/tables/     small, reviewable summary tables
docs/               methods, validation and reproduction notes
figures/            lightweight explanatory figures

## Reproducing the paper summaries

The benchmark drivers use prepared inputs and are not intended to download
large public datasets automatically. Set the relevant source roots and run:

    Rscript benchmarks/summarize_main_results.R

See docs/REPRODUCE_MAIN.md for expected source directories, parameter locks,
and the distinction between paired and SpatialESS-only runs.

## Scope

SpatialESS is a computational implementation and optimization of a
SpatialCellChat-compatible model, not a claim that one CCC model is biological
ground truth. The release supports radius-graph single-sample analysis.
Segmentation-polygon leakage correction, morphology-aware models and external
experimental validation are separate extensions and are not required by the
main SpatialESS benchmark.
