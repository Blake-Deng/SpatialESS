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

result <- spatialess(
  expression = expression,
  coordinates = coordinates,
  group = cell_type,
  lr = CellChatDB.human$interaction,
  complex = CellChatDB.human$complex,
  cofactor = CellChatDB.human$cofactor,
  tol = 5, interaction.range = 30,
  contact.range = 10,
  nboot = 100,
  seed.use = 1
)

head(result$records)
```

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

Across completed paired comparisons, active-record Jaccard was 1, probability
differences were at floating-point scale, p-values agreed exactly where
reported, and no significance threshold calls differed.

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
- `ok`, `singular = TRUE`: converged boundary model with patient variance estimated as zero;
- `convergence_warning`: coefficients are retained for diagnosis but excluded from FDR;
- `insufficient_zero_support`: zero-score fraction exceeded the configured support threshold;
- `insufficient_variation`: response variation was below the configured threshold;
- `numerical_failure`: all configured optimizers failed with a recognized numerical error;
- `fit_failed`: both optimizers failed and no inferential result is reported;
- `insufficient_nonzero`: too few slices express the communication feature.

A singular fit is not a software failure. It means that the data do not support a nonzero patient random-intercept variance for that feature. The archived baseline GSE250346 run had 39,242 features with status `ok` (21,934 non-singular and 17,308 singular), 2,299 convergence warnings, 305 numerical failures and 8,194 features below the nonzero-slice threshold. The robust implementation adds support screens and a third optimizer; its final cohort counts must be regenerated with the release script before being used as confirmatory manuscript numbers. The archived LMM identified 5,158 FDR-significant features versus 4,172 for the slice-level GLM; 4,109 overlapped (Jaccard 0.787), and effect estimates had Spearman rho 0.9967 across valid fits. Full status counts and the GLM-versus-LMM comparison are stored in `results/tables/`. A reviewer-oriented explanation of model diagnostics versus software failures is in [`docs/REVIEWER_GUIDE.md`](docs/REVIEWER_GUIDE.md).

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
