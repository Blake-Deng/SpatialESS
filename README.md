# SpatialESS

SpatialESS is an independent, high-performance implementation for scalable
spatial cell-cell communication analysis. Its main paper evidence has three
parts:

Release snapshot: see docs/RELEASE_0.1.0.md for frozen parameters, environment, and final benchmark status.

1. Numerical fidelity: spatial_v3_compatible() preserves the documented
   SpatialCellChat V3 probability and permutation definitions while replacing
   full cell-by-cell intermediate matrices with streamed sparse operations.
2. Scalability: spatialess() uses compact CSR radius graphs, expression
   support pruning, exact sparse type-7 summaries, and streamed permutation
   inference to reduce runtime and peak memory.
3. Patient-level biology: sample-level communication records are aggregated
   with the patient as the statistical unit and analyzed with FDR-controlled
   GLMs for conserved and disease-associated programs.

The package does not claim that official agreement is biological ground truth.
Official agreement is reported as numerical fidelity; biological support is
assessed separately at the patient and pathway levels.

The package has two deliberately separate entry points:

- `spatial_v3_compatible()` preserves the SpatialCellChat V3 probability and
  permutation definitions and is the numerical-fidelity benchmark path.
- `spatialess()` is the integrated single-sample workflow with separate
  contact/diffusion graphs and spatial block-permutation inference.

The integrated workflow reports only fixed permutation inference or a guarded
hybrid mode in which every significant result is confirmed by fixed
permutation. Direct analytic and adaptive p-values remain experimental.

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
  contact_radius = 15,
  diffusion_radius = 35,
  block_size = 100,
  inference = "permutation",
  nperm = 100,
  seed = 1
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

### Patient-level GLM

The primary biological dataset is GSE250346. The prepared metadata contains 45
sample records from 35 patients, with 9 controls and 26 disease patients.
Communication is first summarized within sample and then analyzed with the
patient as the unit of inference. The current outputs include:

- disease-associated LR, sender-receiver and pathway summaries;
- FDR-adjusted conserved and disease-associated programs;
- leave-one-patient-out analysis;
- mean, median and cell-weighted aggregation sensitivity;
- minimum-cell and mixed-model sensitivity analyses.

These results support robustness of the patient-level communication programs;
they do not turn an inferred communication score into experimental truth.

### Scale and resource evidence

The CosMx ladder provides paired runtime and peak-RSS measurements. The
million-cell Xenium run is retained as a stress test for the current integrated
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

The integrated workflow completed:

- 443,515 CosMx cells, 373 LR and 100 spatial permutations in 46.60 s with
  1.95 GiB peak RSS;
- 1,156,091 Xenium cells, 627 LR, 2,357 groups and 20 permutations in 79.78 s
  with 5.61 GiB peak RSS.

The Xenium run is a scalability stress test: with 20 permutations the minimum
finite-corrected p-value is `1/21`, so it is not presented as a definitive
high-resolution biological analysis.

See [the Stage 2 report](docs/STAGE2_INTEGRATED_RELEASE_20260804.md) and
[the publication validation protocol](docs/PUBLICATION_VALIDATION_20260804.md)
for exact parameters, interpretation boundaries, and evidence paths.

## Repository layout

R/                  R interfaces and workflows
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
