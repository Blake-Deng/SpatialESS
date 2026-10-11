# SpatialESS

SpatialESS 0.1.4 implements spatial cell-cell communication inference using a
CSR radius graph and streamed ligand-receptor scoring and permutation tests.
It provides patient-aware multi-sample analysis with a random-intercept linear
mixed model. Numerical validation targets the tested SpatialCellChat V3 workflow.

## Installation

Requires R >= 4.1 and a C++17 compiler. From this repository:

```sh
Rscript -e 'install.packages(c("Rcpp", "Matrix", "lme4", "testthat"))'
R CMD INSTALL .
Rscript examples/lr_panel_example.R
```

Restart R after installation and confirm `packageVersion("SpatialESS")` is
`0.1.4`. The example uses synthetic data and does not require CellChat.

## Single-sample analysis

`expression` is a non-negative gene-by-cell sparse matrix with unique gene
and cell names. Coordinates and cell-type labels must follow its cell order.
Coordinate units must agree with the spatial parameters.

```r
library(SpatialESS)

# db is one compatible CellChatDB object, used for all three database tables.
coverage <- filter_cellchat_lr(
  lr = db$interaction,
  genes = rownames(expression),
  complex = db$complex,
  cofactor = db$cofactor
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
  tol = 5,
  interaction.range = 30,
  contact.range = 10,
  nboot = 100,
  seed.use = 1
)
head(result$records)
```

All mandatory ligand/receptor subunits must be measured. Excluded interactions
produce a warning and an explicit audit; missing genes are not imputed.
Alternatively use `lr_missing = "filter"` in `spatialess()` and inspect
`result$lr_filter`. The default is `lr_missing = "error"`.
Coverage filtering does not select spatially variable or significant LR.
See [input requirements](docs/INPUTS.md).

## Multi-patient analysis

The primary model retains individual slices and accounts for shared patient
identity:

```text
log1p(score) ~ condition + (1 | patient)
```

```r
prepared <- prepare_multisample_communication(
  results = sample_results,
  sample_metadata = sample_metadata,
  unit = "sample"
)
fit <- fit_multisample_communication_lmm(
  prepared = prepared,
  design = ~ condition,
  coefficient = "conditiondisease",
  patient_col = "patient_id",
  min_units = 6,
  min_nonzero = 2,
  REML = FALSE
)
```

The coefficient must match the actual levels in your metadata. Estimates are
reported on the original log1p-score scale. BH correction applies only to fits
with `status == "ok"`; convergence-warning fits are excluded. A singular fit
is a boundary estimate of patient variance and is recorded explicitly.
See [model interpretation](docs/MULTISAMPLE_LMM.md).

## Results

The [result index](docs/RESULTS.md) links to the current inference outputs,
resource measurements, official-comparator metrics, spatial graph summaries
and cohort model tables. Every table records its measurement scope and source.
Runtime and memory are measured observations; dense matrix size is theoretical
storage and is not peak RSS.

```text
R/, src/, man/       package interfaces, C++ kernels and reference documentation
tests/               unit tests and official-reference fixtures
examples/            runnable synthetic example
benchmarks/          configurable runners
results/             current tables, compact records and parameter contracts
figures/             figures generated from the included tables
docs/                inputs, model interpretation and result provenance
```

To reproduce a prepared sample:

```sh
Rscript benchmarks/run_prepared.R fixture.rds output_directory
```

The fixture is an R list of named arguments for `spatialess()`. Large raw data
are not included. [Reproduction notes](docs/REPRODUCE.md) describe the input
manifest and checksum verification.

## SPARKLE

The separate SpatialESS-SPARKLE repository bundles this same engine and a
sparse H5AD bridge. SPARKLE is independent upstream correction software and
must be cited separately when used.

## Scope and checks

Local installation, examples and R CMD check are recorded in
`results/checks.tsv`. GitHub cloud Actions are **NOT RUN** for this source ZIP.
Numerical agreement is not biological ground truth, and the tested execution
strategy does not establish compatibility with every communication model.

License: GPL-3.
