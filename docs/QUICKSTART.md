# Analyze your own data

Install SpatialESS and run the [synthetic example](../examples/lr_panel_example.R)
first. The examples below use your own expression data, spatial coordinates,
cell-type labels and a compatible CellChatDB object named `db`.

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
See [input requirements](INPUTS.md).

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
See [model interpretation](MULTISAMPLE_LMM.md).


The [reproduction guide](REPRODUCE.md) covers saved manuscript data and figures.
