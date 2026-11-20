# Inputs and statistical contract

Use a non-negative sparse gene-by-cell expression matrix with unique gene and
cell identifiers. Align coordinates, group labels and factor levels to the
matrix columns by barcode. Freeze normalization, gene rows, coordinates,
groups, database tables, selected LR, parameters and RNG settings for a paired
comparison. The supplied expression maximum enters scoring; different gene
panels can change scores even when the LR set is the same.

## Ligand-receptor coverage

`filter_cellchat_lr()` requires every measured single gene or mandatory complex
subunit. It preserves input order and reports original LR rows, missing genes
and excluded interactions. Measured all-zero rows are valid genes. Cofactors
use their documented available-gene treatment and are audited separately.
Never apply an LR mask to a different database object.

The strict default rejects unresolved LR. With `lr_missing = "filter"`, any
exclusion produces a warning and remains in `result$lr_filter`. To filter
explicitly, pass `coverage$lr` to inference.

## Group expression filtering

`percent_filter = "numeric"` (default since 0.1.4.9004) uses the fraction of
cells with positive ligand expression in the sender group and positive receptor
expression in the receiver group. Each denominator is the total number of cells
in its group. Both fractions must be `>= min.percent` (default `0.1`). For example,
8 of 82 cells falls below 10%, while 8 of 80 meets it. SpatialESS recomputes
both fractions using the group labels in each permutation. `min.percent = 0`
disables this filter. `min.cells.sr` still applies.

`percent_filter = "spatialcellchat_v3"` restores the previous comparison. V3
formats the ligand and receptor fraction columns with `format(..., digits = 1)`
and then compares the resulting strings to the threshold. The outcome can
depend on other fractions in the same column and on R's formatting settings.
For example, a fraction of `0.00005` can be printed as `"5e-05"` and pass the
character comparison against `"0.1"`. Record the R version, locale and
formatting options when reproducing that behavior.

The separate `aggregation` argument controls group score arithmetic. The default
`"reference"` uses Matrix/BLAS and checks `max_reference_entries` before
allocating dense cell-by-group or group-by-group workspaces. `"stream"` uses
native sums and can give different tail counts at floating-point ties. Both
aggregation methods accept either filter. For a saved run, inspect
`result$parameters[c("percent_filter", "min.percent", "aggregation")]`.
Older results can identify the legacy rule by
`parameters$percentage_filter == "SpatialCellChat_V3_format_digits1"`.

## Official comparator

The tested SpatialCellChat V3 reference is package version 0.1.0, commit
`be88f300464cf970b36a6165ea7d2e41a47a7511`. V3 names the workflow generation.
The target uses distance constraints, mean group aggregation
(`avg.type="avg"`), `colocalization.use=FALSE` and global label permutations.
Original and dimensions-only reference builds have distinct method labels.
Failed official jobs have no inferred numerical-fidelity result.

For paired V3 comparisons with the current package, set
`percent_filter = "spatialcellchat_v3"` and `aggregation = "reference"`.
Match R/RNG settings and numerical libraries as well as inputs and parameters.
The [0.1.4.9003 comparison](validation/pvalue_reference_20261113/README.md)
records the environment and tested inputs.

ALRA and MERINGUE are upstream preparation steps and are not run internally
by `spatialess()`. A fixed-LR comparison passes the same selected table through
official `LR.use` and SpatialESS `lr`. A database row count does not establish
which LR are actually used after official variable-feature selection.
Official group-level probability and permutation p-values require
`computeAvgCommunProb()` after `computeCommunProb()`.

## Sampling

For scaling experiments, define groups on a fixed parent dataset and extract
each subset's labels by barcode. Subsetting does not require reclustering.
The GSE306130 1k pilot labels derive from PCA/k-means on the 5k parent inputs;
the sample-specific preparation seeds are included in the sampling table.

Central-region sampling retains a local tissue neighborhood. Global
proportional sampling preserves group composition but can disperse sampled
cells and reduce spatial degree. Both are represented in the cervical tests;
they answer different geometry questions. The exact prepared input is
identified by its checksum, cell count, gene count and parameter contract.

## Output interpretation

`records` contains active sender-receiver-LR scores and raw permutation
p-values. Counts at p<0.05 and p<=0.05 are separate fields, not BH FDR calls.
The cohort LMM reports its own BH-adjusted q-values with fit-status checks.
Empty but completed results and missing inference stages are distinct states.
