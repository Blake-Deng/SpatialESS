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

## Official comparator

The tested SpatialCellChat V3 reference is package version 0.1.0, commit
`be88f300464cf970b36a6165ea7d2e41a47a7511`. V3 names the workflow generation.
The target uses distance constraints, mean group aggregation
(`avg.type="avg"`), `colocalization.use=FALSE` and global label permutations.
Original and dimensions-only reference builds have distinct method labels.
Failed official jobs have no inferred numerical-fidelity result.

Percentage filtering uses the pinned V3 group-table formatting
`format(..., digits=1)` before threshold comparison, in both observed and
permuted scores. Record R version, formatting options and locale. This is not
equivalent to applying `round(x, 1)` to unformatted fractions.

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
