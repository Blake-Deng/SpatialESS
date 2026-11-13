# Reference-aligned permutation inference (development 0.1.4.9003)

The main `spatialess()` entry point defaults to `aggregation="reference"`.
It forms one sparse cell-probability matrix per LR and reproduces the reference
`Matrix::crossprod(X, probability %*% X)` numerator and count operations, with
`X` the dense cell-group one-hot matrix. Observed and permuted scores share
that path. It does not build a dense cell-by-cell matrix or retain all LR cell
matrices, but it needs dense cell-by-group and group-by-group workspaces.

`max_reference_entries=1e8` guards each such matrix before allocation. Multiple
matrices can coexist: this limit is not a total-RAM bound. Exceeding it fails
explicitly. The old `aggregation="stream"` remains available for workloads
where the previous memory profile is required, with no claim of bitwise equality
at strict permutation boundaries. Historical streaming performance measurements
must not be relabeled as measurements of the new reference mode.

The RNG sequence reproduces the tested SpatialCellChat V3 future RNG advance;
future.apply is a hard dependency rather than an optional discovery branch.
Fix the R version, RNGkind (including sample.kind), dependency versions, input
cell order and group levels for paired comparisons.

Real gene names bypass the geometric mean; named complexes retain the reference
mean(log(expression)) formula, including one-subunit complexes. Sparse triMean
uses the reference type-7 interpolation and R mean steps. Cofactor products and
arithmetic association match the corresponding reference formulas.

Euclidean distances preserve separate double-precision squared differences
without compiler multiply-add contraction. Distances are calculated in raw
coordinate units and then multiplied by `ratio`, as in the reference.

The empirical convention is unchanged: strict Pboot>Pobs, count/nboot.
No epsilon or greater-equal substitution is used. Main sparse records still
omit zero/filtered entries; dense reference comparisons must match the domain.
Experimental spatial-block, adaptive and hybrid nulls remain separate methods.

This pre-release was validated against original installed SpatialCellChat 0.1.0
(source commit be88f300464cf970b36a6165ea7d2e41a47a7511) in a matched server
environment. The harness captured actual permutations/labels and high-precision
observed/null values, and compared records by LR/sender/receiver keys.
See [the validation summary](validation/pvalue_reference_20261113/README.md)
for exact scope, historical-reference limitations and the preserved statistical rule.
