# Spatially stratified analytic and hybrid validation

q^Date: 2026-07-31

> Historical prototype record. The current one-way conservative hybrid and its
> zero-inflation guard are validated in
> `H1_CONSERVATIVE_HYBRID_VALIDATION_20260731.md`.
^
## Scope

The first homogeneous analytic prototype was extended to the actual SpatialESS
null structure: complete cell profiles are permuted only within fixed
sample-compartment-spatial blocks. The new experimental entry point is
`experimental_stratified_analytic_cellchat_group_support()`.

For each target group and gene, the null center is computed from a group-weighted
empirical distribution assembled from the block-specific group capacities.
The variance of a linearized triMean is the sum of finite-population block
terms. Covariance between different sender and receiver groups uses the
negative covariance of disjoint allocations within each block. The influence
correlation guard is computed after block-centering, so shared block baselines
are not incorrectly treated as cell-level ligand-receptor dependence.

The experimental hybrid entry point is
`experimental_hybrid_cellchat_group_support()`. It retains analytic p-values
only when the record is eligible and farther than `borderline_width` from
`alpha`. All other active records use the block-preserving permutation result.
The current implementation evaluates the complete supported score matrix when
any fallback is needed, so it is a correctness dispatcher, not yet an
adaptive-runtime optimization.

## Validation design

The calibration uses a 20 x 20 spatial grid, 2 alternating groups, 16 blocks
of 5 x 5 cells and a radius graph with 2,964 directed edges. Each block has an
independent spatial baseline shared by the ligand and receptor, while cell
residuals are independent under the null. The alternative multiplies ligand
expression in group A and receptor expression in group B by 1.35.

Both methods use the same supported group-pair domain, plus-one correction and
`greater_equal` upper tail.

## Corrected results

| Reference | Analytic type-I | Permutation type-I | Threshold agreement | Null p Spearman | Analytic / permutation power | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| 200 replicates x 999 permutations | 0.070 | 0.060 | 99.0% | 0.9884 | 1.00 / 1.00 | 237,952 kB |
| 20 replicates x 9,999 permutations | 0.050 | 0.050 | 100% | 0.9744 | 1.00 / 1.00 | 231,824 kB |

The 200-replicate analytic type-I 95% binomial interval was 0.0388-0.1147;
the permutation interval was 0.0314-0.1025. Analytic type-I was not
significantly distinguishable from the permutation reference at this sample
size, but the observed difference is retained rather than rounded away. The
high-permutation run is a precision reference with only 20 replicates and is
not by itself a definitive type-I study.

The first weighted-CDF interpolation implementation was deliberately retained
as an audit record because it produced an anti-conservative 0.085 analytic
type-I rate versus 0.060 for permutation. It is stored under
`benchmarks/results/h1_stratified_analytic_rep200_perm999_v1_cdf_failed/` and
is not used by the current package. The corrected implementation uses an
equal-weight type-7 path and weighted plotting positions.

## Decision

The stratified analytic candidate is suitable for continued validation and for
guarded experimental hybrid dispatch. It is not yet a default inference mode.
Remaining gates are a larger calibration matrix over group size, zero fraction,
ties and ligand-receptor dependence; complex/cofactor fallback checks; FDR
control across many LR records; and an adaptive permutation implementation that
can avoid evaluating the full matrix when fallback is sparse.
