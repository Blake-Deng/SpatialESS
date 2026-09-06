# Experimental simple-LR analytic calibration

Date: 2026-07-31

## Scope

This experiment evaluates the first Stage 1B analytic candidate. It is not a
general replacement for spatial permutation. The implemented domain is limited
to one homogeneous finite-population label-permutation stratum, single-gene
ligand and receptor records, paracrine group pairs, continuous positive
quartiles, sufficiently large groups and moderate ligand-receptor dependence.

The implementation is exposed as
`experimental_analytic_cellchat_group_support()`. Every unsupported record is
returned with `inference_mode = "permutation_fallback_required"` and a specific
`fallback_reason`; it never silently applies the approximation outside this
domain.

## Method

The null triMean for each gene is linearized with the joint influence function
of Q1, Q2 and Q3. The three order statistics are therefore not treated as
independent marginals. For ligand in sender group A and receptor in receiver
group B, the finite-population covariance includes:

- covariance of ligand and receptor influence values in the same cell profile;
- negative covariance induced because A and B are disjoint samples from the
  same finite cell population;
- finite-population corrections for both group triMeans.

The log ligand-receptor product is then evaluated with a first-order delta
method. Because the CellChat Hill function is monotone, its upper-tail decision
is the same as that of the underlying positive ligand-receptor product.

Automatic fallback currently covers spatially stratified blocks, autocrine
pairs, groups smaller than 30 cells, complexes, cofactors, agonists,
antagonists, zero quartiles, excessive ties, unstable quantile-density
estimates, strong influence correlation and oversized influence workspaces.

## Reproducible commands

```bash
SPATIALESS_ANALYTIC_REPLICATES=200 \
SPATIALESS_ANALYTIC_NPERM=999 \
Rscript benchmarks/run_h1_analytic_simple_calibration.R

SPATIALESS_ANALYTIC_REPLICATES=20 \
SPATIALESS_ANALYTIC_NPERM=9999 \
Rscript benchmarks/run_h1_analytic_simple_calibration.R
```

Each replicate contains a homogeneous 400-cell, two-group continuous-expression
null and an independently simulated implanted A-to-B signal. Empirical
p-values use `tail = "greater_equal"` and the plus-one finite correction.

## Results

| Reference | Analytic type-I | Permutation type-I | Null threshold agreement | Null p Spearman | Analytic / permutation power | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| 200 replicates x 999 permutations | 0.040 | 0.035 | 99.5% | 0.9905 | 1.00 / 1.00 | 238,012 kB |
| 20 replicates x 9,999 permutations | 0.050 | 0.050 | 100% | 0.9910 | 1.00 / 1.00 | 237,976 kB |

For the 200-replicate run, the analytic type-I 95% binomial interval was
0.0174-0.0773 and the permutation interval was 0.0142-0.0708. Mean absolute
null p-value difference was 0.0287. For the high-permutation run it was 0.0279.
The larger individual p-value differences did not materially affect the 0.05
decision in these experiments, but they preclude claiming numerical identity.

Machine-readable records are in:

- `benchmarks/results/h1_analytic_simple_rep200_perm999/`
- `benchmarks/results/h1_analytic_simple_rep20_perm9999/`

## Gate decision

The candidate passes the initial homogeneous continuous simple-LR calibration
gate. It does not yet pass the full Stage 1B gate. Required next work includes
calibration across group size, zero fraction, ties and dependence; a joint
stratified influence model for multiple spatial blocks; complex/cofactor
fallback validation; and adaptive permutation around borderline decisions.
FFT convolution remains a possible acceleration for a validated distributional
representation, not a substitute for retaining joint triMean dependence.
