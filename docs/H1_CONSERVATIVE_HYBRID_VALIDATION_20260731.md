# Conservative hybrid validation and Stage 1B decision

Date: 2026-07-31

## Scope

This validation covers the fixed-permutation conservative hybrid implemented
by `experimental_hybrid_cellchat_group_support()`. It does not promote the
adaptive Monte Carlo p-values to BH/FDR use.

The decision rule is one-way:

- eligible analytic `p > alpha + screening_margin`: assign hybrid `p = 1`;
- all other active records: run fixed block-preserving permutation;
- every final discovery must have `permutation_confirmed = TRUE`;
- complex, cofactor, regulator, autocrine, small-group, insufficient-block,
  excessive-tie, excessive-zero and strong-dependence records require
  permutation.

The formal runs use `alpha=0.05`, `screening_margin=0.01`,
`max_zero_fraction=0.10`, the plus-one correction and the
`greater_equal` empirical tail.

## Why the first expanded matrix failed

The first 108-scenario matrix varied grid/group size, zero fraction, ties and
LR dependence. It completed all scenarios and was conservative relative to
the full permutation reference, but it did not preserve every discovery:

- minimum null record-level agreement: 0.984;
- minimum alternative record-level agreement: 0.948;
- 12 scenarios had at least one null rejection-set mismatch;
- 34 scenarios had at least one alternative rejection-set mismatch.

All mismatches were false-negative screening, not analytic-only discoveries.
They occurred only in simulated zero-inflated settings. The previous guard
rejected zero quartiles but could still accept records with positive
quartiles and substantial zero mass. Those records could have a large
analytic p-value and a significant permutation p-value.

Historical result directory:

`benchmarks/results/h1_hybrid_parameter_matrix_rep25_perm499/`

## Correction

Both homogeneous and stratified analytic engines now compute zero-expression
mass for each candidate gene. The stratified implementation uses the same
group-capacity weights as its weighted null CDF. A gene/group profile with
zero mass above `max_zero_fraction=0.10` receives
`excessive_zero_fraction` and is sent to fixed permutation.

This guard changes only inference dispatch. It does not change sparse type-7
triMean, CellChat molecular scoring, the Hill transform, spatial support, the
observed communication probability or the fixed-permutation null.

Regression tests include moderate 20 percent zero mass with positive
quartiles, in addition to the existing zero-quartile and tie guards.

## Main multi-LR results

The main design contains 400 cells, two groups, 16 spatial blocks and 10 LR
records: eight simple LR records, one ligand complex and one receptor
cofactor. Two A-to-B signals are implanted under the alternative.

| Metric | 100 x 499 | 20 x 9,999 |
|---|---:|---:|
| Analytic completion | 0.80 | 0.80 |
| Analytic-only null BH FWER | 0.04 | 0.10 |
| Conservative hybrid null BH FWER | 0.03 | 0.00 |
| Full permutation null BH FWER | 0.03 | 0.00 |
| Hybrid implanted-signal power | 1.00 | 1.00 |
| Full permutation implanted-signal power | 1.00 | 1.00 |
| Hybrid mean FDP | 0.0233 | 0.0167 |
| Full permutation mean FDP | 0.0233 | 0.0167 |
| Null rejection-set agreement | 1.00 | 1.00 |
| Alternative rejection-set agreement | 1.00 | 1.00 |
| Every hybrid discovery permutation-confirmed | yes | yes |
| Wall time | 138.40 s | 523.57 s |
| Peak RSS | 0.229 GiB | 0.224 GiB |

Machine-readable results:

- `benchmarks/results/h1_multilr_conservative_zero_guard_rep100_perm499/`
- `benchmarks/results/h1_multilr_conservative_zero_guard_rep20_perm9999/`

## Expanded parameter matrix

The corrected matrix contains 108 scenarios:

- grid sides 8, 12 and 20, corresponding to 64, 144 and 400 cells;
- zero fractions 0, 0.25, 0.50 and 0.75;
- tie quantization steps 0, 0.25 and 1;
- LR shared-dependence levels 0, 0.50 and 0.85;
- 25 null/alternative replicate pairs and 499 permutations per scenario.

All 108 scenarios completed. Across 2,700 null and 2,700 alternative
replicates:

- hybrid FWER never exceeded the matched full-permutation FWER;
- every hybrid discovery was permutation-confirmed;
- every hybrid/full-permutation rejection set was identical;
- minimum null and alternative record-level agreement were both 1.00.

For zero-free scenarios, analytic completion ranged from 0 to 0.80 depending
on group size, ties and dependence; its mean was 0.351. Thus the method still
uses analytic screening inside its declared domain. Scenarios with simulated
zero fraction at least 0.25 used permutation fallback.

Machine-readable results:

`benchmarks/results/h1_hybrid_parameter_matrix_zero_guard_rep25_perm499/`

The 16-job matrix completed in 153.10 s. Driver peak RSS was 0.214 GiB;
individual worker RSS is retained in every scenario directory.

## Zero-threshold boundary scan

A separate boundary scan used zero fractions 0, 0.02, 0.05, 0.10, 0.15, 0.20
and 0.25; grid sides 12 and 20; dependence 0 and 0.85; 100 replicate pairs and
499 permutations for each of 28 scenarios.

| Simulated zero fraction | Mean analytic completion | Hybrid null FWER | Permutation null FWER | Rejection-set agreement |
|---:|---:|---:|---:|---:|
| 0.00 | 0.7958 | 0.0375 | 0.0375 | 1.00 |
| 0.02 | 0.7990 | 0.0375 | 0.0375 | 1.00 |
| 0.05 | 0.7878 | 0.0300 | 0.0300 | 1.00 |
| 0.10 | 0.0615 | 0.0375 | 0.0375 | 1.00 |
| 0.15 | 0.0000 | 0.0575 | 0.0575 | 1.00 |
| 0.20 | 0.0000 | 0.0375 | 0.0375 | 1.00 |
| 0.25 | 0.0000 | 0.0425 | 0.0425 | 1.00 |

All 28 scenarios had exact null and alternative rejection sets and all hybrid
discoveries were permutation-confirmed. The scan supports 0.10 as the upper
boundary of the currently declared analytic domain; it is a conservative
software guard, not a universal biological constant.

Machine-readable results:

`benchmarks/results/h1_hybrid_zero_boundary_rep100_perm499/`

## Statistical interpretation

For each record, the hybrid p-value equals the full permutation p-value when
permutation is run and equals 1 when analytic screening excludes the record.
It is therefore coordinatewise no smaller than the conceptual full
permutation vector. A coordinatewise monotone correction cannot produce an
extra hybrid rejection. Exact equality is an empirical power gate and was
verified over the declared matrix; it is not assumed as a theorem outside
that domain.

The absolute FWER estimates remain simulation estimates with binomial
uncertainty. The guarantee established by implementation is relative
conservativeness and permutation confirmation, not universal FDR control
under arbitrary LR dependence.

## Stage 1B decision

The guarded fixed-permutation hybrid passes the declared Stage 1B gate and
can be frozen for Stage 2 comparisons. The following remain experimental and
must not be presented as validated publication inference:

- analytic p-values used directly as final significant p-values;
- adaptive/sequential p-values used in BH/FDR;
- analytic inference outside the calibrated group-size, zero, tie and
  dependence domain.

## Comparison with original R CellChat

The original R CellChat comparison must separate two questions.

1. Observed statistic compatibility: use identical expression, groups,
   CellChatDB records, triMean settings and Hill parameters. SpatialESS
   compatibility probabilities should agree to floating-point tolerance.
2. Significance model: original CellChat globally permutes group labels,
   whereas SpatialESS preserves spatial exchangeability blocks. These
   p-values answer different null hypotheses and should not be claimed to be
   numerically identical.

The benchmark should report observed probability differences, aggregate
count/weight agreement, runtime and peak RSS for both methods. Spatial
significance should additionally be compared with a spatial method such as
SpatialCellChat V3, with the null model stated explicitly.
