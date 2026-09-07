# Reviewer guide: patient-aware LMM results

## Analysis unit and inference

The GSE250346 analysis contains 45 spatial slices from 35 patients (9 Control,
26 PF), not 45 independent patients. Each sender-receiver-LR feature uses:

```text
LMM: log1p(score) ~ condition + (1 | patient)
GLM comparator: log1p(score) ~ condition
```

Multiple slices remain as rows; the LMM includes a patient random intercept.
Cells are not independent patient replicates. The models use the same score
inputs and a PF-versus-Control coefficient. The main model does not adjust
for TMA, run or affected region; retaining these metadata does not mean they
were fitted as covariates.

P-values use a Wald normal approximation, followed by BH among valid fits.
These are screening associations. Model convergence and numerical checks do
not establish finite-sample calibration, remove possible condition/batch
confounding, or validate a biological mechanism. Small-sample inference may
require Satterthwaite/Kenward-Roger methods or a suitable bootstrap.

## Completed robust rerun

All 50,040 features were rerun using the distributed robust source. The
table below partitions the entire feature set; categories do not overlap.

| Status | Features | Treatment |
|---|---:|---|
| `ok`, non-singular | 22,262 | Included in primary BH |
| `ok`, singular | 17,308 | Included in primary BH; boundary flag retained |
| `convergence_warning` | 1,984 | Diagnostic estimates only; q-value is NA |
| `numerical_failure` | 292 | No p-value or q-value |
| `insufficient_nonzero` | 8,194 | Not fitted; fewer than two nonzero slices |

There are no additional `insufficient_variation`, `insufficient_zero_support`,
`rank_deficient` or unclassified `fit_failed` records in this rerun. Support
screens were present, but recovery in this cohort must not be attributed to
their excluding extra features.

## Change from the original run

Original statuses are preserved under `results/tables/gse250346_lmm_baseline/`.
In that run, 305 features had the actual status `fit_failed`, not a retroactive
`numerical_failure` label. With the robust implementation:

- 328 previous convergence warnings became valid `ok` fits.
- 13 previous failures returned estimates with warnings, not valid fits.
- 292 previous failures remain numerical failures after all three optimizers.
- All 39,242 previously valid effect estimates are exactly unchanged.

Thus valid fits increased to 39,570. BH-significant features changed from
5,158 to 5,147; 5,140 are shared. BH uses the new valid-fit universe, so
unchanged coefficients need not imply unchanged q-values or discovery counts.
No previous `fit_failed` feature became an `ok` fit.

## Singular fits and convergence warnings

`isSingular(tol = 1e-4)` flags a patient variance component at or near the
boundary. This is not a software crash, proof that its variance is exactly
zero, or proof that patients are interchangeable. Among valid nonsingular
fits, 3,266 are significant using primary global BH; recomputing BH within
only the 22,262 nonsingular fits gives 3,476. These are distinct sensitivity
analyses with different testing universes.

Warnings are classified once per feature, using the precedence negative
Hessian eigenvalue, numerically singular Hessian, combined gradient/scale,
gradient alone, then scale alone:

| Diagnostic class | Features |
|---|---:|
| Negative Hessian eigenvalue | 468 |
| Numerically singular Hessian | 59 |
| Gradient and scale warning | 1,227 |
| Gradient warning | 19 |
| Scale warning | 211 |

This classification differs from the baseline summary, so class-by-class
counts are not directly comparable. Full messages and optimizer attempts are
retained in the quoted feature TSV, not discarded or called successes.

## Residual numerical failures

All 292 residual failures returned the following factorization error for
`bobyqa`, `nloptwrap` and `Nelder_Mead`:

```text
Downdated VtV is not positive definite
```

The recorded error identifies a numerical factorization failure, not its
ultimate cause for each feature. We have not established that all failures
are caused by low variation. They remain excluded from inference; they are
not silently replaced by slice-level GLM results. These are feature-specific
model failures, not a failure of the upstream SpatialESS communication run.

## GLM comparison and audit

The slice GLM has 4,172 BH-significant features; the LMM has 5,147, including
4,108 shared (Jaccard 0.788332). Across 39,570 paired valid fits, estimate
Spearman rho is 0.996623 and direction agreement is 99.2823%. This assesses
sensitivity to the dependence model; neither model is experimental truth.

The audit verifies 50,040 unique features, 45 slices/35 patients, global BH
recomputation, excluded q-values, failed p-values, retention of all optimizer
attempts and unchanged slice-GLM estimates/statuses. All eight checks passed.

See [LMM_RESULTS_20260907.md](LMM_RESULTS_20260907.md) for the table index,
source hash, exact parameters and reproduction commands. This cohort did not
undergo SPARKLE correction; copies in the SPARKLE repository validate only
its bundled downstream engine.
