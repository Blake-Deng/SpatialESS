# Reviewer guide: patient-level model diagnostics

This page explains how to read the patient-aware multi-sample analysis without
confusing a statistical boundary case with a software failure.

## Analysis unit and model

SpatialESS first computes one communication profile for each spatial slice.
The cohort analysis uses slices as rows and patient as the clustering variable:

```text
log1p(score) ~ condition + (1 | patient)
```

Cells are never treated as independent patient replicates. Multiple slices
from the same patient remain in the model, while `(1 | patient)` accounts for
their shared baseline. The slice-level model
`log1p(score) ~ condition` is reported as a transparent comparator.

## Frozen GSE250346 status table

The counts below are the previously frozen baseline run, before the additional
`Nelder_Mead` retry and the new support-screen fields were added. The robust
release keeps this baseline for traceability, but a manuscript using the new
implementation must regenerate the table and report the new status counts.

The frozen run contains 50,040 sender-receiver-LR features from 45 slices and
35 patients. The status counts are:

| Status | Features | How it is used |
|---|---:|---|
| `ok`, non-singular | 21,934 | Primary valid fit; included in BH correction |
| `ok`, singular | 17,308 | Valid boundary fit; included in BH correction and reported separately |
| `convergence_warning` | 2,299 | Coefficients retained for diagnosis; excluded from BH correction |
| `numerical_failure` | 305 | All configured optimizers failed; no inferential result |
| `insufficient_nonzero` | 8,194 | Not fitted because the feature was present in too few slices |

The primary LMM result therefore contains 39,242 fits with `status = "ok"`.
The complete counts are in
`results/tables/gse250346_lmm_fit_status.tsv`.

## What `singular` means

`singular = TRUE` means that the estimated patient random-intercept variance
is on the boundary, usually zero. In practical terms, that feature does not
provide evidence that patients need an additional random-intercept variance
component after the fixed effects are included. The fixed-effect estimate is
still estimable, so this is a successful boundary fit, not a crash and not a
claim that patients are identical. The non-singular subset is reported as a
sensitivity analysis.

## What `convergence_warning` means

The optimizer returned coefficient estimates, but `lme4` reported a numerical
diagnostic such as a large gradient, a non-positive Hessian eigenvalue, a
scale warning or a numerically singular Hessian. These warnings indicate that
the uncertainty calculation for that feature is not sufficiently reliable for
confirmatory FDR reporting. The estimates are preserved in the output so the
diagnostic is auditable, but the feature is conservatively excluded from
multiple-testing correction.

The 2,299 warnings were classified as:

| Diagnostic | Features |
|---|---:|
| Gradient warning | 1,447 |
| Negative Hessian eigenvalue | 543 |
| Scale warning | 241 |
| Numerically singular Hessian | 68 |

These categories are diagnostic labels, not additional biological findings.
The detailed table is
`results/tables/gse250346_lmm_convergence_summary.tsv`.

## What `numerical_failure` means

For 305 features, both the primary `bobyqa` fit and the `nloptwrap` fallback
failed with:

```text
Downdated VtV is not positive definite
```

This is a numerical failure for an individual mixed-model design, typically
caused by near-zero or highly redundant feature variation in the slice-level
response. The robust implementation now also tries `Nelder_Mead` and records
all optimizer attempts. If all configured optimizers fail, the feature is
marked `numerical_failure`. It is not an assertion that the input dataset or
SpatialESS engine failed. No p-value is reported for these features. The exact
message is retained in
`results/tables/gse250346_lmm_failure_summary.tsv`.

`fit_failed` is reserved for an unrecognized error after all configured
optimizers fail.

## Multiple-testing rule and comparison

BH-adjusted values are calculated only for `status = "ok"` fits. This prevents
uncertain numerical fits from being presented as confirmatory discoveries.
The primary LMM found 5,158 FDR-significant features; the slice-level GLM
found 4,172. Of these, 4,109 were shared (significant Jaccard 0.787), effect
estimates had Spearman rho 0.9967 and directions agreed for 99.28% of valid
comparisons. This comparison is a robustness and dependence-structure check,
not a claim that either model is experimental ground truth.

## Reproduce and audit

Run the cohort driver after supplying the prepared sample-level inputs:

```bash
Rscript benchmarks/run_gse250346_lmm_comparison.R
```

The driver records the model formula, feature-level status, optimizer,
variance components, maximum gradient, convergence messages and fit errors.
The exact input locations and parameter locks are listed in
`docs/REPRODUCE_MAIN.md`.
