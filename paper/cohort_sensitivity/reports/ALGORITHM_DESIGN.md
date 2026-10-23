# SpatialESS cohort robustness layer: design and tested prototype

The computational engine remains SpatialESS 0.1.4. The proposed addition consumes
its saved group-level communication scores and makes cohort inference and its
limitations explicit. It does not redefine CellChat probabilities, permutations,
spatial graphs or self-interactions.

## 1. Input contract and identification checks

- One row per slice; one column per full sender/receiver/LR key.
- Align scores, metadata and model tables by full keys, never by LR names alone.
- Preserve absent-feature zeros and all original cell annotations.
- Check condition consistency within patient; record independent patient counts.
- Audit design-matrix rank before fitting. Refuse non-identifiable adjustments;
  do not silently remove a covariate or call it adjusted.
- In this cohort, condition + run and condition + TMA are identifiable separately.
  condition + sample_affect is not identifiable, and run + TMA is redundant.
  sample_affect is an affectedness label, not verified anatomical-region metadata.

## 2. Preserve the specified primary analysis

Fit the original maximum-likelihood random-intercept model to log1p(slice score).
Return the original-unit condition coefficient, SE, fit state, singular flag,
optimizer diagnostics and global BH q. Retain warning-free singular boundary fits
under the documented primary rule. Exclude convergence warnings and failures from
primary q values. Singularity is a diagnostic; it is not suppressed by increasing
an optimizer limit or automatically deleting an inconvenient feature.

## 3. Require an independent-patient sensitivity result

For patient i and feature j:

    y_ij = mean over retained slices s of log1p(score_sj)

This includes zeros and is not log1p(mean(score)). It gives each of the 35 patients
one observation; it does not increase the number of controls beyond nine.

The tested sensitivity model is:

    y_ij = intercept + beta_j * PF_i
           + gamma2_j * fraction_of_slices_in_Run2_i
           + gamma3_j * fraction_of_slices_in_Run3_i + error_ij

Run fractions retain the mixed-run information for a patient with slices from two
runs. Choosing the first slice's run would discard that information. The reference
is Run1. With 35 patients and four design columns, the residual df is 31.

Use HC3 heteroskedasticity-consistent covariance:

    V_HC3 = (X'X)^(-1) X' diag(e_i^2 / (1 - h_ii)^2) X (X'X)^(-1)

Two-sided t-based p values use n_patient - rank(X) df; BH is computed over the full
eligible feature universe. This is approximate small-cohort sensitivity inference,
not an exact finite-sample guarantee. Unit leverage, rank deficiency and inadequate
support produce explicit unavailable states. No cell-level pseudo-replication is
introduced. A separate patient Welch comparison checks the unadjusted contrast.

## 4. Diagnose imbalance and patient influence

- Refit patient HC3 after each patient omission; recompute full-universe BH each
  time. Save signs, effect ranges and retained discoveries, including unfavourable
  results.
- Refit run-adjusted patient effects in 1,000 9-Control versus 9-PF subsets. Keep
  all nine controls and sample nine PF patients without replacement.
- These subset percentiles summarize effect stability under reduced PF sampling;
  they are not confidence intervals, p values, or independent validation.
- Preserve an audit of zero-inclusive score distributions and condition-specific
  nonzero-slice/patient support. Low support does not become biological absence.

## 5. Report evidence without changing the primary testing universe

The prototype adds descriptive agreement labels:

| Label | Rule |
| --- | --- |
| supported_by_run_and_patient_sensitivity | Primary, run LMM and patient HC3 are valid, all have their own full-universe BH q < 0.05, and signs match |
| primary_hit_sensitivity_dependent | Primary hit has valid sensitivities but does not satisfy that agreement rule |
| primary_hit_sensitivity_unavailable | Primary hit has an unavailable sensitivity fit |
| not_primary_hit | Valid primary feature without primary BH significance |
| primary_fit_unavailable | Primary fit is invalid or unavailable |

These labels are exploratory cross-analysis summaries, not a new statistically
calibrated FDR classifier. Do not call 3,808 agreeing features new discoveries or
independent validation. Do not rename the analysis as a novel biological model.
Original figures and their selected features remain traceable; the additional
results are supplied as Supplementary Figure 5 and Table 2.

## 6. Implemented API and verification

`scripts/robust_cohort_inference.R` provides `audit_cohort_designs()`,
`fit_patient_run_hc3()` and `classify_cohort_evidence()`. This is an experimental
analysis module supplied beside the paper, not a silently changed 0.1.4 package
API. The full 50,040-feature R implementation agrees with the independent Python
implementation on estimates, SE, intervals, p and q within 1e-10; the 28 original
plotted HC3 fits additionally agree with statsmodels.

## 7. What remains outside this algorithm

This layer cannot resolve confounded or missing anatomical-region data, recover
unmeasured genes, separate self/autocrine contributions under the retained model,
establish collagen deposition, or provide causal/experimental validation. More
controls, anatomical-region metadata and independent cohorts remain necessary for
stronger biological claims. A hurdle model would define a different estimand and
requires prespecified support/identifiability rules and separate validation; it is
not introduced here merely to obtain more significant results.
