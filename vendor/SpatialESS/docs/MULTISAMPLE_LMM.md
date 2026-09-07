# Multi-patient communication models

SpatialESS reports communication scores for each spatial sample. The sample,
not an individual cell, is the observational unit for cohort-level inference.
When a patient contributes multiple slices, the primary model keeps those slices as separate rows and accounts for their shared patient identity with a random intercept.

## Two models

The slice-level baseline is:

```text
log1p(score) ~ condition
```

It treats slices as independent observations. It is useful as a transparent
baseline, but it can underestimate uncertainty when multiple slices come from
the same patient.

The patient-aware model is:

```text
log1p(score) ~ condition + (1 | patient)
```

This is implemented by
`fit_multisample_communication_lmm()` using `lme4::lmer()`. Each feature is fit separately. Responses are standardized for numerical conditioning and estimates are returned on the original log1p-score scale. Before fitting, the function records zero-score fraction and response variation, checks the fixed-effect design rank, and screens features with insufficient support. The fitter uses `bobyqa`, `nloptwrap` and `Nelder_Mead` in sequence, preferring a warning-free fit, and records the optimizer, variance components, maximum gradient, convergence message, fit error and singular-fit flag. P-values use the stated Wald normal approximation, and BH correction is applied only to fits with `status = "ok"`.

## Status interpretation

- `ok`, `singular = FALSE`: converged model with nonzero patient variance.
- `ok`, `singular = TRUE`: converged boundary model; patient variance is estimated as zero.
- `convergence_warning`: diagnostic coefficients are retained, but the feature is excluded from FDR.
- `insufficient_zero_support`: zero-score fraction exceeded `max_zero_fraction`.
- `insufficient_variation`: response standard deviation was below `min_response_sd`.
- `numerical_failure`: all configured optimizers failed with a recognized numerical error.
- `fit_failed`: all configured optimizers failed without a recognized numerical error.
- `rank_deficient`: the requested fixed effect was not estimable.
- `insufficient_nonzero`: the feature was present in too few slices to fit.

Singular fits are successful fits, not software failures. Features rejected by
the support screens are also not treated as model failures. Numerical failures
are retained with their complete optimizer messages and never receive p-values
or q-values. The frozen GSE250346 result tables should be regenerated after
changing the optimizer or support parameters; the release records the exact
parameters used for each result.

The reviewer-oriented explanation of these distinctions is in
`docs/REVIEWER_GUIDE.md`.

## Reproduce on GSE250346

The full sample-level result files and metadata are expected at
`/data/dzf/GSE250346/benchmarks/multisample/`. Run:

```bash
Rscript benchmarks/run_gse250346_lmm_comparison.R
```

The script writes:

```text
slice_glm_vs_patient_lmm/slice_level_glm.tsv
slice_glm_vs_patient_lmm/patient_random_intercept_lmm.tsv
slice_glm_vs_patient_lmm/model_summary.tsv
slice_glm_vs_patient_lmm/lmm_status_summary.tsv
slice_glm_vs_patient_lmm/lmm_failure_summary.tsv
slice_glm_vs_patient_lmm/lmm_convergence_summary.tsv
```

The GLM is the slice-level baseline. The LMM is the primary patient-aware model. For confirmatory inference, degrees of freedom can be
replaced by `lmerTest` or a parametric bootstrap; the current implementation
keeps the dependency limited to `lme4` and makes the approximation explicit.
