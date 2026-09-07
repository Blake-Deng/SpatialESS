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
`fit_multisample_communication_lmm()` using `lme4::lmer()`. Each feature is fit separately. Responses are standardized for numerical conditioning and estimates are returned on the original log1p-score scale. The fitter uses `bobyqa`, retries failures with `nloptwrap`, and records the optimizer, variance components, maximum gradient, convergence message, fit error and singular-fit flag. P-values use the stated Wald normal approximation, and BH correction is applied only to fits with `status = "ok"`.

## Status interpretation

- `ok`, `singular = FALSE`: converged model with nonzero patient variance.
- `ok`, `singular = TRUE`: converged boundary model; patient variance is estimated as zero.
- `convergence_warning`: diagnostic coefficients are retained, but the feature is excluded from FDR.
- `fit_failed`: both optimizers failed.
- `rank_deficient`: the requested fixed effect was not estimable.
- `insufficient_nonzero`: the feature was present in too few slices to fit.

Singular fits are successful fits, not software failures. The frozen GSE250346 run contained 39,242 `ok` fits (21,934 non-singular and 17,308 singular), 2,299 convergence warnings, 305 failures after both optimizers and 8,194 features below the nonzero-slice threshold. The LMM identified 5,158 FDR-significant features, compared with 4,172 for the slice-level GLM; 4,109 were shared (Jaccard 0.787). Valid-fit effect estimates had Spearman rho 0.9967 and 99.28% direction agreement. The main report should present all status counts and include a non-singular sensitivity subset.

## Reproduce on GSE250346

The full sample-level result files and metadata are expected at
`/data/dzf/GSE250346/benchmarks/multisample/`. Run:

```bash
Rscript vendor/SpatialESS/benchmarks/run_gse250346_lmm_comparison.R
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
