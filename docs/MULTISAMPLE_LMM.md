# Multi-patient communication models

The bundled SpatialESS engine fits
`log1p(score) ~ condition + (1 | patient)` using `lme4::lmer()` and retains
multiple slices per patient. A slice-level GLM is the comparator, not a
fallback substituted for failed LMM features.

See the complete [model documentation](../vendor/SpatialESS/docs/MULTISAMPLE_LMM.md)
and [diagnostic guide](REVIEWER_GUIDE.md). The robust GSE250346 rerun is
complete and indexed in [LMM_RESULTS_20260907.md](LMM_RESULTS_20260907.md).

From the repository root, use `vendor/SpatialESS` as the source-repository
argument in the checkpointed driver:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
Rscript vendor/SpatialESS/benchmarks/run_gse250346_lmm_checkpointed.R \
  vendor/SpatialESS /path/to/multisample /path/to/new_output 4
```

The required prepared cohort inputs are described in the result index.
GSE250346 is downstream-engine validation, not a SPARKLE-corrected cohort.
