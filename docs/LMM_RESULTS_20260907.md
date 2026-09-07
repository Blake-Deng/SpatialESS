# GSE250346 robust LMM results

The [complete result index](../vendor/SpatialESS/docs/LMM_RESULTS_20260907.md)
documents formulas, parameters, source hashes, statuses and limitations.

- Convenient table location: [benchmark/gse250346_lmm](../benchmark/gse250346_lmm/).
- Identical bundled-engine tables: [vendor/SpatialESS/results/tables](../vendor/SpatialESS/results/tables/).
- Diagnostic interpretation: [REVIEWER_GUIDE.md](REVIEWER_GUIDE.md).

Verify from this repository root with base R:

```bash
Rscript vendor/SpatialESS/benchmarks/verify_gse250346_result_bundle.R \
  benchmark/gse250346_lmm
```

There are 45 slices, 35 patients, 39,570 valid LMM fits, 1,984 warning fits,
292 numerical failures and 8,194 features without sufficient nonzero support.
Only valid fits enter BH correction; 5,147 have q < 0.05 under the stated
Wald-normal approximation.

This is not a SPARKLE-corrected cohort. The result demonstrates downstream
patient-aware analysis using the bundled SpatialESS code; it does not test
whether SPARKLE correction improves biological truth.
