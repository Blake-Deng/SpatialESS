# SpatialESS 0.1.3

- Match fixed official V3 percentage formatting in observed and permutation
  group filters. Outputs near nonzero min.percent boundaries may change from
  0.1.2; the former mismatch is retained as failure evidence.
- Add an official-anchored cervical regression fixture and format/threshold
  tests. CSR geometry and molecular scoring are unchanged.
- Document shared ALRA/MERINGUE preparation, distinct distance-patched official
  results, and the difference between missing p-values and valid zero results.

# SpatialESS 0.1.2

- Add `filter_cellchat_lr()` with complete mandatory-subunit checks, explicit
  exclusion warnings, retained/excluded tables, and missing-gene provenance.
- Add opt-in `spatialess(..., lr_missing = "filter")`; preserve strict default.
- Improve incomplete-LR and empty-input diagnostics and document targeted panels.
- Preserve C++ inference, global expression normalization, permutation RNG, and
  the existing multi-patient GLM/LMM implementation.
