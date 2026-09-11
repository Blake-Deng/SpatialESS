# Changelog

## 0.1.3 - 2026-09-11

- Bundle SpatialESS 0.1.3 with the official formatted-percentage filter in
  observed and permutation inference; retain strict complete-subunit auditing.
- Preserve SPARKLE correction and sparse H5AD conversion, and replay RAW and
  corrected ovarian inputs against 0.1.2 and direct main inference.
- Include the same cervical default-filter acceptance evidence as main and
  check both R packages in the standalone workflow.

## 0.1.2 - 2026-09-10

- Bundled the same SpatialESS 0.1.2 source as main, including GLM/LMM and validation evidence.
- Replaced silent per-LR exception handling with the complete-subunit filter and missing-gene audit.
- Added `lr_missing = "filter"` / `"error"`; retained `filter_unresolved` with conflict checks.
- Preserved H5AD conversion, all-gene normalization, C++ inference and permutation statistics.
- Updated command-line database loading to SpatialCellChat V3 and CI version handling.

## Earlier LMM integration

- Synchronized the bundled SpatialESS engine with the patient-random-intercept `lme4::lmer()` workflow.
- Added response scaling, optimizer fallback and explicit fit, singularity and convergence status fields.
- Added the GSE250346 slice-level GLM versus LMM reproduction script and compact status tables.
- Updated CI to install `lme4` and test both the bundled engine and H5AD bridge.

## 0.1.1 - 2026-09-06

- Rewired the bridge to the standalone SpatialESS package and renamed the public engine API.
- Preserved all-gene CP10K normalization denominators before gene selection.
- Added direct `spatialess_from_h5ad()` execution.
- Added external metadata support for labels and coordinates.
- Added a 24-cell reproducible example and bundle validation tests.
- Recorded the full ovarian SPARKLE official-versus-SpatialESS benchmark.
