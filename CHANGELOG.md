# Changelog

## Unreleased

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
