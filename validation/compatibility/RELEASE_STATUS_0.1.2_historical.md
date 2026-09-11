# SpatialESS + SPARKLE 0.1.2 validation status

Locally validated and packaged for independent testing. GitHub cloud CI: **NOT RUN**.
No GitHub push, tag or release was created. Tested on Linux with R 4.5.3.

- Bundled engine: every file in the frozen main 0.1.2 tree matches its counterpart; see `validation/bridge_interface/engine_sync.tsv`.
- Clean engine and bridge installation: passed.
- Bridge R tests: 16 passed, zero failures/warnings/skips; Python bridge tests: 2 passed.
- Bridge `R CMD build/check --no-manual --no-build-vignettes`: **Status OK**.
- The unchanged main engine already passed its independent source check; its full logs are bundled under `vendor/SpatialESS/validation/`.
- R syntax, workflow YAML and actionlint checks: passed. Cloud execution is not inferred from these checks.
- Synthetic 24-cell demo and real RAW/corrected ovarian 1k mini: passed.
- Command-line workflow with the newer SpatialCellChat database: passed; the `PYTHON` override selects the interpreter when R and Python use separate environments.
- Preserved 50,040-feature LMM result bundle: integrity verification passed; existing warning/failure statuses were not reclassified.

## Numerical regression

RAW and corrected H5AD, each at 1,000 / 5,000 / 16,247 cells: **18/18 isolated processes completed**.
All **12/12 strict comparisons** passed: old bridge versus updated bridge, and updated bridge versus direct main.
Full records, expression/metadata inputs, result contracts and permutation RNG agree exactly.
Maximum probability difference = **0**; exact p-value fraction = **1**; active Jaccard = **1**; significance disagreements = **0**.

Three comparisons to archived full-data SpatialESS / official V3 outputs also passed the 1e-15 probability threshold,
with exact p-values and zero significance disagreements. Maximum absolute probability difference across these
historical comparisons is **1.32923e-16**. Historical text-prepared matrices and the current H5AD export can differ
at floating-point precision; these are not described as bitwise identical.

Tables: `validation/bridge_interface/comparisons.tsv`, `resources.tsv`, `history.tsv`.
The historical ovarian database is frozen (1,939 candidates, 1,828 computable LR); the newer default V3 database
is a separate version. See `docs/LR_INTERFACE_20260910.md` for all parameters and comparison boundaries.

The SPARKLE correction model was not modified or refitted. This update reuses its existing corrected H5AD files.
The Python conversion and all-gene normalization logic are unchanged. Numerical fidelity does not establish
greater biological accuracy of correction. Resource measurements are single runs on a shared server.

The first mini invocation occurred before its database fixture had finished copying and failed with a missing-file
error. Its diagnostic log is retained; the complete bundled mini subsequently passed for both conditions.
The first CLI smoke test selected an R-environment Python without anndata. The CLI now accepts the `PYTHON`
environment override, and a repeat with the intended interpreter passed. This did not change inference code.

Upload the contents of `SpatialESS-SPARKLE/` to the integration repository root, including `.github/`.
The bundled engine avoids a separate main checkout. Python/SPARKLE and R dependencies still need installation.
