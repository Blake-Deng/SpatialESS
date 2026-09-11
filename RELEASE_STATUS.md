# SpatialESS + SPARKLE 0.1.3 Local Validation Status

Locally validated for independent testing. GitHub cloud CI: **NOT RUN**.
No push, tag or release was created. Tested on Linux/R 4.5.3.

- Bundles exactly the repaired main 0.1.3 tree; no separate main checkout needed.
- Shares the same CSR inference, LR coverage interface and patient-random-intercept
  LMM implementation. See vendor/SpatialESS/RELEASE_STATUS.md for main acceptance.
- RAW and corrected ovarian H5AD, each at 1k/5k/16,247 cells: 12/12 new runs
  completed; 12/12 comparisons passed (old vs new, bridge vs direct main).
- Full records, prepared inputs, RNG and contracts (except the filter marker)
  identical. Probability difference=0, p-value exact fraction=1, active Jaccard=1,
  significance disagreements=0. Historical filters remain 0/1.
- Python tests, R tests, minimal demo and both real ovarian 1k minis: passed.
  Python tests emitted a SciPy deprecation warning, not a failed test.
- Engine and bridge R CMD build/check: Status OK. Syntax and actionlint passed.
- The existing LMM result bundle is preserved with all diagnostic statuses.

SPARKLE's upstream model and H5AD conversion were not changed or refitted.
Corrected-versus-RAW biological accuracy is not established by engine equality.
Old benchmark tables keep their original dates/version scope; no new official
full ovarian run or new full-cohort fit is claimed.

See validation/compatibility/bridge_replay.tsv, bridge_runs/, checks/ and logs/.
Standalone dependencies (R/Python/SPARKLE) still require installation.
Upload all contents of SpatialESS-SPARKLE/, including vendor/ and .github/.

Fresh ZIP extraction, checksum validation, source installation and minis/tests passed.
The standalone installer selected its bundled vendor engine. Logs: `validation/compatibility/zip_verification/`.
