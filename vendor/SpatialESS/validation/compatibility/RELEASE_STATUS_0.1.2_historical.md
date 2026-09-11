# SpatialESS 0.1.2 release status

Local validation passed. Ready for user upload and independent testing.
GitHub cloud CI: **NOT RUN**. No push, release or tag was created.
Linux/R 4.5.3 was tested; other operating systems and R versions were not run here.

## Completed checks

- Source installation into an isolated library: passed for 0.1.1 and 0.1.2.
- Full unit suite, including panel coverage and GLM/LMM: passed.
- Final `R CMD build` and `R CMD check --no-manual --no-build-vignettes`: Status OK.
- All R scripts parse; workflow YAML and actionlint v1.7.12: passed.
- The r-lib check action receives an R character vector for `args`; its R expressions were checked locally.
- Real cervical 1k mini (100 permutations): passed against full-precision baseline.
- Optional SpatialCellChat V3 molecular helper mini: passed (not full official inference).
- GSE250346 existing 50,040-feature LMM result bundle: passed its integrity verifier.

## Numerical regression

11 real-data inputs, 33 isolated runs: all completed, reaching 443,515 CosMx cells.
22 old/new comparisons: records, prepared input contracts, parameters/diagnostics
and permutation RNG were identical; maximum probability difference = 0,
exact p-value fraction = 1, active Jaccard = 1, significance disagreements = 0.
All C++ sources and the multi-patient model implementation are unchanged.

6 comparisons to existing official V3 RDS outputs also passed. Maximum
probability difference = 7.8062556418956304e-18; exact p-value fraction = 1,
active Jaccard = 1 and significance disagreements = 0. Original official and
dimension-patched official results are labelled separately. These official
outputs are archived runs, not newly timed 0.1.2 experiments. Official 100k+
failures remain failures; no fidelity claim is made at those official-failed sizes.

## Files and boundaries

- `validation/lr_interface/new_vs_old.tsv`: strict version regression.
- `validation/lr_interface/resources.tsv`: runtime, RSS, completion and parameters.
- `validation/lr_interface/new_vs_archived_official.tsv`: historical V3 comparisons.
- `validation/lr_interface/logs/`: full successful checks and earlier diagnostic logs.
- `validation/lr_interface/provenance.json`: source versions, threads and environment.
- `validation/lr_interface/input_output_checksums.tsv`: full-run input/output hashes.
- `validation/fixtures/`: real-data mini; `validation/baseline/`: old source tarball.
- `SHA256SUMS.txt`: checksums for this complete GitHub tree.

The initial mini compared rounded TSV numbers and failed its 1e-15 threshold;
the RDS baseline confirmed exact equality. This was fixed by using RDS in the
mini, and the initial failure log is retained. Documentation generation first
failed because the conda compiler was not on PATH; generation, clean compilation
and final checks subsequently passed with the correct environment.

The coverage feature has preparation/audit overhead. Timings are single runs on
a shared server; they do not establish unchanged end-to-end speed. Numerical
fidelity is not biological accuracy. LMM statuses remain model diagnostics;
the existing 292 numerical failures and 1,984 convergence warnings were not
relabelled as successful results by this release.

Upload the contents of the ZIP's `SpatialESS/` directory to the repository root,
including `.github/`, then inspect the actual GitHub Actions result. Uploading
only the ZIP file will not update the repository's package code or trigger an
appropriate code check.
