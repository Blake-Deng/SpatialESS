# SpatialESS 0.1.3 Local Validation Status

Local numerical acceptance and package checks passed. Ready for independent testing.
GitHub cloud CI: **NOT RUN**. No push, tag or GitHub release was created.

## Completed in this build

- Frozen delivered 0.1.2 ZIPs retained outside this repository with SHA256 manifests.
- Cervical center 1k/5k and submitted stratified 1k/5k, each with 0.1/5 and 0/1
  filters: **8/8 strict comparisons passed**, 1076 shared LR, 100 permutations.
- Active Jaccard=1; maximum probability difference=1.94289029309402e-16; exact p-value
  fraction=1; significance disagreements=0 at p<0.05 in all eight cases.
- Shared ALRA+MERINGUE 1k/596 LR against frozen official: passed, 232 active records.
- Fresh full official 1k run using the portable mini script: passed, 57 active records.
- Historical CosMx 1k/5k/10k/50k, filters 0/1: 4/4 old/new comparisons passed;
  full records, input, RNG and contracts (except the new filter marker) identical.
- The same four candidate runs also passed against checksum-verified archived
  original official V3 records; these are not freshly timed official runs.
- Source install, full R tests, LR coverage mini, R syntax and actionlint: passed.
- Main and standalone bridge R CMD build/check: **Status: OK**, Linux/R 4.5.3.
- Existing 50,040-feature GSE250346 LMM bundle: integrity check passed. Model code
  and cohort results were not changed or refitted.

## What changed and what remains open

ESS now follows the pinned V3 formatted-fraction filter in both observed and
permutation inference. The old default-filter discrepancy is recorded, not
relabelled as a pass. Original official dimensions fail at 927x927 for the
submitted 1k; only adding dims restores 1000x1000 with unchanged distance entries.
Both original failure and patched successes are retained.

The submitted file contains barcodes only. Tests use public graph clusters,
not the collaborator's unavailable annotation. Their original all-zero session
is **UNRESOLVED** without the actual expression/object and downstream calls.
No claim is made that this release fixes every possible zero-result cause.

Historical 100k+ studies and cohort tables remain version-scoped evidence,
not new full-cohort 0.1.3 reruns. Existing LMM warnings and numerical failures
retain their labels. Numerical fidelity is not biological accuracy.

## Evidence and reproduction

`validation/compatibility/dual_track.tsv`, `historical_replay.tsv`,
`shared_alra_meringue.tsv`, `fresh_official_mini.tsv`, `resources.tsv`,
`source_changes.tsv`, `input_output_checksums.tsv` and `logs/` contain the audit.
Resource logs are diagnostic shared-server measurements, not a clean speedup
experiment. See `docs/V3_COMPATIBILITY.md` for the pinned reference and portable
frozen-reference or fresh-official reproduction commands.

Upload the **contents** of the ZIP's `SpatialESS/` folder, including `.github/`,
to the repository root. Uploading a ZIP alone does not update package code or CI.

Fresh ZIP extraction, checksum validation, source installation and minis/tests passed.
The standalone installer selected its bundled vendor engine. Logs: `validation/compatibility/zip_verification/`.
