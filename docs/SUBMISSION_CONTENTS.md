# Manuscript data and verification records

These materials accompany SpatialESS 0.1.4. The engine and model source match
revision 626ffecdd4cd632a462a037edd903cf67d062c54. No new communication inference,
permutations or cohort fits were run for this update. Manuscript Word drafts are
retained locally for author review and are not part of this public data snapshot.

- `paper/figure_package/`: main data figures, S1-S4, vector PDF/SVG, plotting
  scripts, source workbook and per-panel CSV files.
- `paper/tables/`: Supplementary Table 1 in TSV, CSV and XLSX formats.
- `paper/provenance/`: archived ovarian pairs, bridge comparisons, RCTD labels
  and scripts, lung metadata mapping evidence, input hashes and verification logs.
- `results/`: current 0.1.4 inference records and archived cohort fits, with
  their original configuration and measurement scopes.

Install the R package from the repository root:

```sh
R CMD INSTALL .
```

Regenerate figures from frozen tables (no inference):

```sh
python -m pip install -r paper/figure_package/requirements.txt
python paper/figure_package/scripts/run_all.py
```

Verify saved historical pairs and current ovarian linkage, writing this run
into a separate output directory:

```sh
Rscript paper/provenance/verify_archived_records.R . /tmp/SpatialESS-record-verification
python paper/provenance/ovarian/recheck_label_transfer.py
python paper/verify_manifest.py
```

The label verifier requires pandas. The R record verifier uses base R. The
figure package's original staging scripts document assembly from the original
project; they are archival provenance, not portable entry points. The commands
above operate on the included records and source tables.

`SUBMISSION_SHA256.tsv` checks file sizes and SHA256 hashes, including the current
`SHA256SUMS.txt`. To avoid circular checksums, `SHA256SUMS.txt` covers packaged
files except itself and `SUBMISSION_SHA256.tsv`. The original 1 October source
manifest is preserved at `paper/provenance/release_checks/SHA256SUMS_20261001.txt`.

The 66 current official comparisons exclude ovarian and use the 99-gene lung
signaling preparation, not the 343-gene cohort preparation. The cohort fit date
is 7 September 2026. Current resource measurements and archived reference values
retain their separate scopes. Valid singular models remain eligible for the
primary BH correction; warnings and failures are retained in the result tables.

Fresh macOS R 4.5.2 checks and saved-record verification are documented in
`paper/provenance/editorial_checks_20261012/`. The original Linux check and
verification environments remain archived separately. The figure pipeline was
also rerun: 64 input hashes, 50 data checks and 58 export checks passed, with the
known Figure 1A lettering warning retained.

The original RCTD annotation sessionInfo was not saved. Figure 1A remains an
unfinished reference composition: minimum lettering is approximately 2.59 pt
in the compact slot, and the data master leaves A empty. Neither is a finished
submission figure. A persistent archive DOI has not yet been assigned.
