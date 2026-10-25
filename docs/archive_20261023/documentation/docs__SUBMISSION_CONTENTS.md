Current public package mapping and publication status: [REPRODUCE_20261023.md](REPRODUCE_20261023.md). The historical inventory below describes earlier versions and remains useful for path provenance.

# Manuscript data and verification records

These materials accompany SpatialESS 0.1.4. The engine and model source match
revision 626ffecdd4cd632a462a037edd903cf67d062c54. The 2026-10-22 update assembles
additional saved experiments and current figure masters, audited against public
main `f92dff98f9019300f8fe5086d8a011eed4f089d8`; it does not rerun inference,
permutations or cohort fits. Manuscript Word drafts are retained locally for
author review and are not part of this public data snapshot.

Start with [the current manuscript-materials index](MANUSCRIPT_MATERIALS_20261022.md).

- `paper/figure_package/`: main data figures, S1-S4, vector PDF/SVG, plotting
  scripts, source workbook and per-panel CSV files.
- `paper/tables/`: Supplementary Table 1 in TSV, CSV and XLSX formats.
- `paper/provenance/`: archived ovarian pairs, bridge comparisons, RCTD labels
  and scripts, lung metadata mapping evidence, input hashes and verification logs.
- `paper/atera_validation/`: Atera preview-data preparation, exact sampled IDs,
  database/parameter contracts, scale comparisons, successful and failed attempts.
- `paper/ablation/`: component variants, 170 validated technical runs, timings,
  peak RSS, numerical verification and source snapshots.
- `paper/lung_cohort/`: complete 343-gene primary models/scores/metadata and a
  separate 99-gene official-reference archive.
- `results/`: current 0.1.4 inference records and archived cohort fits, with
  their original configuration and measurement scopes.

Install the R package from the repository root:

```sh
R CMD INSTALL .
```

Regenerate current main figures from their frozen tables (no inference):

```sh
python -m pip install -r paper/figure_package/updates/Figure1_Atera_20261016/requirements.txt
python paper/figure_package/updates/Figure1_Atera_20261016/scripts/redraw.py --figure1-only
python paper/figure_package/updates/Figure1_Atera_20261016/scripts/compose_figure1.py
python paper/figure_package/updates/Figure2_20261012/scripts/redraw.py --figure2-only
```

Verify saved historical pairs and current ovarian linkage, writing this run
into a separate output directory:

```sh
Rscript paper/provenance/verify_archived_records.R . /tmp/SpatialESS-record-verification
python paper/provenance/ovarian/recheck_label_transfer.py
python paper/provenance/github_update_20261022/verify_addendum.py
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
also rerun after the Figure 1 presentation correction. The revised workflow and
export checks are documented in `paper/figure_package/reports/QA_REPORT.md`.

The current Figure 1 proportionally places the original author workflow without
internal edits; older workflow replacements are historical exports. Main Figure
2D uses I15/I06/I19, whereas S4D retains I03/I21/I19. Use the current figure
labels and full feature-key mappings.

The original RCTD annotation sessionInfo was not saved. The archived DESCRIPTION
identifies spacexr 2.2.1, which does not substantiate the draft's 2.2.16. See the
ovarian environment note before naming an exact original-run version. A
The historical 2026-10-22 Zenodo archive was published as version DOI 10.5281/zenodo.23567248;
concept DOI 10.5281/zenodo.23567247.
