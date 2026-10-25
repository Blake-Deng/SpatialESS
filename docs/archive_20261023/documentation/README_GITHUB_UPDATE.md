# SpatialESS GitHub update bundle (2026-10-22)

This is an add-on for `Blake-Deng/SpatialESS`, audited against public `main`
commit `f92dff98f9019300f8fe5086d8a011eed4f089d8` (2026-10-12; package
version 0.1.4). It is staged for review and manual upload; it does not modify
the remote repository or the local SpatialESS working tree.

## What is included

- `paper/figure_package/updates/Figure1_Atera_20261016/`: the latest Figure 1
  master, panels, source workbook, plotting code, settings and QA reports. The
  author-provided workflow asset is retained unchanged. Panel I shows the
  200k and 700k categories; the 700k SpatialCellChat run is an attempted
  allocation failure, not an unattempted reference. The 400k and nominal 450k
  categories are omitted from panel I as specified, while their source records
  remain documented.
- `paper/figure_package/updates/Figure2_20261012/`: the revised main Figure 2,
  including the separated panel C annotation headers, I15/I06/I19 panel D
  examples, condition-specific nonzero-section counts, plotting code and sources.
- `paper/atera_validation/`: Atera preview-data provenance, QC and sampling
  rules, exact barcode lists, fixed LR inputs, parameter contracts, scripts,
  tables, environment records, successful runs and failed attempts. Raw
  expression matrices are not included; source URLs and file checksums are.
- `paper/lung_cohort/primary_343gene37LR/`: the complete saved 50,040-feature
  model, nonzero slice scores, fixed feature universe, metadata, cell-type
  counts, pDC counts, original model code, checksums and session information.
  Missing score rows are zero only within the fixed 45-section by 50,040-feature
  universe, as documented in the included cohort README.
- `paper/lung_cohort/reference_99gene_validation/`: the separate official
  reference comparison using 99 signaling genes. It must not be merged with
  the 343-gene primary disease model.
- `paper/ablation/`: the reproducible spatial-hashing, graph-representation,
  edge-score-reuse and dense-baseline ablation tables, configurations,
  scripts, run records and validation reports.
- `paper/provenance/ovarian/`: the archived `spacexr` DESCRIPTION evidence
  and a note on the missing original RCTD session information.

No inference, model fit, significance test or source-data transformation was
run to assemble this bundle. Files are copied from the frozen analysis and
figure outputs; checksums are provided in `SHA256SUMS.txt`.

## Upload map

Copy the contents of this bundle into the repository root, preserving paths.
The new Atera, ablation and lung directories are additive. Updated README/index,
ignore and checksum files are also supplied for the exact audited public commit.
Main Figure 2 on that commit still has I03/I21/I19 examples, so the revised main
Figure 2 is included. Supplementary Figure 4 already matches the review package
and keeps its I03/I21/I19 interval examples; its files and the existing cohort
sensitivity analysis should be retained.

The supplied `paper/figure_package/README.md` identifies the dated Figure 1 and
Figure 2 directories as current masters, with direct links to figures, source
data and reproduction commands. Earlier figure files remain historical outputs;
the new index prevents treating the old `figures/Figure2.pdf` or the original
all-figure runner as the current main figure. Do not overwrite the author workflow
source asset.

Add the included `paper/provenance/ovarian/RCTD_ENVIRONMENT_NOTE.md` beside the
existing ovarian provenance README. The repository already contains a 2.2.1
DESCRIPTION and currently documents that no original RCTD `sessionInfo()` was
saved; do not replace this with 2.2.16 without a matching archived package
record.

## Already present in GitHub

The public repository already contains SpatialESS 0.1.4, Supplementary Figure 4
assets, the existing cohort-sensitivity package,
the patient-level model summaries, ovarian label-transfer artifacts, and the
statement that the original RCTD session information is unavailable. Those
files are not duplicated here.

No engine source, package version or statistical method is changed by this
update. The archived ablation variants are experimental sources under
`paper/ablation/source/`, not replacements for the installable engine.

## Verify after merging

Use an independent checkout of the audited commit and copy the bundle contents
into it. The updated `.gitignore` admits the compact contracts/database RDS
files in these specific new paper directories; large expression inputs remain
excluded. Then run from the repository root:

```sh
python3 paper/provenance/github_update_20261022/verify_addendum.py
python3 paper/verify_manifest.py
sha256sum --check --quiet SHA256SUMS.txt
git status --short
```

The root manifests cover the merged audited repository. The addendum's own
`paper/provenance/github_update_20261022/ADDENDUM_SHA256SUMS.txt` can also be
checked before merging. If public main changes after the audited commit,
reconcile its changes before replacing README/index/ignore/checksum files.

The TSTA/SecAct exploration is not included because it is outside the current
Figure 1/Figure 2 manuscript evidence. Add it only if the manuscript makes a
claim based on that analysis. The submission Word documents are also kept out
of this code-and-data addendum; publish them separately only when intended.

## Verified interpretation boundaries

- Atera is 10x WTA FFPE **preview** data. The package documents fixed-LR
  comparisons, not the full overexpression-based LR selection workflow and not
  a biological validation cohort.
- The Atera reference failure at 715,512 cells is a saved allocation failure
  for a 96.1 GiB dense intermediate under the recorded 96 GiB virtual-address
  guard; it does not establish a universal SpatialCellChat limit.
- The lung 99-gene reference set has 42 evaluable pairs: 41 strict passes and
  one `VUILD142MA` permutation p-value difference of 0.05 without a changed
  significance call. Three other reference attempts failed. These are not
  343-gene per-feature comparisons.
- For the selected SPP1-CD44 feature, saved nonzero support is 29/35 PF and
  0/10 control sections. Under the descriptive rule “at least one nonzero
  section”, patient support is 22/26 PF and 0/9 control patients. These counts
  are support summaries, not a prespecified condition-consensus definition;
  control zeros do not establish biological absence.
- The ovarian archive records `spacexr` 2.2.1, not 2.2.16. The original
  annotation `sessionInfo()` is unavailable.
