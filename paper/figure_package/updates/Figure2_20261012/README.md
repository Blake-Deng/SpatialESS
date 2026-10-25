# Figure 2: patient-aware pulmonary-fibrosis score analysis

Panel C displays 28 selected features with support counts, patient ICC and
singular-fit annotations. Panel D uses examples I15, I06 and I19, matched by
complete feature keys. Supplementary Figure S4D uses I03, I21 and I19.
Model estimates and q values come from the saved primary cohort analysis.

## Files

Download and mount the companion data package using the
[reproduction guide](../../../../docs/REPRODUCE.md) to access:

- `figures/Figure2.pdf`, `.svg`, `.png`: figure masters and preview.
- `source_data/Figure2_Source_Data.xlsx`: source tables and feature mapping.
- `reports/QA_REPORT.md`: verification and visual inspection records.

## Reproduce

From the repository root, install `requirements-reproduce.txt`, then run:

```sh
python3 tools/redraw.py Figure2 --materials ../SpatialESS_reproduction_20261023
```

To extract example matrices from a complete saved score table, run from this
figure directory:

```sh
python3 scripts/prepare_example_sources.py --score-table /path/to/slice_scores.tsv.gz
```

The score input is a section-by-feature gzip TSV with `sample_id` and complete
feature keys as columns. A score table is included in
`paper/lung_cohort/condition_support/inputs/` at the repository root.
The scripts plot saved values without LR inference, primary model fitting
or multiple-testing correction. The
[condition-support analysis](../../../lung_cohort/condition_support/README.md)
describes support labels separately from condition-association statistics.
