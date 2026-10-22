# Revised Figure 2

Only Figure 2 is updated. Original package figures and manuscript files remain intact.

- `figures/Figure2.pdf`: editable PDF master.
- `figures/Figure2.svg`: editable SVG, including text.
- `figures/Figure2.png`: 600 dpi preview.
- `source_data/Figure2_Source_Data.xlsx`: exact per-panel source tables and feature mapping.
- `reports/QA_REPORT.md`: verification and visual inspection notes.

Regenerate the supplied saved sources:

```sh
python3 -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/redraw.py --figure2-only
```

The three example matrices can be re-extracted from a full saved score table with:

```sh
python scripts/prepare_example_sources.py --score-table /path/to/slice_scores.tsv.gz
```

That input is a slice-by-feature gzip TSV with `sample_id` and the complete feature
keys as columns. Exact input checksums are recorded in reports. No LR inference,
primary model fitting or multiple-testing correction is performed by these scripts.

Run commands from this dated figure directory. The `preview/` directory is
provided for generated previews. The original export environment is in `logs/`;
the dependency list is shared with the current Figure 1 environment.

Supplementary Figure 4D retains I03, I21 and I19 and is a separate archived
interval display. It should not be described as showing all three updated
main Figure 2D examples (I15, I06, I19).
