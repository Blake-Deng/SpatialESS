# SpatialESS cohort sensitivity and robustness-layer prototype

Saved-score analyses of 45 GSE250346 slices from 35 patients; no LR inference,
original manuscript, original benchmark or primary model was overwritten.

Start with `reports/SENSITIVITY_RESULTS.md` and `reports/ALGORITHM_DESIGN.md`.
The primary engine remains 0.1.4. `figures/FigureS5` and
`data/Supplementary_Table_2_Cohort_Sensitivity.xlsx` contain the added evidence.

## Reproduce from the supplied score matrix

Requires R (lme4, data.table) and Python (numpy, scipy, pandas, matplotlib,
statsmodels, openpyxl, tabulate). Arial must be installed for identical typography.
Exact versions are in `logs/R_session_info.txt` and `logs/python_environment.json`.
Linux/macOS fork workers are used for mixed-model fits; keep native/BLAS threads
at one. Choose a new output path; refitting deliberately refuses an existing one.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript scripts/refit_from_prepared.R . /tmp/spatialess_sensitivity_reproduced 16
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/run_patient_sensitivity.py --out /tmp/spatialess_sensitivity_reproduced --figure-source data
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript scripts/verify_robust_layer.R /tmp/spatialess_sensitivity_reproduced
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/summarize_and_plot.py --out /tmp/spatialess_sensitivity_reproduced
```

The original extraction runner `scripts/run_adjusted_models.R` additionally accepts
the repository and directory of 45 saved original slice profiles, with signatures
and resumable checkpoints. Their checksums and the exact model source snapshot are
supplied. Full input profiles and raw Seurat data are not required to reproduce the
new fits from the included prepared matrix.

Raw full-universe tables are gzip TSV. The workbook presents summaries, complete
identities of the original 28 plotted features, zero-inclusive support, influence,
fit states, design ranks and all-feature descriptive agreement labels. `q` always
refers to the relevant model's full eligible testing universe. Balanced sampling
ranges are descriptive, not confidence intervals.
