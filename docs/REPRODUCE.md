# Reproduce SpatialESS analyses

The repository provides the R package, analysis and plotting scripts, parameter
contracts, source snapshots and validation records. Figure exports, workbooks
and most saved results are stored in a companion data package. The
[condition-support analysis](../paper/lung_cohort/condition_support/README.md)
is self-contained and can be run without that package or expression matrices.

## Package and synthetic example

Requires R ≥ 4.1 and a C++17 compiler. From the repository root:

```sh
Rscript -e 'install.packages(c("Rcpp", "Matrix"))'
R CMD INSTALL .
Rscript examples/lr_panel_example.R
```

Install `lme4` for cohort modelling and `testthat` to run the package tests:

```sh
Rscript -e 'install.packages(c("lme4", "testthat"))'
Rscript -e 'testthat::test_local(".")'
```

For your own data, follow the [analysis guide](QUICKSTART.md) and
[input contract](INPUTS.md).

## Download and verify saved materials

Download the [software source ZIP](https://zenodo.org/records/23591306/files/SpatialESS-materials-2026.10.23-r1.zip),
[companion reproduction tar](https://zenodo.org/records/23591306/files/SpatialESS_reproduction_20261023_r1.tar.gz)
and [download checksums](https://zenodo.org/records/23591306/files/ZENODO_SHA256SUMS_20261023_r1.txt).
Verify the downloads with `sha256sum -c` on Linux or `shasum -a 256 -c` on macOS.
Extract the software and companion beside each other. Run the following from
the software directory:

```sh
python3 tools/verify_software.py
python3 tools/reproduction.py verify --materials ../SpatialESS_reproduction_20261023
python3 -m pip install -r requirements-reproduce.txt
```

The companion contains one object per SHA-256 under `objects/`.
`docs/archive_20261023/MATERIALS.tsv` maps logical file paths to these objects.
To inspect the result paths used by the analysis documentation:

```sh
python3 tools/reproduction.py mount --materials ../SpatialESS_reproduction_20261023
```

Mounting creates relative symlinks, without duplicating data. Use an extraction
tool that preserves symlinks on Linux/macOS. Files documented under
`results/tables/`, `results/records/` and most figure/data directories become
available after mounting; GitHub links to the scripts and READMEs work directly.
Raw expression matrices and manuscript Word drafts are not distributed.

The software manifest verifies repository files; the companion manifest
verifies data objects. Fixed figure-input hashes and the condition-support
manifest verify their respective input sets. Source paths and logs document
the recorded experiments and are not portable runtime configuration.

## Redraw manuscript figures

Each command mounts only the inputs required by that figure:

```sh
python3 tools/redraw.py Figure1 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py Figure2 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S1 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S2 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S3 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S4 --materials ../SpatialESS_reproduction_20261023
```

The input sets are listed in `docs/archive_20261023/FIGURE_INPUTS.json`.
Figure 1 requires `pdftocairo` from Poppler for vector composition. Its workflow
panel is placed proportionally without internal edits and is verified by a
fixed SHA-256. Arial can be supplied through `SPATIALESS_FONT_DIR`; Liberation
Sans is the fallback. Rendering metadata and fonts can affect export bytes.

The commands plot saved values without inference, refitting or FDR
recalculation. Run logs go to `local_runs/`. New exports replace mounted output
links rather than modifying companion objects. Figure 2D uses I15/I06/I19;
S4D uses I03/I21/I19. Match examples by their full feature keys.

## Condition-support analysis

Requires Python ≥ 3.11, NumPy, pandas and openpyxl. All inputs and outputs are
included under `paper/lung_cohort/condition_support/`. From that directory:

```sh
python3 -m pip install -r requirements.txt
python3 scripts/analyze_condition_support.py
python3 scripts/validate_outputs.py
python3 scripts/write_report.py
```

These commands write to that directory's `results/` and `reports/`; use a copy
to retain the supplied output files. The analysis reuses saved scores for
50,040 features and all 35 patients. It computes descriptive support classes,
1,000 seeded section selections and exact probabilities across all 768
one-section-per-patient choices. The independent validator uses a
Poisson-binomial calculation. Primary model estimates, fit states and
statistical calls are retained without new fitting or BH correction.

Verify the supplied files before regenerating them, using
`sha256sum -c SHA256SUMS.txt` on Linux or `shasum -a 256 -c SHA256SUMS.txt` on
macOS. See the [analysis README](../paper/lung_cohort/condition_support/README.md)
for definitions and interpretation limits.

## Rerun inference on prepared data

The included runner accepts an R list of named arguments for `spatialess()`:

```sh
Rscript benchmarks/run_prepared.R /path/to/input.rds /path/to/output
```

Prepared expression fixtures must be obtained separately. After mounting,
`results/tables/inputs.tsv` records their filenames and SHA-256 values, and
`results/tables/parameters.tsv` records settings. Full inference experiments
require the documented public inputs, the relevant R environment and path
adaptation in their configurations. Scripts and recorded paths do not by
themselves establish a portable end-to-end rerun.

Whole-process resource measurements include loading and writing; internal
inference timings do not. Single-fixture RSS and full-parent-fixture RSS have
different scopes. Reference and SpatialESS runtime measurements retain their
recorded execution settings and should not be treated as a paired timing experiment.

## Supplement mapping

| Item | Logical path after mounting |
| --- | --- |
| Supplementary Data 1 | `paper/submission_materials/Supplementary_Data_1_Ablation.xlsx` |
| Supplementary Data 2 | `paper/submission_materials/Supplementary_Data_2_TSTA.xlsx` |
| Supplementary Data 3 | `paper/submission_materials/Supplementary_Data_3_Cohort_Sensitivity.xlsx` |
| Supplementary Data 4 | `paper/submission_materials/Supplementary_Data_4_Lung_Support.tsv.gz` |
| Source Data Figure 1 | `paper/submission_materials/Source_Data_Figure1.xlsx` |
| Source Data Figure 2 | `paper/submission_materials/Source_Data_Figure2.zip` |

The [materials index](SUBMISSION_CONTENTS.md) links the matching methods and
experiment records. Condition-support tables are directly available in the
repository, including `results/all_50040_feature_classification.tsv.gz` within
the analysis directory.

## Evidence and interpretation limits

The primary lung model is the saved September 2026 fit. The 343-gene cohort
analysis and 99-gene reference validation use distinct preparations. Valid
singular fits remain eligible for primary BH correction; convergence-warning,
failed and ineligible fits retain their diagnostic states.

Saved cohort sensitivity includes primary/run/TMA/HC3/Welch tables, patient
omission summaries, balanced patient selections and scores. The all-feature
balanced-effect table is unavailable; selected balanced summaries do not
replace it. A score TSV is included with the condition-support analysis.
Zero-inclusive denominators and support labels are described there; observed
support does not establish biological absence or causal signalling.

TSTA includes 1,258 sample summaries from 132 studies, 37 cancer labels and
3,648,329 retained spots. It is feasibility evidence without paired reference
validation. Of the two representative record folders, TSTA_0001 is usable;
the copied TSTA_0117 records stream is truncated and is retained as received.
Aggregate completion does not certify the integrity of every copied output.

The original RCTD sessionInfo is unavailable. An archived spacexr DESCRIPTION
identifies 2.2.1; an exact original-run version of 2.2.16 is unverified.
Hardware snapshots collected after experiments are not per-run hardware
records. Annotation, run/region confounding and score-zero mechanisms retain
the limits documented with each analysis.
