# SpatialESS manuscript figures

## Current main figures

The dated directories below are the current main-figure masters. They supersede
the main figures in the original `figures/` directory; historical exports and
their sources remain archived.

| Figure | Editable masters and preview | Source data and reproduction |
| --- | --- | --- |
| Figure 1: numerical agreement, scaling and local-neighbor sparsity | [PDF](updates/Figure1_Atera_20261016/figures/Figure1.pdf), [SVG](updates/Figure1_Atera_20261016/figures/Figure1.svg), [PNG](updates/Figure1_Atera_20261016/figures/Figure1.png) | [Workbook](updates/Figure1_Atera_20261016/source_data/Figure1_Source_Data.xlsx), [panel tables](updates/Figure1_Atera_20261016/source_data/panels/), [code and QA](updates/Figure1_Atera_20261016/README.md) |
| Figure 2: patient-aware pulmonary-fibrosis score analysis | [PDF](updates/Figure2_20261012/figures/Figure2.pdf), [SVG](updates/Figure2_20261012/figures/Figure2.svg), [PNG](updates/Figure2_20261012/figures/Figure2.png) | [Workbook](updates/Figure2_20261012/source_data/Figure2_Source_Data.xlsx), [panel tables](updates/Figure2_20261012/source_data/panels/), [code and QA](updates/Figure2_20261012/README.md) |

Figure 1A uses the author-supplied vector asset, placed proportionally without
changing its internal content. Figure 1F-H retain the measured radius-graph
analysis. Atera is a fixed-LR benchmark on preview inputs; inference and failure
scopes are documented in the experiment notes and source data.

Figure 2C retains all 28 selected features and separates its support, patient ICC
and singular annotation headers from the data. Figure 2D uses I15, I06 and I19,
with nonzero-slice counts for each condition. All original model estimates,
q values and stable feature keys are retained.

## Supplementary figures

The existing [Figure S1](figures/FigureS1.pdf), [Figure S2](figures/FigureS2.pdf),
[Figure S3](figures/FigureS3.pdf) and [Figure S4](figures/FigureS4.pdf) remain
current. Matching SVG/PNG files, per-panel CSV files, original code and the
original `source_data/Source_Data.xlsx` workbook are retained in this directory.
Definitions and labels are documented with the source data and verification reports.

S4D retains I03, I21 and I19. It does not show intervals for all three revised
Figure 2D examples; use the current main Figure 2 labels and the saved model
tables to avoid mixing their identities. No supplementary data or statistics
were altered for the GitHub update.

## Redraw saved data

Follow the commands in each dated main-figure README. For Figure 2:

```sh
cd paper/figure_package/updates/Figure2_20261012
python3 -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 scripts/redraw.py --figure2-only
```

For Figure 1:

```sh
cd paper/figure_package/updates/Figure1_Atera_20261016
python3 -m pip install -r requirements.txt
python3 scripts/redraw.py --figure1-only
python3 scripts/compose_figure1.py
```

These commands redraw frozen values, not inference or model fits. Rendering
requires the dependencies and system PDF utilities documented in the source
scripts. Arial is used where installed, with a recorded sans-serif fallback.
Fonts are not redistributed. Dense point layers may be rasterized; labels,
axes, legends and the other data primitives remain vector elements.

The original `scripts/run_all.py` is the historical S1-S4/all-figure pipeline.
It does not apply the current dated Figure 1/Figure 2 revisions. Its staging
steps that refer to original server paths are archival provenance.

The [manuscript materials index](../../docs/MANUSCRIPT_MATERIALS_20261022.md)
links experiment methods, complete model results and verification commands.
