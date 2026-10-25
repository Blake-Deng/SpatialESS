# Manuscript figures

Code and method documentation are available in this repository. Figure exports,
workbooks and panel tables are in the companion data package. Follow the
[reproduction guide](../../docs/REPRODUCE.md) to download and mount the inputs.

| Figure | Code and definitions | Logical export path after mounting |
| --- | --- | --- |
| Figure 1: numerical agreement, scaling and spatial sparsity | [Figure 1](updates/Figure1_Atera_20261016/README.md) | `updates/Figure1_Atera_20261016/figures/Figure1.pdf` |
| Figure 2: pulmonary-fibrosis score analysis | [Figure 2](updates/Figure2_20261012/README.md) | `updates/Figure2_20261012/figures/Figure2.pdf` |
| Figures S1–S4 | [Plotting code](scripts/redraw.py) | `figures/FigureS1.pdf` through `figures/FigureS4.pdf` |

Paths in the final column are relative to this directory. Matching SVG/PNG
exports, source workbooks and per-panel CSVs are listed in the companion
manifest. Figure 1A preserves the supplied workflow panel without internal
edits. Atera comparisons use fixed LR records and preview inputs; measurements
and reference failures retain their documented execution scopes.

Figure 2C displays 28 selected features. Figure 2D uses I15, I06 and I19;
S4D uses I03, I21 and I19. Match panels to saved model tables by full feature
keys. The [condition-support analysis](../lung_cohort/condition_support/README.md)
provides descriptive support labels for these examples separately from their
condition-association statistics.

## Redraw saved values

From the repository root, after installing `requirements-reproduce.txt`:

```sh
python3 tools/redraw.py Figure1 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py Figure2 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S1 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S2 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S3 --materials ../SpatialESS_reproduction_20261023
python3 tools/redraw.py S4 --materials ../SpatialESS_reproduction_20261023
```

These commands plot saved values without inference, model fitting or FDR
recalculation. Figure 1 composition requires Poppler. Font and output
handling are described in the [reproduction guide](../../docs/REPRODUCE.md).
