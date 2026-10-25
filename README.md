# SpatialESS

[![R-CMD-check](https://github.com/Blake-Deng/SpatialESS/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/Blake-Deng/SpatialESS/actions/workflows/R-CMD-check.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23591306.svg)](https://doi.org/10.5281/zenodo.23591306)

SpatialESS is an R package for spatial cell–cell communication inference. It
uses sparse radius graphs and streamed ligand–receptor scoring and permutation
tests, with numerical validation against the tested SpatialCellChat V3 workflow.
For cohorts with repeated tissue sections, it also provides patient-aware
mixed-effects analysis.

## Get started

Requires **R ≥ 4.1** and a **C++17 compiler**.

```sh
git clone https://github.com/Blake-Deng/SpatialESS.git
cd SpatialESS
Rscript -e 'install.packages(c("Rcpp", "Matrix"))'
R CMD INSTALL .
Rscript examples/lr_panel_example.R
```

The example uses synthetic data and needs no external dataset or CellChat
installation. For cohort modelling, also install `lme4`.

To analyze your own data, follow the [analysis guide](docs/QUICKSTART.md) and
[input requirements](docs/INPUTS.md). Expression matrices are gene-by-cell;
coordinates and cell-type labels must follow the same cell order. Required
ligand and receptor subunits must be measured; excluded interactions are
reported in an explicit coverage audit.

## Patient-aware analysis

SpatialESS retains individual sections and accounts for repeated sampling:

```text
log1p(score) ~ condition + (1 | patient)
```

Estimates are reported on the log1p-score scale. BH correction includes fits
with `status == "ok"`, including valid singular fits. See the
[model definitions and diagnostics](docs/MULTISAMPLE_LMM.md).

## Research materials

| Material | Start here |
| --- | --- |
| Figures, source data and experiment overview | [Materials index](docs/SUBMISSION_CONTENTS.md) |
| Data download, verification and figure reproduction | [Reproduction guide](docs/REPRODUCE.md) |
| Inference, numerical agreement and resource measurements | [Result definitions](docs/RESULTS.md) |
| Pulmonary-fibrosis cohort: 343 genes, 37 LR records | [Primary cohort](paper/lung_cohort/primary_343gene37LR/README.md) |
| Observed condition support and section-selection sensitivity | [Condition-support analysis](paper/lung_cohort/condition_support/README.md) |
| Atera validation and component ablations | [Atera](paper/atera_validation/README.md) · [Ablations](paper/ablation/README.md) |
| Demonstration and measured setup requirements | [Demo](SpatialESS-demo/README.md) |

The repository contains software, analysis scripts, documentation and the
self-contained condition-support tables. Manuscript figures and larger saved
results are supplied through the companion data package described in the
reproduction guide. Raw expression matrices and manuscript Word drafts are
not included.

Observed support and condition association answer different questions. A
feature shared across conditions can still have a significant condition
association; a zero score does not establish biological absence. Numerical
agreement with a reference workflow does not establish biological validity.

## Citation and license

Citation metadata are in [CITATION.cff](CITATION.cff). SpatialESS is distributed
under **GPL-3**. SPARKLE is separate upstream correction software and should be
cited separately when used with SpatialESS.
