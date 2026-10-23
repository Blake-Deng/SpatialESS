# Current manuscript materials

This index supplements the SpatialESS 0.1.4 public snapshot audited at commit
`f92dff98f9019300f8fe5086d8a011eed4f089d8`. It adds saved experiment materials;
the installable engine, model definition and original numerical results are
unchanged.

| Material | Repository location | Use |
| --- | --- | --- |
| Current Figure 1 | [Figure1_Atera_20261016](../paper/figure_package/updates/Figure1_Atera_20261016/README.md) | Current main figure; original workflow preserved, graph sparsity and Atera scaling included |
| Current Figure 2 | [Figure2_20261012](../paper/figure_package/updates/Figure2_20261012/README.md) | Main panel D uses I15, I06, I19; panel C headings separated |
| Supplementary figures | [Figure package](../paper/figure_package/README.md) | Existing S1-S4; S4D retains I03, I21, I19 |
| Atera methods and scale validation | [Atera archive](../paper/atera_validation/README.md) | Preprocessing, exact sampled IDs, fixed LR/database, contracts, source, measurements and attempt ledger |
| Component ablations | [Ablation archive](../paper/ablation/README.md) | Variant definitions, five technical repetitions, validated timing/resource summaries and logs |
| Primary PF cohort | [343-gene cohort](../paper/lung_cohort/primary_343gene37LR/README.md) | All 50,040 fits, fixed feature universe, nonzero scores and zero-fill rule, metadata and original model provenance |
| Official lung comparison | [99-gene validation](../paper/lung_cohort/reference_99gene_validation/README.md) | 42 evaluable pairs, one retained p-value discrepancy, three failed reference attempts |
| Cohort sensitivity | [Existing sensitivity material](../paper/cohort_sensitivity/) | Existing saved sensitivity analyses; not rerun for this update |
| Ovarian annotation | [Environment evidence](../paper/provenance/ovarian/RCTD_ENVIRONMENT_NOTE.md) | Archived spacexr DESCRIPTION and limits of original annotation environment recovery |

The primary model has 39,570 valid fits, including 17,308 singular fits, and
5,147 BH-significant condition associations. These are the archived September
cohort results. They are distinct from the 99-gene engine-reference comparisons.
Descriptive beta/response-SD values are not fold changes or independent
biological validation. Zero scores are not proof of biological absence.

The Atera datasets are 10x WTA FFPE preview inputs. Ten completed paired scales
(1k, 5k, 10k, 25k, 50k in each tissue) met the saved strict agreement criterion.
Larger-scale reference failures are recorded under their actual memory/time
budgets. Fixed-LR inference is not a test of the default overexpression-selection
pipeline. Technical repeats and multiple scales are not biological replicates.

## What can be reproduced directly

The dated main-figure folders contain their frozen plotting inputs and scripts.
Follow each README to redraw the saved results. Regenerating the original
all-figure pipeline does not apply the dated main-figure revisions.

Full inference reruns require the public raw datasets or prepared expression
fixtures, an appropriate R environment and path adaptation in the archived
experiment configuration. Raw expression matrices are intentionally not included.
Original absolute server paths in run logs/configuration are provenance, not
claims of portable execution without setup. The Atera barcode lists, input
checksums, preprocessing scripts and exact database/parameter records document
how the tested inputs were made.

The 2026-10-20 CPU/physical-memory snapshot was collected after the experiments
and is not a per-run historical hardware record. The original RCTD sessionInfo
is unavailable; an archived DESCRIPTION identifying 2.2.1 does not substantiate
the draft's 2.2.16 or a complete original R package environment.

## Verification

From the repository root:

```sh
python3 paper/provenance/github_update_20261022/verify_addendum.py
python3 paper/verify_manifest.py
sha256sum --check --quiet SHA256SUMS.txt
```

The 2026-10-23 candidate adds saved TSTA feasibility summaries and representative outputs under `paper/tsta/`, and the saved cohort-sensitivity analysis under `paper/cohort_sensitivity/`. The companion package must be mounted before inspecting result paths. See [the current distribution guide](REPRODUCE_20261023.md) for scope, missing items and the exact supplement mapping. No manuscript Word drafts are distributed.
