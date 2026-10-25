# Research materials

SpatialESS combines an installable R package with scripts, parameter contracts,
validation records and saved analyses. Start with the [analysis guide](QUICKSTART.md)
for your own data, or the [reproduction guide](REPRODUCE.md) for manuscript results.

The **Location** column links to code and documentation in this repository.
Most figure exports, workbooks and large result tables are in the companion
data package; mount that package before using their documented file paths.
The condition-support analysis includes its own inputs and results.

| Material | Location | Scope |
| --- | --- | --- |
| Figure 1 | [Figure package](../paper/figure_package/README.md) | Numerical agreement, resource measurements, spatial sparsity and Atera scaling |
| Figure 2 | [Figure package](../paper/figure_package/README.md) | Pulmonary-fibrosis score analysis; panel D uses I15, I06 and I19 |
| Supplementary figures | [Figure package](../paper/figure_package/README.md) | S1–S4; S4D uses I03, I21 and I19 |
| Primary pulmonary-fibrosis cohort | [343-gene cohort](../paper/lung_cohort/primary_343gene37LR/README.md) | 45 sections, 35 patients, 50,040 feature fits and zero-inclusive support |
| Condition-support classification | [Condition support](../paper/lung_cohort/condition_support/README.md) | Post hoc descriptive classes, 40%/50%/60% thresholds and all 768 one-section-per-patient selections |
| Cohort model sensitivity | [Sensitivity analyses](../paper/cohort_sensitivity/README.md) | Patient/run/TMA analyses, patient omission and balanced selections |
| Official lung comparison | [99-gene reference validation](../paper/lung_cohort/reference_99gene_validation/README.md) | Numerical fidelity; separate from the 343-gene cohort model |
| Atera validation | [Atera](../paper/atera_validation/README.md) | Preparation, fixed-LR comparisons, sampled IDs, measurements and failed attempts |
| Component ablations | [Ablations](../paper/ablation/README.md) | Variant definitions, 170 technical runs and numerical verification |
| TSTA | [Feasibility evidence](../paper/tsta/README.md) | 1,258 sample summaries; SpatialESS feasibility rather than paired reference validation |
| Supplementary Table 1 | [Table definitions](../paper/tables/README.md) | Dataset and benchmark summary |
| Annotation and record verification | [Ovarian provenance](../paper/provenance/ovarian/README.md) | Label-transfer evidence, source identities and environment limits |

## Reading the cohort results

The primary cohort contains 39,570 valid fits, including 17,308 singular fits,
and 5,147 condition associations at primary BH q < 0.05. These fits use the
saved September 2026 model; the package version does not imply a refit.
The 99-gene reference comparisons use a distinct preparation.

Condition-support labels describe how often saved nonzero scores occur in
each condition. They retain all 50,040 features regardless of fit status and
do not change the primary estimates, P values or q values. Shared support is
compatible with a significant difference in score. Zero does not establish
biological absence.

Figure 2D and S4D display different example sets. Use their full feature keys
when matching plots to model or support tables. A beta divided by response SD
is a descriptive standardized effect, not a fold change.

## Verification

From the repository root:

```sh
python3 tools/verify_software.py
```

This checks the current software manifest. The companion data, fixed figure
inputs and condition-support outputs have separate checks described in the
[reproduction guide](REPRODUCE.md).
