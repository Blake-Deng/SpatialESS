# Results

For tables and saved records, download and mount the companion data package
as described in the [reproduction guide](REPRODUCE.md). The paths below are
logical paths after mounting, rather than directly browsable GitHub files.
The self-contained [condition-support analysis](../paper/lung_cohort/condition_support/README.md)
includes its own inputs, classifications and validation records.

The inference and resource measurements below used **SpatialESS 0.1.4** on
2026-09-30 through 2026-10-01, with streaming aggregation and the historical
V3 percentage rule. Each configuration has its own input checksum, LR table
and parameter contract. The series predates both the default reference
aggregation in 0.1.4.9003 and the default numeric filter in 0.1.4.9004.

| Table | Contents |
|---|---|
| inference_summary.tsv | Saved records, cells, genes, LR, permutations, spatial edges and raw-p counts |
| inference_resources.tsv | Each completed single-fixture process, inference time and peak RSS |
| parameters.tsv / inputs.tsv | Inference settings and prepared-input SHA256 |
| official_fidelity.tsv | Full-precision comparisons to the identified official references |
| cosmx_scaling.tsv | SpatialESS 0.1.4 original-harness medians and source-labeled official measurements |
| cosmx_resource_repetitions.tsv | Individual SpatialESS 0.1.4 harness measurements |
| spatial_sparsity.tsv | Radius-graph edges, degree, density and storage definitions |
| bridge_results.tsv | RAW/SPARKLE H5AD pipeline measurements |
| gse306130_sampling.tsv | Fixed-parent cluster provenance and preparation seeds |
| gse250346_lmm_summary.tsv | Completed primary cohort model summary |
| lmm_provenance.tsv | Model source, fitting date and cohort input scope |

All table paths below refer to the mounted `results/tables/` directory. `results/records/<case>.rds` stores each saved result. `results/contracts/<case>.rds` stores its LR, parameters and diagnostics.

## Dataset coverage

| Configuration family | Inputs |
|---|---:|
| cervical | 8 |
| cosmx | 9 |
| gse306130 | 8 |
| gse313006 | 2 |
| lung_paper | 45 |
| lung_signaling | 45 |
| lung_teacher | 45 |
| ovarian | 6 |

The lung panel/cohort input contains 343 measured genes; the official-comparator configuration contains 99 signaling genes. The `lung_teacher` input uses raw counts and the recorded 250-micrometre interaction range. These configurations are not pooled.

## Measurement interpretation

Completed official comparisons report all measured metrics and the strict fidelity flag. Official jobs without a completed result do not contribute a fidelity measurement. Dimensions-only references have a distinct method label.

The CosMx table uses medians across three ESS repetitions at 50k and 443,515 cells; other displayed harness scales have one observation. The server was shared and native threads were limited to one. Wall time includes loading and writing; single-fixture worker RSS and full-parent-fixture RSS have different scopes.

LMM tables represent the completed September fit using the included model implementation; see the fitting date in lmm_provenance.tsv. Convergence-warning fits do not receive BH significance calls.

Significance counts in single-sample tables use raw permutation p-values. Patient-level q-values use the documented BH adjustment. Sparse graph observations describe the tested inputs and radius settings, not universal linear algorithmic complexity or biological correctness.
