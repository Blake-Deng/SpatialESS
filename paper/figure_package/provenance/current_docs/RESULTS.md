# Results

SpatialESS version: **0.1.4**. Inference and resource measurement series: 2026-09-30 through 2026-10-01. Each configuration has its own input checksum, LR table and parameter contract.

| Table | Contents |
|---|---|
| inference_summary.tsv | Current records, cells, genes, LR, permutations, spatial edges and raw-p counts |
| inference_resources.tsv | Each completed single-fixture process, inference time and peak RSS |
| parameters.tsv / inputs.tsv | Inference settings and prepared-input SHA256 |
| official_fidelity.tsv | Full-precision comparisons to the identified official references |
| cosmx_scaling.tsv | Current ESS original-harness medians and source-labeled official measurements |
| cosmx_resource_repetitions.tsv | Individual current ESS harness measurements |
| spatial_sparsity.tsv | Radius-graph edges, degree, density and storage definitions |
| bridge_results.tsv | RAW/SPARKLE H5AD pipeline measurements |
| gse306130_sampling.tsv | Fixed-parent cluster provenance and preparation seeds |
| gse250346_lmm_summary.tsv | Completed primary cohort model summary |
| lmm_provenance.tsv | Model source, fitting date and cohort input scope |

All table paths are relative to `results/tables/`. `results/records/<case>.rds` stores each current result; `results/contracts/<case>.rds` stores its LR, parameters and diagnostics.

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
