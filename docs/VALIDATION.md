# Validation and interpretation

## Questions tested

The release records two distinct comparisons.

1. **Implementation fidelity:** Does SpatialESS reproduce the tested official
   SpatialCellChat V3 output when both receive the same SPARKLE-corrected
   expression, labels, coordinates, LR database, parameters and seed?
2. **Input sensitivity:** How does applying SPARKLE before SpatialESS alter the
   inferred communication set relative to RAW expression?

Only the first comparison is an exactness test. The second comparison measures
the consequence of changing the input expression.

## Ovarian validation design

- platform: 10x Genomics Visium HD human ovarian cancer;
- analyzed cells: 16,247;
- shared RAW-derived RCTD labels: 15 groups;
- signaling genes after matching: 967;
- computable ligand-receptor interactions: 1,828;
- permutations: 100;
- threads: 1;
- identical coordinates, labels, LR order, seed and inference parameters.

A tissue-matched 10x scFFPE ovarian reference and official barcode annotation
were used for RCTD/label transfer in this benchmark. These files are not runtime
dependencies and are not redistributed here.

## Official V3 versus SpatialESS

| Metric | Official V3 | SpatialESS |
|---|---:|---:|
| Core inference time | 6,541.848 s | 35.765 s |
| Independent process wall time | 6,549.0 s | 40.6 s |
| Peak RSS | 4.008 GiB | 0.590 GiB |
| Active records | 37,808 | 37,808 |
| Significant records | 8,417 | 8,417 |

Derived results:

- core-computation speedup: 182.9-fold;
- process-wall speedup: 161.3-fold;
- peak RSS reduction: 85.3%;
- active-record Jaccard: 1.000;
- maximum absolute probability difference: 4.987e-18;
- bitwise-exact probability fraction: 0.7848;
- maximum absolute p-value difference: 0;
- exact p-value fraction: 1.000;
- significance disagreements at p < 0.05: 0.

The probability values are numerically equivalent to floating-point precision;
not every value is bitwise identical because the implementations use different
storage and traversal orders.

## RAW versus SPARKLE-corrected input

Both inputs were analyzed by SpatialESS with the same cells, labels,
coordinates, LR set, parameters and seed.

| Metric | RAW | SPARKLE |
|---|---:|---:|
| Active records | 41,294 | 37,808 |
| Significant records | 10,942 | 8,417 |
| Mean scFFPE reference correlation | 0.4748 | 0.5301 |

Additional comparisons:

- retained input counts: 51.1%;
- aggregate LR Spearman rho: 0.9864;
- top-20 LR overlap: 16/20;
- top-50 LR overlap: 41/50;
- top-100 LR overlap: 85/100;
- shared significant records: 7,227;
- lost after SPARKLE: 3,715;
- gained after SPARKLE: 1,190;
- significant-record Jaccard: 0.5957.

These changes are expected because SPARKLE modifies the expression input. They
must not be described as an error between inference implementations.

## Claim boundary

The official comparison supports **numerical fidelity** for the tested
SpatialCellChat V3-compatible workflow. It does not establish universal
biological accuracy.

The RAW-versus-SPARKLE comparison demonstrates integration behavior and input
sensitivity. A stronger claim that SPARKLE improves biological truth would
require independent truth or controls across additional samples, such as
simulation, orthogonal single-cell references, marker-retention checks,
parameter sensitivity and experimental LR evidence.

## Reproducibility assets

- `source_data/official_vs_spatialess.tsv`: exactness summary.
- `source_data/official_vs_spatialess_records.tsv.gz`: matched active records.
- `source_data/resource_summary.tsv`: process and core resource measurements.
- `source_data/comparison_metrics.tsv`: RAW-versus-SPARKLE metrics.
- `source_data/lr_aggregate_comparison.tsv`: LR-level aggregate scores.
- `source_data/significant_transition.tsv`: shared, lost and gained records.
- `source_data/expression_reference_summary.csv`: reference comparison.
- `scripts/export_figure_data.R`: compact source-data exporter.
- `../../scripts/make_publication_figures.py`: figure generator.

Large expression matrices, raw 10x files and downloaded annotation references
are deliberately excluded from the repository.
