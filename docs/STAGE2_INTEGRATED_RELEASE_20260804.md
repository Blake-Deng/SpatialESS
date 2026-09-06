# Stage 2 Integrated SpatialESS release gate

Date: 2026-08-04

## Decision

Stage 2 meets the minimum publishable scope defined in
`MultiSpatialESS_implementation_roadmap.docx`. This decision applies to the
single-sample Integrated SpatialESS release. It does not include Stage 3
multi-sample inference or Stage 4 morphology-aware modeling.

## Integrated method

The exported `spatialess()` workflow combines:

1. separate contact and diffusion CSR radius graphs;
2. projection to LR-, complex- and cofactor-referenced genes;
3. exact sparse R type-7 triMean summaries;
4. expression-supported sender-receiver group pairs;
5. CellChat-compatible complex, cofactor and Hill scoring; and
6. fixed spatial block permutation or guarded hybrid inference.

Production inference exposes only `permutation` and `guarded_hybrid`. The
analytic component of guarded hybrid can exclude clearly non-significant
records, but every reported discovery must be confirmed by fixed spatial
permutation. Direct analytic and adaptive p-values remain experimental.

## SpatialCellChat V3 fidelity and scaling

The primary comparator is SpatialCellChat V3 because the
`spatial_v3_compatible()` path preserves its probability and permutation
definitions. The comparison used nested cells from the CosMx CancerousLiver
sample, all 373 computable LR interactions, 20 permutations and one thread.

From 250 through 50,000 cells:

- active-record Jaccard was 1 at every scale;
- probability Pearson correlation was 1;
- maximum absolute probability difference was at most `7.81e-18`;
- p-values were identical;
- significant-record Jaccard was 1; and
- no result crossed the 0.05 threshold differently.

At 50,000 cells, SpatialESS completed in 10.48 s with 0.802 GiB peak RSS,
whereas SpatialCellChat V3 required 23,919 s with 60.78 GiB peak RSS.
SpatialCellChat V3 terminated by signal 11 at 100,000 cells after 21,655 s and
86.53 GiB peak RSS. At 200,000 and 443,515 cells it failed while requesting
298.0 GiB and 1,465.6 GiB vectors, respectively. The compatible SpatialESS path
completed 443,515 cells in 101.75 s with 1.97 GiB peak RSS.

These results establish formula fidelity for the compatibility path. They must
not be presented as equivalence for the separate integrated group-first score.

## Integrated real-data application

The Integrated API was run on the QC-filtered CosMx CancerousLiver sample:

- 443,515 cells, 419 signaling genes and 12 cell groups;
- 373 LR interactions, including 86 contact and 287 diffusion interactions;
- contact radius 15 um, diffusion radius 35 um and 100 um spatial blocks;
- 100 fixed spatial permutations; and
- one thread.

The complete process required 46.60 s and 1.953 GiB peak RSS. It built
1,456,470 contact edges and 8,619,362 diffusion edges, emitted 651 active
records and identified 435 records at `p <= 0.05`. All 435 were
permutation-confirmed.

Dominant significant pathways included COLLAGEN, MIF, SPP1, MHC-II, FN1, APP,
VTN, GALECTIN and MHC-I. Examples requiring biological follow-up include
tumor-boundary/core to macrophage or T-cell `MIF_CD74_CXCR4`, tumor to
macrophage `APP_CD74`, cholangiocyte to tumor-core
`SPP1_ITGAV_ITGB1`, and tumor to gamma-delta T-cell `HLA-A_CD8A`.
These are computationally supported hypotheses, not causal validation.

## Million-cell gate

The Integrated API completed the 10x Xenium Prime fresh-frozen human ovary
stress test with:

- 1,156,091 cells, 648 signaling genes and 2,357 spatial grid groups;
- all 627 LR interactions, including 129 contact and 498 non-contact LR;
- contact radius 15 um and diffusion radius 35 um;
- 17,966 spatial blocks of 100 um;
- 20 fixed spatial permutations; and
- one thread.

It built 5,781,482 contact and 33,297,328 diffusion directed edges. The run
emitted 430,406 active records, of which 57,026 had `p <= 0.05`; no significant
record lacked permutation confirmation. Internal end-to-end time was 75.52 s.
Complete-process wall time was 79.78 s and peak RSS was 5.61 GiB.

With 20 permutations, the minimum finite-corrected p-value is `1/21`, so this
run is a computational stress test rather than a final high-resolution
biological analysis.

## Comparator boundary

The formal comparator is SpatialCellChat V3 because the compatible path
preserves its probability and permutation definitions. Independent spatial CCC
methods use different LR databases, neighborhood graphs, score definitions and
null models and are therefore discussed as related work rather than treated as
numerical ground truth in this release.

## Additional release validation

Five independent 5,000-cell runs with identical settings were exactly
reproducible. Method-time, process-wall and peak-RSS coefficients of variation
were 1.78%, 4.91% and 0.28%, respectively. Across five distinct
100-permutation seeds, all runs returned the same 148 significant records and
p-value rank correlation was at least 0.999864. A 20-200 permutation
sensitivity series retained exactly identical observed probabilities; the
100-permutation significance set had Jaccard 0.993 with the 200-permutation
reference and one threshold disagreement. Full results and interpretation
boundaries are in `docs/PUBLICATION_VALIDATION_20260804.md`.

## Evidence

- V3 resources: `/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead/ladder_method_results.tsv`
- V3 accuracy: `/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead/ladder_accuracy_results.tsv`
- Integrated CosMx: `/data/dzf/SecAct_2026/results/spatialess_integrated_stage2/stage2_integrated_results.tsv`
- CosMx biology: `/data/dzf/SecAct_2026/results/spatialess_integrated_stage2/cells443515_full_permutation_perm00100/`
- Million-cell Xenium: `/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/spatialess_integrated_stage2/cells1156091_grid50_permutation_perm00020/`
- Publication validation: `/data/dzf/SecAct_2026/results/spatialess_publication_validation/`

## Remaining work outside Stage 2

- segmentation-polygon contact graphs when boundary data are available;
- external biological validation beyond computational pathway plausibility;
- multi-sample, sample-as-replicate inference in Stage 3; and
- morphology-aware modeling in Stage 4.
