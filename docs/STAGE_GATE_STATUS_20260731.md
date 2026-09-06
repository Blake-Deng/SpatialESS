# SpatialESS roadmap gate status

Last updated: 2026-08-04

This file evaluates the current implementation against
`MultiSpatialESS_implementation_roadmap.docx`. A working prototype is not
reported as a completed release gate when required calibration or comparator
experiments remain missing.

| Stage | Status | Completed evidence | Remaining release gate |
|---|---|---|---|
| 0, frozen CellChatESS | Substantially complete | Frozen v0.2.1 baseline; CellChat, dense Rcpp, AutoZyme and ESS comparisons; 216-object and million-cell records | Preserve release artifacts and rerun only revision-required analyses |
| 1A, sparse exact triMean | Complete | Exact sparse zero-aware R type-7 Q1/Q2/Q3 kernel; shared multi-select production path; 30-condition runtime matrix; million-cell exact regression; boundary tests | Preserve regression and benchmark records |
| 1B, analytic/hybrid significance | Complete for guarded fixed-permutation mode | One-way analytic exclusion; every discovery permutation-confirmed; 100 x 499 and 20 x 9,999 references; 108-scenario matrix; 28-scenario zero-boundary scan; record-sliced C++ fallback | Keep adaptive p-values out of BH/FDR; treat larger LR panels as sensitivity analysis |
| 1C, edge-first SpatialESS | Complete for radius-graph release | CSR radius/kNN graphs; sparse support; CellChat complexes/cofactors/regulators/Hill score; V3-equivalent compatibility path through 50,000 cells; million-cell CSR execution | Segmentation-polygon contact graph remains an optional extension when boundaries are available |
| 2, integrated SpatialESS | Complete for roadmap minimum scope | Exported integrated API; calibrated production inference boundary; CosMx application; SpatialCellChat V3 fidelity/scaling; 1,156,091-cell Xenium stress test | External biological replication; Stage 3 sample-as-replicate inference |

## Current verified implementation

- The package compiles with R 4.5.3 and GCC 15.2.0.
- All current test assertions pass, including sparse triMean boundaries, homogeneous and stratified analytic guards, and hybrid fallback concordance.
- The full Xenium compatibility run includes 1,156,091 cells, 2,357 groups,
  627 LR records and 220 complex LR records.
- Molecular probabilities agree with the independent CellChat helper reference
  to a maximum absolute difference of `1.11e-16`.
- The Integrated API completed 1,156,091 Xenium cells, 627 LR interactions
  and 20 spatial permutations in 79.78 s complete-process wall time with
  5.61 GiB peak RSS. All 57,026 significant records were permutation-confirmed.
- On the CosMx ladder, the V3-compatible path matched SpatialCellChat V3
  through 50,000 cells with maximum probability difference `7.81e-18`,
  identical p-values and zero significance-threshold disagreements.
- The 443,515-cell CosMx biological run used 100 spatial permutations and
  completed in 46.60 s with 1.95 GiB peak RSS.
- The earlier full CSR compatibility workflow remained below 5 GiB peak RSS.
- The SecAct CosMx adapter passes a complete synthetic RDA integration test,
  including peak-RSS logging and compressed significance output.
- The narrow homogeneous simple-LR analytic candidate achieved type-I 0.040
  versus 0.035 for 999-permutation inference across 200 null replicates, with
  99.5% threshold agreement. A 20-replicate 9,999-permutation reference gave
  type-I 0.050 for both methods and 100% threshold agreement.

## Stage 1B statistical boundary

FastCCC's FFT framework is a design reference, not a drop-in null model for
CellChat triMean. CellChat averages Q1, Q2 and Q3 computed from the same sample,
so these order statistics are correlated. Convolving three marginal quantile
distributions as if independent would generally be miscalibrated. Spatial
blocks also create stratified, non-identically distributed sampling. The
first two-sided hybrid failed the high-permutation gate and the first
one-way expanded matrix exposed false-negative screening under zero inflation.
The final guard forces records above 10% weighted zero mass to fixed
permutation. Across the corrected 108-scenario matrix and 28-scenario
zero-boundary scan, every hybrid/full-permutation rejection set was identical
and every discovery was permutation-confirmed. Direct analytic p-values and
adaptive BH/FDR remain experimental.

## Next-stage sequence

1. Freeze the Stage 2 API, benchmark inputs and statistical boundary.
2. Add external biological replication without changing the score definition.
3. Begin Stage 3 with biological samples, not cells, as inferential replicates.
4. Keep segmentation-polygon and morphology-aware models as separately
   validated extensions.

The full Stage 2 decision and evidence paths are recorded in
[STAGE2_INTEGRATED_RELEASE_20260804.md](STAGE2_INTEGRATED_RELEASE_20260804.md).
