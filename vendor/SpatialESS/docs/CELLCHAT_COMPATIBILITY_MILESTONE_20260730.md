# SpatialESS CellChat-compatibility milestone

Date: 2026-07-30

## Method boundary

This implementation is an independent SpatialESS compatibility path. It uses
the sparse-neighborhood and structure-preserving-null concepts studied from the
kidney interaction framework, CytoSignal and FastCCC, but does not copy or
relabel any of those methods. The implemented observed-score order is:

1. exact CellChat type-7 group triMean on referenced genes only;
2. geometric-mean ligand and receptor complex expression;
3. activating and inhibitory coreceptor adjustment;
4. CellChat Hill activation;
5. sender and receiver agonist/antagonist factors;
6. evaluation only on group pairs supported by the cell-level spatial graph.

The raw CellChat molecular probability and spatially weighted probabilities
are stored separately. This milestone does not yet include calibrated spatial
permutation p-values or an analytic/hybrid spatial null.

## Numerical validation

The deterministic small-data validation directly calls CellChat's
`triMean()`, `computeExpr_LR()`, `computeExpr_coreceptor()`,
`computeExpr_agonist()` and `computeExpr_antagonist()` helpers. The sparse
triMean difference was zero and the maximum complex/cofactor probability
difference was `2.78e-17`.

All 67 package assertions pass. The full Xenium run was independently checked
against the same CellChat helper formulas on every emitted LR-group record.
The maximum difference was `1.11e-16`.

## Full Xenium result

Canonical input:

- 1,156,091 cells;
- 648 signaling genes;
- 534 total referenced genes: 504 ligand/receptor genes plus 30 genes used
  only by agonist, antagonist or coreceptor definitions;
- 627 LR records, including 220 complex LR records;
- 2,357 spatial grid groups;
- 24,417,532 directed cell-level spatial edges.

Sparse output:

- 18,597 spatially supported group pairs;
- 468,314 positive LR-group-pair records;
- 26 checkpointed LR chunks;
- maximum chunk object size 12.43 MiB.

Cold-stage timings:

- referenced-gene sparse triMean: 6.730 s;
- cell-edge to group-support collapse: 0.313 s;
- CellChat-helper validation reference: 4.480 s;
- SpatialESS 627-LR core scoring: 0.428 s;
- measured end-to-end checkpoint-validation process: 33.798 s;
- peak RSS: 4,561,848 kB (4.35 GiB).

The end-to-end timing includes object and graph loading, cold preprocessing,
construction of the independent CellChat-helper reference, checkpoint loading,
full numerical validation and output bookkeeping. Core scoring is reported
separately because the end-to-end process deliberately performs extra
validation work that is not required for routine inference.

## Reproducibility

- Full script:
  `benchmarks/run_cellchat_compat_all627_prepared1156091.R`
- Small direct CellChat validation:
  `benchmarks/validate_cellchat_compat_small.R`
- Full results:
  `benchmarks/results/cellchat_compat_all627_prepared1156091`
- Primary summary: `summary.tsv`
- Peak-memory log: `run.log`
- Initial all-chunk scoring log: `scoring_run.log`

## Next statistical gate

The next required method-development stage is a sample-, compartment- and
spatial-block-preserving permutation engine. Analytic or hybrid significance
must not be claimed until type-I error and FDR are calibrated against a
high-permutation reference. A biologically annotated spatial dataset and true
segmentation-contact polygons are also still required.
