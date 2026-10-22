# Atera WTA supplementary validation

Isolated SpatialESS 0.1.4 experiments on newly downloaded breast and cervical
Atera WTA preview specimens. Original data, manuscript and package are not edited.

## Queue

Run `scripts/run_queue.py` with the Python executable in `config/design.json`.
The queue runs one process on one available CPU with single-thread BLAS settings,
a 96 GiB address-space cap, no core dumps, and per-task timeouts. A filesystem
lock prevents duplicate queues. Successful checkpoints resume; failures remain
recorded and are not automatically retried. Use `scripts/status.py` for live PID
and log checks. Background session: `spatialess_atera_20261013`.

Preparation is followed by restricted-LR 250-cell paired smoke tests, full-LR
1,000-cell pilots with 20 permutations, formal ESS central subsets (1k, 5k, 25k,
50k, 100k) and the full eligible tissue with 100 permutations, then three
independent proportional-group 50k samples and three deterministic 50k blocks.
Reference 1k/5k formal jobs are bounded to one hour each. Smoke probability/key
disagreement stops dependent tissue experiments; strict p-value discrepancies
are retained explicitly. Two tissues reach preflight before large jobs.

## Inputs and definitions

Only 18,028 measured Gene Expression rows are normalized as
`log1p(10000 * count / measured-GEX library size)`. Cells require >=20 measured
counts, >=20 measured genes, and an unambiguous official annotation. All
exclusions, original indices, barcodes, group quotas and seeds are saved.
Central selection takes nearest centroids to the mean eligible coordinate;
block selection uses nearest observed anchors at x quartiles and median y.
Stratified selection uses largest-remainder quotas and R seeds 20261013–15.

LR rows are selected by complete mandatory measured ligand/receptor component
coverage only; absent optional cofactors remain audited. No outcome-driven
selection or overexpression LR screening is performed. SpatialCellChat receives
the same complete expression panel in `data.signaling` to preserve normalization
scope. These are fixed-LR engine benchmarks, not an unchanged standard
overexpression-selection workflow. Interaction/contact/tolerance are 30/10/5
micrometres. Per-tissue distance scaling is fixed at twice the reciprocal of
the smallest eligible nonzero centroid distance and shared between engines.

Whole-process wall time and peak RSS are recorded by GNU time. Inference-only
time is separately measured in R. Directed CSR graph edges exclude graph self
edges, while the established self-score convention remains unchanged. Array
storage and theoretical dense numeric storage do not represent measured RSS.
Reference failures contribute neither successful runtimes nor agreement.

## Outputs and limitations

`tables/Atera_supplementary_results.xlsx` and TSV/CSV tables refresh after every
run. Full feature-key agreement, probability errors, p-value exceptions,
coverage, graph invariants, descriptive sampling comparisons, source hashes,
package environments, logs and run contracts are retained. `figures/` contains
editable vector performance outputs using completed formal measurements only.

The two tissues are single specimens, not patient cohorts. Sampling repetitions
are not biological replicates. Broad gene coverage and newly inferred scores
do not establish biological discovery or causal signaling. Subset/full-tissue
comparisons change geometry/composition and estimands; full tissue is not a
biological ground truth. No algorithmic O(N) or cross-model equivalence claim is
made. The same-block cervical WTA and earlier Xenium data are not paired cells.
SecAct downloads and patient-level validation are not included in this queue.
