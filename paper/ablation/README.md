# SpatialESS component ablation, 2026-10-12

A component experiment using SpatialESS 0.1.4. The installable package and
manuscript source are separate from these experimental variants. All variants retain the object
protection fixes, molecular formula, group filters, self-edge convention and
permutation tail criterion. No new biological model or significance filtering
is introduced.

Saved outputs are in the companion data package; follow the
[reproduction guide](../../docs/REPRODUCE.md) to download and mount it.

## Inputs and runs

Inputs are byte-matched to the submission input manifest, copied from the saved
prepared fixtures. Core inputs: CosMx CancerousLiver 5k/10k/25k and cervical
Xenium central/cluster-stratified 5k (default filters). CosMx keeps 419 genes,
373 LR, 20 permutations and seed 20260803; Xenium keeps its saved 730 genes,
1076 LR, 100 permutations, seed 1 and default filtering settings.
The full 50k CosMx fixture is the scale extension. CosMx 10k also tests 100/500
permutations for Full and Without edge-score reuse. Its existing 20-permutation
core rows are the first point of that diagnostic.

Five independent-process technical repetitions per setting: 125 core runs,
20 additional permutation-diagnostic runs, and 25 extension runs (170 total).
Run order is deterministically shuffled within repetition blocks, seed 20261012.
All workers are sequential and pinned to the same available logical CPU.
OMP/OpenBLAS/MKL/BLIS thread counts are 1. Each worker has the same 64 GiB virtual
address-space limit and 2-hour timeout. Failure and timeout remain explicit.

## Variant definitions

| Variant | Change | Preserved |
| --- | --- | --- |
| Full | Instrumented current CSR streaming | All optimizations |
| Without spatial hashing | Two-pass, C++ blockwise exhaustive pair search, no full distance matrix | Exact radius predicate, CSR, score reuse |
| COO instead of CSR | Direct sender/receiver/distance edge arrays; linear edge traversal | Hashing, local edge set, score reuse, streaming |
| Without edge-score reuse | Recompute Hill/molecular probability for every active edge in each label permutation | CSR, hashing, invariant molecular-factor preprocessing and active-edge support |
| Dense neighborhood baseline | Actual N x N float64 distance matrix; per-LR dense-column scan | Hashing, same local edge set and probability cache, one LR at a time |

The graph builder uses temporary O(N) degree/cursor arrays for COO filling but
does not retain CSR offsets for its scoring. Ancillary graph arrays are reported
separately from primary representation bytes. Dense is a sparse-versus-dense
representation/traversal test, not a pure CSR-format measurement, nor an exact
recreation of the SpatialCellChat reference implementation. No-reuse isolates
the edge-probability cache; it does not disable all invariant preprocessing.
Component effects are conditional on the other optimizations being retained.

## Verification and measurements

Freeze actual global label-permutation indices once per input, including RNG
handling used by the saved workflow. Their prefix is reused across permutation
counts. Observed and bootstrap aggregation uses canonical sender/receiver order,
with identical thresholds and strict `bootstrap > observed` tail decisions.

Before formal timing, check the unmodified release against the instrumented
Full baseline, then each variant against Full on a tiny filtered fixture,
CosMx 5k and both Xenium 5k fixtures. Verify full active keys, maximum probability
error <= 1e-15, exact raw p-values, both p<.05 and p<=.05 significance calls,
active edge counts and exact graph endpoint/distance sets.

Every formal replicate is compared to its Full replicate on the same input,
number of permutations and frozen indices. Only numerical-validation PASS
measurements enter the summaries. Failed and pending repetitions remain visible.
Graph equivalence is audited directly in preflight; formal runs also check E.

The core ablation table contains graph, observed, permutation and inference
times (seconds), whole-process peak RSS (GiB), Time / Full and consistency.
Values are median [25th, 75th percentile], five technical repetitions.
Whole-process wall time includes loading and writing and is provided in the
exact summary/run sheets. Internal inference timing excludes input loading and
output serialization. Actual numeric-array graph bytes are measured separately;
they are not RSS. Ratios are finalized only once both groups are complete.

## Source and execution

From the repository root, install the supplied vanilla and instrumented sources:

```sh
Rscript tools/install_ablation.R .
```

The queue script `scripts/run_queue.py` uses GNU time, Linux CPU affinity,
resource limits and the recorded experimental paths. Rerunning inference
requires the prepared expression fixtures and configuration for your system.
The companion includes the saved tables, validation records, parameters and
logs. Technical repetitions measure performance under these settings; they
are not biological replicates.
