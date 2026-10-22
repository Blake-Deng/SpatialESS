# SpatialESS component ablation, 2026-10-12

An independent experiment from the SpatialESS 0.1.4 submission snapshot. Original
package and manuscript files are not modified. All variants retain the object
protection fixes, molecular formula, group filters, self-edge convention and
permutation tail criterion. No new biological model or significance filtering
is introduced.

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

## Background execution and outputs

Tmux session: `spatialess_ablation_20261012`.

```sh
tmux attach -t spatialess_ablation_20261012
```

Resume safely from this directory:

```sh
python3 -u scripts/run_queue.py
```

The queue updates `tables/SpatialESS_ablation_table.xlsx`, `.docx`, `.md`,
`ablation_summary.tsv`, `.csv` and `ablation_runs.tsv` after each completed job.
`reports/progress.json`, `reports/table_status.json` and `logs/queue.log` expose
progress. At the end it writes `SpatialESS_ablation_results.zip`, containing
the final table, scripts, source variants, design, validation records and logs.
Prepared expression inputs remain local and are not put in the downloadable ZIP.

To rebuild in an R environment with Rcpp and Matrix:

```sh
Rscript -e 'Rcpp::compileAttributes("source/instrumented")'
R CMD INSTALL --preclean -l lib_vanilla source/vanilla
R CMD INSTALL --preclean -l lib source/instrumented
```

This experiment has not been added to the submission manuscript or GitHub
package. Performance conclusions must follow the completed measurements.
