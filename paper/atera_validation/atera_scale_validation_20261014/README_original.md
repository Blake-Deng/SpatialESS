# Atera reference scale validation, 2026-10-14

The reference ladder extends the saved Atera fixed-LR benchmark from 5k to
10k, 25k, 50k, 100k and full eligible tissue (breast 169,940; cervical 715,512).
All inputs retain 18,028 measured genes, 3,081 fixed LR rows, 100 permutations,
seed 1, full measured-panel normalization and the saved radius/contact settings.
The original overexpression LR selection is bypassed, as in the saved benchmark.
The reference percentage-formatting behavior is preserved. No new model or FDR.

10k is the exact first 10,000 entries of the saved central 25k input, itself
ordered by squared distance to the full eligible coordinate mean with original
indices breaking ties. Its expression, coordinates, group labels, LR definitions
and constants are checked against the original saved input and the 5k prefix.
The saved ESS 25k/50k/100k/full results are reused read only. Comparisons are
written to `comparisons/`, never into the original results directory.

Workers run sequentially on logical CPU 0. OMP/OpenBLAS/MKL/BLIS thread counts
are one. Reference timeout is 12 hours per task; ESS timeout is 24 hours.
The virtual-address-space limit is 96 GiB; this is not an RSS cap or a measured
memory value. The finite plan has at most ten new reference runs, so total
reference wall budget is at most 120 hours. Each tissue stops at its first
reference failure or agreement discrepancy; the other tissue continues.

Failures and timeouts are recorded without automatic retries. Allocation errors,
segfaults, implementation errors, unknown signal-9 terminations and timeouts
have separate statuses; signal 9 alone is not evidence of kernel OOM.
The early package-path setup errors are archived and excluded as pre-inference
setup failures. Checkpoints validate the saved input MD5 before reuse.

Tables refresh after each task. The output reports maximum completed reference
N, maximum strict-agreement N and first uncompleted reference N. A running,
unattempted or failed-by-budget scale is not an intrinsic engine performance
limit. Completed outputs and numerical agreement are separate statuses.
Single measurements have no invented uncertainty intervals. Successful speedups
exclude failed measurements. Old 1k/5k baseline runs used a one-hour timeout.

Tmux session: `spatialess_atera_scale_20261014`.

```sh
tmux attach -t spatialess_atera_scale_20261014
/home/dzf/miniforge3/bin/python scripts/status.py
/home/dzf/miniforge3/bin/python scripts/queue_scale_v2.py
```

The last command resumes checkpoints; the queue lock prevents duplicate runners.
`tables/Atera_scale_validation.xlsx`, TSV/CSV, `reports/`, GNU time logs and
per-comparison exception tables persist on disk. At bounded-plan completion,
`Atera_reference_scale_validation_results.zip` is generated automatically.
The earlier delivery ZIP and original analyses remain unchanged.
