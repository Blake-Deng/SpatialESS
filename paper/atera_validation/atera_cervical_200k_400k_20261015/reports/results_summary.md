# Cervical 200k / 400k supplemental runs

Updated: 2026-10-15T23:08:07.073554+08:00

| N | Engine | Status | Whole-process seconds | Peak RSS (decimal GB) |
|---|---|---|---|---|
| 200000 | SpatialESS | completed | 1319.8 | 3.529719808 |
| 200000 | SpatialCellChat | allocation_failure_under_memory_budget | 36.49 | 33.293033472 |
| 400000 | SpatialESS | completed | 2440.91 | 6.170734592 |
| 400000 | SpatialCellChat | allocation_failure_under_memory_budget | 84.85 | 65.585041408 |

Inputs: 18,028 genes, 3,081 fixed LR records, 100 permutations, seed 1.
Both methods use the identical saved input and one CPU/thread, sequentially.
Per-worker budget: 96 GiB virtual address space; SpatialCellChat 12 h, SpatialESS 24 h.
Resource measurements for unsuccessful attempts are failure observations, not successful runtime results.
Compare numerical agreement only if both methods complete. No overexpression-based LR selection is run.
The delivered figure remains a pre-run snapshot with 200k/400k marked not attempted.
