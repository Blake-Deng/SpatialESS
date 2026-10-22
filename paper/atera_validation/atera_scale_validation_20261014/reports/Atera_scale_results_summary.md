# Atera reference scale validation

Updated: 2026-10-15T12:06:04.788493+08:00

The table distinguishes successful execution, numerical agreement, and budget-limited failures.

| Tissue | Maximum reference completed N | Maximum strict N | First uncompleted N | State |
| --- | ---: | ---: | ---: | --- |
| breast | 50000 | 50000 | 100000 | allocation_failure_under_memory_budget |
| cervical | 50000 | 50000 | 100000 | allocation_failure_under_memory_budget |

Unattempted or running scales do not establish an engine limit. Each tissue stops at its first failed reference run or agreement discrepancy. The other tissue continues its predefined ladder.

Old 1k/5k results used a one-hour reference timeout; extension attempts use twelve hours. CPU, thread count, virtual address cap and scientific inputs are held fixed. ESS timings reused from the original experiment retain their source.

Raw permutation p < 0.05 is descriptive; no new FDR correction or biological validation is added.
