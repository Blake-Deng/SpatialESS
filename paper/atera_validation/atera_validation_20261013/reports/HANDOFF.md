# Durable queue handoff — 2026-10-13

Authorized objective: run supplementary Atera WTA experiments unattended.

- tmux session: `spatialess_atera_20261013`, windows `queue` and `status`.
- Runner initially started at 13:09 China time; PID 3821002. Always inspect
  `progress.json` and actual processes rather than relying on this historical PID.
- Isolated SpatialESS 0.1.4 installed successfully from the repaired 20261012
  submission source. Source hashes and reference implementation log are retained.
- Raw breast/cervical downloads and validation reside under
  `/home/dzf/cellchat_acceleration/external/atera_wta_20261013/`.
- Breast preparation completed: 169,940 eligible cells, 18,028 measured gene
  rows, 3,081 mandatory-component-complete LR records; 117 cell exclusions.
- Breast 250-cell/32-LR/5-permutation paired smoke passed: 30 active features,
  Jaccard 1, maximum probability error 0, maximum p-value error 0.
- Breast full-LR 1k/20-permutation ESS pilot completed: inference 7.223 s,
  whole process 8.591 s, measured RSS 0.273 GB. These are pilot measurements,
  not the formal 100-permutation benchmark. Its reference was still running
  at this handoff; do not assert full-LR agreement until comparison exists.
- Plan: 36 inference jobs including smoke/pilots; 24 formal ESS jobs and four
  formal small reference jobs. Inspect `config/run_plan.tsv` for exact jobs.
- Source data, parameters, source code snapshot, exclusions and barcode files
  are isolated here; no original data or manuscript files are overwritten.

Check `scripts/status.py`, `logs/queue.log`, `reports/queue_error.json` if present,
and `tables/Atera_supplementary_results.xlsx`. Tables refresh after every task.
Failed tasks are retained rather than automatically retried. Diagnose failures
before relocating their exit markers to an attempt archive for a deliberate retry.
Do not run a second queue while the current lock is held.

The task uses full measured-panel normalization with a fixed measurable LR set
for both engines. This bypasses overexpression-based selection and must be
described as a fixed-LR computational benchmark. No equivalence is asserted for
ESS-only large inputs. See README for sample-stability and biological limits.
