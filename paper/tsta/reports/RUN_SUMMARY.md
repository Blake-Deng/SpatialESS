# TSTA SpatialESS run

Updated: 2026-10-15T16:40:17+08:00

All 1258 samples are included in the queue. Status: {'COMPLETE': 1258}.

Completed retained spots: 3,648,329.

Saved array IDs are converted with the archived SecAct calWeights() formula. These are conventional Visium spots, not segmented single cells. Groups are per-sample expression graph clusters, not cell-type annotations. All complete-coverage human database LR records are tested with 100 global label permutations and a 200 um spot radius; cell-contact restriction is disabled.

Each worker is CPU-pinned and single-threaded. GNU time records complete process wall time and peak RSS; inference time is separately measured inside R. Sixteen concurrent sample workers measure batch throughput under contention, not isolated speedup against SpatialCellChat. No reference-equivalence or new cell-type-specific biological finding is established by this feasibility run.
