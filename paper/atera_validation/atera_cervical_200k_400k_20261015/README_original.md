# Atera cervical central 200k / 400k extension

This independent extension uses the archived SpatialESS 0.1.4 installation and
the saved full cervical preparation. Original inputs, libraries and results are
read only. New fixtures and results are written here.

Central sampling takes the nearest cells to the arithmetic mean of all 715,512
eligible centroids, breaking distance ties by original order. There is no random
sampling. The old 100k barcodes must equal the first 100k selected barcodes; 200k
must be nested within 400k. The full measured expression matrix, normalization,
annotation, LR ordering, scale.distance and all inference parameters are preserved.

Run order: prepare both fixtures, SpatialESS 200k, SpatialESS 400k,
SpatialCellChat 200k, SpatialCellChat 400k. Reference attempts are explicitly
requested beyond the previous 100k allocation failure. A reference failure does
not cancel the other explicitly requested scale. There are at most two reference
attempts, each with a 12 h timeout. Each SpatialESS run has a 24 h timeout.
Every worker uses one CPU/thread and a 96 GiB virtual address-space budget
(not an RSS cap). GNU time records whole-process runtime and peak RSS, including
loading and output. Inferential runtime is also saved separately.

All cases use 18,028 measured genes, 3,081 fixed LR records, 100 permutations,
seed 1, effective interaction radius 35 micrometres, and contact.range 10.
This is a fixed-LR computational comparison, not the complete reference
overexpression-selection pipeline. Strict comparisons are made only for
completed identical-input pairs.

Queue: `tmux attach -t spatialess_cervical_200k400k`

Status: `reports/progress.json` and `logs/queue.log`.
Results: `tables/performance.tsv`, `tables/agreement.tsv`, and
`reports/results_summary.md`. Task logs, GNU time outputs, per-run R sessions,
input barcode files and hashes are retained. Restart the same queue to resume;
completed and failed attempts are retained, with no automatic retry of failures.

The revised Figure 1 is a pre-run snapshot and retains not-attempted markers
for cervical 200k and 400k. This queue does not automatically modify that figure
or insert incomplete timing values into it.
