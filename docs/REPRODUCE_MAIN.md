# Reproduce Main Summaries

The repository contains code and compact tables, not the large expression
objects. Original benchmark outputs remain on the analysis server and can
regenerate the compact tables.

## Source roots

    GSE250346: /data/dzf/GSE250346/benchmarks/multisample_v3_exact
    GSE306130: /data/dzf/GSE306130/benchmarks
    GSE313006: /data/dzf/GSE313006/benchmarks/spatial_v3_headtohead
    CosMx: /data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead
    Xenium: /home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results

These are provenance locations, not runtime defaults for a new user. A public
release should replace them with user-supplied paths or archived Zenodo files.

## Parameter locks

For exact comparisons, use the same nested cells, gene panel, LR table, group
labels, distance settings, permutation count, seed and single-thread setting
for both methods. The principal integrated settings are:

    contact radius: 15 um
    diffusion radius: 35 um
    spatial block size: 100 um
    inference: fixed permutation
    same seed within a pair
    native threads: 1
    BLAS threads: 1

## Summary command

    Rscript benchmarks/summarize_main_results.R

The script reads small TSV summaries and does not load million-cell matrices.
Full reruns use dataset-specific drivers under benchmarks/ and should use one
independent process at a time when peak RSS is being measured.


The final CosMx status reconciliation is recorded in
results/tables/cosmx_queue_reconciliation.tsv. The paired comparison table
uses NA for fidelity metrics whenever the official process did not produce
communication records.
