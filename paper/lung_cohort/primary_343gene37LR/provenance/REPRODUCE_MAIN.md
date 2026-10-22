# Reproduce Main Summaries

The repository contains code and compact tables, not the large expression
objects. Original benchmark outputs remain on the analysis server and can
regenerate the compact tables.

## Source roots

    GSE250346 fidelity: /data/dzf/GSE250346/benchmarks/multisample_v3_exact
    GSE250346 LMM: /data/dzf/GSE250346/benchmarks/multisample/slice_glm_vs_patient_lmm
    GSE306130: /data/dzf/GSE306130/benchmarks
    GSE313006: /data/dzf/GSE313006/benchmarks/spatial_v3_headtohead
    CosMx: /data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead
    Xenium: /home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results

These are provenance locations, not runtime defaults for a new user. A public
release should replace them with user-supplied paths or archived Zenodo files.

## Parameter locks

For exact comparisons, use the same nested cells, gene panel, LR table, group
labels, distance settings, permutation count, seed and single-thread setting
for both methods. The principal SpatialESS settings are:

    interaction range: 30 um
    contact range: 10 um
    tolerance: 5 um
    inference: label permutation
    same seed within a pair
    native threads: 1
    BLAS threads: 1

## Summary command

    Rscript benchmarks/summarize_main_results.R

The script reads small TSV summaries and does not load million-cell matrices.
Full reruns use dataset-specific drivers under benchmarks/ and should use one
independent process at a time when peak RSS is being measured.


Run the patient-aware model and its slice-level comparator with:

    Rscript benchmarks/run_gse250346_lmm_comparison.R

The LMM output records `status`, `singular`, optimizer, variance components, maximum gradient, convergence messages and fit errors. Only `status = "ok"` fits receive BH-adjusted q-values.

The final CosMx status reconciliation is recorded in
results/tables/cosmx_queue_reconciliation.tsv. The paired comparison table
uses NA for fidelity metrics whenever the official process did not produce
communication records.
