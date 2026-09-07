# Completed robust results

This directory contains the complete 50,040-feature LMM/GLM outputs, exact
feature transitions, metadata, source snapshot, input checksums and logs.
Compact summaries are one directory above. The full documentation is in
`docs/LMM_RESULTS_20260907.md` relative to the SpatialESS repository root.

Use a structured TSV reader: diagnostic messages may contain quoted newlines.
`patient_random_intercept_lmm.tsv.gz` contains all statuses; filter to `ok`
for valid inference, and then use the supplied globally adjusted q-values.
Do not recompute BH on just a selected significant or nonsingular subset
unless explicitly reporting that different testing universe.

The finalization runtime is not fitting runtime. `elapsed_seconds` is NA
because completed checkpoints were merged during a resumed invocation.
This is downstream GSE250346 model validation, not SPARKLE correction and
not a new communication-inference runtime or peak-memory benchmark.
