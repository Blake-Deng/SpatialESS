# Compatibility Repair Audit

Read RELEASE_STATUS.md at the repository root for current acceptance.
The 0.1.2 status is historical, not the current code's claim.
before_fix_diagnostics.tsv separates original/dims-only fidelity failures from
numeric-filter causal diagnostics. The numeric-filter diagnostic is not a
valid official pass. Only dims is patched in all acceptance comparisons.

The 1k stratified cell-stage cache was reused after full input-contract equality;
group inference was newly run. The portable official mini separately regenerated
that cell stage and its permutations. Other dual-track cell stages were fresh.

Sampling preserves submitted barcode order. It cannot establish the original
biological stratification rule without the collaborator's labels. Public graph
clusters are used consistently by both engines. Center samples are nested ROIs.

scripts/ preserves server-specific full audit orchestration with absolute input
paths. Use the portable benchmarks/ mini commands for independent reproduction.
All full input/output checksums are retained; large intermediate official objects
remain on the server and are not duplicated into the GitHub tree.

The original-distance audit first had a harness-only Matrix attachment error;
the final run reproduces the actual original official dimension failure. A
transient unit test formatting-option warning was corrected before final tests.
Python bridge tests pass with a SciPy upstream default-change deprecation warning.
These are distinguished from inference failures; full logs are preserved.
