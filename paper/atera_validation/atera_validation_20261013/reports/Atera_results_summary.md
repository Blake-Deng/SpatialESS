# Atera supplementary experiment — live summary

Completed inference runs: 36; completed formal runs: 28; failed inference attempts: 0.

Two independent Atera WTA preview tissues are evaluated. There are 18,028 measured Gene Expression rows per input; control/probe features are excluded. Missing official annotations and QC exclusions are listed, and exact barcode files are retained.

Both engines receive the same normalized full panel, spatial coordinates, official groups, mandatory-component-complete LR database rows, spatial settings, seeds and permutation counts. This is a fixed-LR computational benchmark that bypasses overexpression-based LR selection; it is not the complete standard LR-selection workflow.

Performance: whole-process wall time includes R startup, input loading and output writing. Inference time is separately measured inside R. Peak RSS is measured by GNU time; CSR array bytes and theoretical dense matrix bytes are separate quantities. Failed/timeout runs contribute no completed performance values or numerical-agreement claims.

Sampling comparisons describe stability under changed geometry/composition. The full tissue is not biological ground truth. Replicate subsets are not independent patients; these two specimens do not establish cohort-level, causal or general biological validation.

Tables are refreshed after each task. Raw permutation p < 0.05 is descriptive and is not pathway-enrichment FDR. No new model, significance rule or biochemical validation is introduced.

Strictly passing reference comparisons: 8/8. Strict pass requires full active-key agreement, maximum absolute probability difference <= 1e-15, identical p values and no decision disagreements. Exceptions are retained per job.
