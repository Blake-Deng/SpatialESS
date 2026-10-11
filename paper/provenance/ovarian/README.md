# Ovarian evidence and comparison scopes

The spatial data are the 10x Visium HD Human Ovarian Cancer FF feature slice,
not the scFFPE single-cell reference. The scFFPE matrix and FLEX annotations
provide the reference and its labels. `rctd/` contains the archived annotation
script, converter, spacexr DESCRIPTION, original completed RAW run log, labels
and actual mapping. The recovered spacexr version is 2.2.1. The script uses
doublet mode. Of 16,247 cells, 16,213 map to first_type; 34 remain Unknown.
RAW and corrected expression inputs share these labels.

The original RCTD annotation run did not save sessionInfo or a lockfile.
DESCRIPTION plus the contemporaneous log are the available version evidence;
the record-verification session is a later environment, not that original one.
Absolute paths in files ending `_original` are preserved archival evidence.
The portable `recheck_label_transfer.py` verifies the saved labels/mapping.

`official_pair/` contains the archived 4 September full-input official/ESS
comparison shown in Fig. S1c-e (37,808 records; 8,417 raw p < 0.05 calls).
`bridge_regression/` contains six input cases and 12 pairwise comparisons shown
in S1b, with original record, contract, RNG and session artifacts. These are
historical compatibility results, not 12 independent biological datasets.
The six current-release bridge summaries are in results/tables/bridge_results.tsv.

None of these ovarian comparisons is part of the 66 rows in
results/tables/official_fidelity.tsv. Run the packaged R record verifier to
check the archived pairs and separately compare current 0.1.4 with the saved
official full-input result. Full keys use sender, receiver and interaction.
