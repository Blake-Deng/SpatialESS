# Atera preview-data validation

This archive documents the 10x Genomics Atera WTA FFPE breast and cervical
**preview** datasets. It supports a computational validation and scale test;
it is not a biological validation cohort or a test of a released commercial
workflow.

Code and method records are available here. Saved data and result tables
are supplied in the companion package; use the
[reproduction guide](../../docs/REPRODUCE.md) to download and mount them.

## Inputs and preprocessing

The prepared inputs contain 18,028 measured Gene Expression rows. Per-cell
counts and detected features were calculated on those measured rows. Cells
were retained when counts >= 20, detected features >= 20, and an official 10x
cell-group label was available. This yielded 169,940 breast cells from 170,057
and 715,512 cervical cells from 717,576. Exclusion combinations and overlapping
missing-label/QC exclusions are in `qc_exclusion_combinations.tsv` and the
cohort-specific preparation metadata.

Expression was normalized as `log1p(10,000 * count / total measured-GEX
counts per cell)`. Group labels came directly from the official 10x
`*_cell_groups.csv` files. Coordinates were already in micrometres and were
used with a coordinate conversion ratio of 1.

Central subsets select the N eligible cells nearest the arithmetic mean of all
eligible coordinates; distance ties retain source row order. Group-stratified
subsets allocate quotas from full-cohort group proportions using largest
remainders, sample without replacement within groups, then shuffle the result.
The three stratified seeds are 20261013, 20261014 and 20261015. Block subsets
use observed centroids near the x quartiles and y median as anchors, then select
the nearest N cells. Exact cell IDs, source-order fields, scripts and checksums
are retained in `barcode_lists/` and `preparation_metadata/`.

Both engines used identical saved fixtures, expression, coordinates, group
labels, LR order, parameters, 100 permutations and `seed.use = 1`. This was a
fixed-LR comparison: 3,081 mandatory-component-complete records from 3,234
database interactions; the default overexpression-based LR-selection workflow
was not run. The database tables actually used are in `database/`.

## Spatial and runtime settings

The saved settings use `interaction.range = 30` µm, `tol = 5`,
`contact.range = 10` µm, contact-dependent filtering enabled, forced contacts
disabled, `use.AGAN = TRUE`, `min.percent = 0.1`, `min.cells.sr = 5`, `Kh =
0.5`, and `n = 1`. The distance scale is tissue-specific and equals twice the
minimum positive inter-cell distance (approximately 2.04351143947 for breast
and 2.29750351419 for cervical tissue). Full parameter contracts are included.

Each worker was bound to one CPU with OMP/OpenBLAS/MKL/BLIS threads set to one.
The virtual-address guard was 96 GiB. Final extended SpatialCellChat attempts
were allowed 12 h and SpatialESS attempts 24 h; earlier one-hour attempts and
retries remain separately identified in the attempt ledger. GNU `time` reports
whole-process elapsed time and peak RSS, including loading and output writing.
Recorded values are individual measurements, not repeated-run medians. The
CPU/physical-memory file is a 2026-10-20 collection-time snapshot, not a
per-run hardware record.

## Results and interpretation

SpatialCellChat and SpatialESS completed paired comparisons at 1k, 5k, 10k,
25k and 50k cells in both tissues. All ten pairs met the saved strict numerical
agreement criterion. At 100k, SpatialESS completed in both tissues while
reference attempts failed allocation under the recorded budget. SpatialESS
completed the full 169,940-cell breast and 715,512-cell cervical inputs. The
full breast reference was not attempted. Cervical reference attempts at 200k,
400k and 715,512 cells failed allocation; the 715,512-cell attempt requested a
96.1 GiB dense intermediate under the 96 GiB virtual-address guard. These
budget-specific outcomes do not establish a universal limit for the reference
implementation.

`tables/`, `resource_summary_all_source_rows.tsv`, and
`all_attempts_ledger.tsv` preserve successful measurements, failures,
unattempted states, retries and their measurement provenance. Failed attempts
are not treated as successful runtimes or paired validation results. The
complete original method and audit notes are in
`METHODS_AND_AUDIT_中文.md`.

