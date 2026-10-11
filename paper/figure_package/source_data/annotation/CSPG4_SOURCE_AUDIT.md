# CSPG4 molecular-identity audit

## Decision

Retain the computed molecular definition CSPG4-(ITGA3+ITGB1) and all primary
statistics. Correct human-readable labeling and preserve the original interaction
and feature keys. No feature is excluded or replaced.

## Trace

- Archived input: `/data/dzf/GSE250346/benchmarks/multisample/computable_lr.rds`.
- `interaction_name`: `CSPG4_ITGA2_ITGB1`.
- `interaction_name_2`: `CSPG4 - (ITGA2+ITGB1)`.
- `ligand`: `CSPG4`; `receptor`: `ITGA3_ITGB1`.
- The receptor complex resolves to ITGA3 and ITGB1. ITGA2_ITGB1 is a separate
  complex. See `CSPG4_ORIGINAL_LR.tsv` and `CSPG4_COMPLEX_DEFINITIONS.tsv`.
- The historical CellChat 1.6.1 database reproduces the same field mismatch.
- `prepare_cellchat_lr_components()` in `lr_components.R` builds component
  lists from the ligand/receptor columns. `spatial_v3_exact.R` passes those
  indices into the inference kernel, then attaches `interaction_name` as a
  record identifier. Copies and hashes are in `source_snapshot/` and
  `CSPG4_CODE_SOURCES.tsv`.
- All 45 archived per-slice record files contain the original key with ligand
  CSPG4 and receptor ITGA3_ITGB1. All pass; see `CSPG4_45_SLICE_RECORD_AUDIT.tsv`.
- The LR row cites PMID 30213051. Its full text explicitly describes CSPG4
  association with alpha3beta1 integrin. The source XML, bibliographic response
  and relevant paragraphs are retained as `CSPG4_PMID30213051_fulltext.xml`,
  `CSPG4_reference.json` and `CSPG4_REFERENCE_EVIDENCE.txt`.

Reference: Schiffer et al. Int J Mol Sci. 2018;19(9):2724.
doi:10.3390/ijms19092724. https://pmc.ncbi.nlm.nih.gov/articles/PMC6164575/

## Statistical implications

The audit supports an identifier/display-name discrepancy. It does not indicate
that the archived computation used ITGA2 instead of its documented ITGA3 receptor.
There is therefore no molecular-input replacement and no reason to recompute the
primary scores, fits, intervals, p-values or global BH correction for this display
repair. The 28 heatmap rows and 12 forest rows retain their original selection.
The evidence concerns molecular annotation, not lung-specific biological activity.

`CSPG4_DISPLAY_MAPPING.tsv` maps the original fields to the corrected display.
`FIGURE2_CSPG4_DISPLAY_AUDIT.tsv` identifies the retained highlighted features.
Original source tables are preserved byte-for-byte.

## Counterfactual sensitivity, audit only

Before the evidence audit was resolved, excluding this LR family was evaluated
as a sensitivity check. Of 1,029 related model features, 633 were valid and 68
were original BH hits. Excluding them would leave 38,937 valid fits with 5,079
hits under the original correction; recomputing BH in that reduced universe
would yield 5,081 hits (two added, none lost). Those files are explicitly named
`LMM_BH_SENSITIVITY_EXCLUDING_CSPG4.tsv.gz` and
`CSPG4_BH_SENSITIVITY_SUMMARY.tsv`. This hypothetical exclusion was NOT applied
to any manuscript figure, reported primary result or main source table.
