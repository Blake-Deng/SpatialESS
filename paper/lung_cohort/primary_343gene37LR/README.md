# Pulmonary-fibrosis primary cohort (343 genes, 37 LR records)

This is the primary manuscript analysis: 45 lung Xenium sections from 35
individuals (10 control sections from 9 individuals; 35 PF sections from 26
individuals). Communication was inferred per section using the 343-gene panel
and 37 computable ligand-receptor records. The patient-level model was
`log1p(score) ~ condition + (1 | patient)` with global Benjamini-Hochberg
correction over valid fits.

## Files and score reconstruction

- `complete_lmm_all_50040_features.tsv.gz`: all 50,040 feature fits,
  including fit status, singular flags, coefficients, standard errors,
  confidence intervals, p values and original q values.
- `slice_scores_nonzero_long.tsv.gz`: saved nonzero, untransformed scores
  keyed by `sample_id` and `feature_id`.
- `complete_feature_universe.tsv.gz`: the fixed feature universe.
- `slice_metadata.tsv` and cell-count tables: condition, patient and
  per-section cell-type information.
- `pDC_counts_by_slice.tsv`: per-section pDC counts.
- `original_primary_model/`: original model source, saved primary LMM/OLS
  outputs, session information and input/source checksums.

For score reconstruction, form the Cartesian product of the 45 sections and
the 50,040 features and fill only absent keys with zero. This yields 2,251,800
section-feature combinations; the saved nonzero table contains 692,338 rows.
The zero-fill rule applies only within this fixed universe. Scores are raw
communication scores; model responses use `log1p(score)`. A zero is not proof
of biological absence.

The saved analysis has 39,570 valid fits, of which 17,308 are singular, and
5,147 PF-associated features after global BH correction. Valid singular fits
remain in the primary correction set. Failed and ineligible fits are retained
with their statuses and are not included in the valid-fit BH family. The model
is the frozen 2026-09-07 primary analysis; the later SpatialESS 0.1.4 release
does not imply that this cohort model was rerun under 0.1.4.

## Example support counts

For `22:Activated Fibrotic FBs:12:CD8+ T-cells:9:SPP1_CD44`, the saved score is
nonzero in 29/35 PF sections and 0/10 control sections. Counting a patient as
supported when at least one of that patient's sections has a nonzero score
gives 22/26 PF patients and 0/9 control patients; 20/26 PF patients have
nonzero scores in all of their available sections. These are descriptive
support summaries, separate from the mixed-effects condition association.
They are not a prespecified condition-consensus test and do not establish
biological absence in control tissue.

The associated beta and q value are approximately `2.3033e-4` and
`9.6331e-5`. The fitted coefficient is a condition association on the
log1p-score scale, not a direct measurement of receptor activation or a causal
effect.

## Limits and related results

The archived `final_CT` labels were transferred unchanged from the saved
manuscript input; this archive does not contain the original manual/marker
annotation procedure. A pDC count table is available, but an independent pDC
marker-expression validation figure and same-cell versus different-cell
homotypic decomposition were not found. A pDC-to-pDC feature should therefore
not be described as proven intercellular crosstalk or autocrine signaling.

Existing patient/run/TMA and other sensitivity analyses are already in the
repository under `paper/cohort_sensitivity/`; they are not duplicated here.
The distinct 99-gene official-reference fidelity analysis is in the sibling
`reference_99gene_validation/` directory. Do not combine that fidelity set
with this 343-gene biological model.
