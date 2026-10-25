# SpatialESS: condition-support classification

This post hoc descriptive analysis uses saved scores from the 343-gene / 37-LR pulmonary-fibrosis cohort: **50,040 complete sender → receiver × LR features, 45 sections and 35 patients (26 PF, 9 Control)**. SpatialESS inference, primary model fits and BH corrections are not repeated or modified.

## Definitions

A saved section score greater than zero is counted as observed support. In the main analysis a patient is supported if at least one available section has a nonzero score. The denominators remain all 26 PF and 9 Control patients, rather than a post hoc subset with detectable cell groups.

| Label | Rule at the primary 50% threshold |
|---|---|
| PF-restricted | At least 13/26 PF patients supported, and 0/9 Control patients supported |
| Control-restricted | At least 5/9 Control patients supported, and 0/26 PF patients supported |
| Cross-condition shared | At least 13/26 PF and 5/9 Control patients supported |
| Other | All remaining features |

These four labels are mutually exclusive and exhaustive. "Restricted" means restricted **observed support in this dataset**, not biological absence in the other condition. "Shared" is the operational cross-condition consensus category, not equality of score distributions or an equivalence test. Within-condition recurrence is separately quantified by patient support fractions. A shared feature may have a significant condition association.

Threshold checks use 40% (11 PF / 4 Control patients) and 60% (16 PF / 6 Control patients). These descriptive thresholds were selected post hoc; they were not prespecified for the primary cohort analysis. Every original feature is retained irrespective of fit status. Original `estimate`, `std_error`, `statistic`, `p_value`, `q_value`, confidence intervals, `status` and `singular` are copied from SD4. The separate `primary_association_label` uses original status-ok primary q < 0.05 and the original beta direction; no new P or q is computed.

## Files to read first

- [Summary workbook](results/Condition_support_summary.xlsx): counts, threshold transitions, one-section sensitivity, patient units and the manuscript examples.
- [Full classification table](results/all_50040_feature_classification.tsv.gz): full source table with all SD4 fields, all three labels and per-feature sampling probabilities. Decompress to TSV for spreadsheet software, or read directly with R/Python.
- [Heterotypic subset](results/heterotypic_feature_classification.tsv.gz): convenience subset with sender group different from receiver group; no new statistical filtering or BH correction.
- [Results in Chinese](reports/RESULTS_zh.md): results and interpretation in Chinese.
- [Methods and results](reports/METHODS_and_RESULTS.md): manuscript-ready English definitions, methods and results.
- `reports/validation.tsv`, `reports/independent_validation.json`: checks against SD4 and independent validation of all exact sampling probabilities.

## One-section-per-patient sensitivity

In each random repetition, one existing section is selected independently and uniformly for each patient. There are 1,000 repetitions with NumPy PCG64 seed **20261024**. A section can appear again in another repetition. All 35 patients remain in each repetition; patients are not bootstrap-resampled, models are not refitted, and primary q values are not reassessed.

There are only **768** possible one-section-per-patient selections: one patient has 3 sections, eight have 2, and 26 have 1 (3 × 2^8). All 768 equally weighted combinations are additionally enumerated. Exact per-feature category probabilities and exact ranges therefore have no Monte Carlo sampling error conditional on the existing saved sections. They are not estimates of biological reproducibility or confidence intervals across a patient population.

`one_slice_exact_retention_50pct` is the proportion of all 768 choices retaining a feature's original 50% class. A value of 1 means every available choice retains the label; 0.5 means half do. Analogous values are supplied for 40% and 60%. `one_slice_exact_p_*` gives probabilities for each of the four classes. `one_slice_mc_*` reports the 1,000 random repetitions at 50%. Modal class ties are explicitly flagged; the category order in the script only resolves the displayed modal label, not the underlying probabilities. Pooled category-count percentiles are descriptive selection distributions, not inferential confidence intervals.

`exact_slice_selections.tsv.gz` and `random_slice_selections.tsv.gz` retain every patient/section choice. `one_slice_counts_by_selection.tsv.gz` retains per-selection category counts at every threshold. `one_slice_expected_transitions.tsv` shows expected transitions from the full-section patient-any baseline to one-section classifications; "features_ever_in_target" can overlap across target classes and must not be summed as exclusive counts.

## Reproduce

Inputs are included as small saved score/metadata tables; no expression matrices or inference RDS are required. Original source paths, SHA-256 values and environment are in `config/design.json` and `reports/run_provenance.json`.

The scripts require Python ≥ 3.11. Install dependencies from `requirements.txt`.

Validated environments:

| Platform | Python | NumPy | pandas | openpyxl |
| --- | --- | --- | --- | --- |
| Ubuntu Linux x86_64 | 3.13.13 | 2.4.6 | 3.0.5 | 3.1.5 |
| macOS arm64 | 3.11.2 | 1.26.0 | 2.2.3 | 3.1.5 |

The [reproduction validation record](reports/reproduction_validation.json)
checks all 12 saved result tables against a complete analysis rerun.

From this analysis directory (the commands write into `results/` and `reports/`):

```sh
python3 -m pip install -r requirements.txt
python3 scripts/analyze_condition_support.py
python3 scripts/validate_outputs.py
python3 scripts/write_report.py
```

The independent validator derives patient-level support-count distributions using a Poisson-binomial dynamic program, rather than enumerating section choices, and checks all four category probabilities for all features at all thresholds. The main script verifies agreement with all four SD4 support-count fields and preservation of all original columns. Tiny text-serialization differences in floating-point values are recorded separately; original inputs are supplied byte-for-byte.

## Interpretation limits

Zeros may reflect missing cell populations, missing local adjacency, targeted-panel coverage or filtering, as well as lower scores. This analysis does not separate those mechanisms, adjust their denominators, validate annotations or establish causal signaling. It does not resolve remaining run/tissue-region confounding. The single-section check addresses unequal section contributions within these saved data, not all cohort-design limitations. The 50% Control-restricted class is small and is not a basis for a symmetric disease-versus-control biological claim.

[Primary cohort](../primary_343gene37LR/README.md) ·
[Materials index](../../../docs/SUBMISSION_CONTENTS.md) ·
[Reproduction guide](../../../docs/REPRODUCE.md)
