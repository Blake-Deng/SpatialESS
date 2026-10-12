# Figure provenance

All values originate in the packaged source tables. Absolute original paths
and SHA256 hashes are in `provenance/INPUT_MANIFEST.tsv`. Source-data CSVs
preserve full keys and saved precision; display abbreviations are mapped separately.

| Panel | Source file | Columns / key | Filtering and statistical rule | Display transform | n (unit) |
| --- | --- | --- | --- | --- | --- |
| Figure1B | source_data/current/cosmx_50000_paired.tsv | probability_official;probability_ess; key: sender_group_name;receiver_group_name;interaction_name | Full-key pairing; count nonpositive, NA and unpaired records; Pearson on original probability | log10 for strictly positive values | 18484 (records) |
| Figure1C | source_data/current/official_fidelity.tsv | N;probability_pearson;max_probability_difference;exact_pvalue_fraction;active_jaccard;significant_disagreements_le005; key: case | Six completed CosMx pairs; calls use raw p <= 0.05 | none | 6 (records) |
| Figure1D | source_data/current/cosmx_scaling.tsv | cells;ESS_wall_seconds;official_wall_seconds;*_status; key: cells | Curves: completed only; median of saved repeats, no inferred error bars | log10 x and y | 9 (records) |
| Figure1E | source_data/current/cosmx_scaling.tsv | cells;ESS_peak_rss_GiB;official_peak_rss_gib;*_status; key: cells | Completed only; failed resources not plotted | log10 x and y | 9 (records) |
| Figure1F | source_data/current/spatial_sparsity.tsv | N;N_squared;E;fixture;sampling_design; key: fixture | All 14 PASS spatial fixtures | log10 x and y | 14 (records) |
| Figure1G | source_data/current/spatial_sparsity.tsv | N;E;mean_out_degree;fixture; key: fixture | All 14 fixtures; original fixture order within platform; no duplication by LR filter | horizontal bars + endpoint + numeric label; fixture labels include N | 14 (records) |
| Figure1H | source_data/current/spatial_sparsity.tsv | N;graph_density;fixture; key: fixture | All fixtures; only nested CosMx values connected | log10 x and y | 14 (records) |
| Figure1I | source_data/current/cosmx_scaling.tsv | cells;ESS_status;official_status;official_failure_reason; key: cells | Paired filled; ESS-only open; official failed cross | none | 9 (records) |
| Figure2A | source_data/archived/lmm.tsv.gz;source_data/archived/volcano_labels.tsv.gz | feature_id;estimate;q_value;status; key: feature_id | status == ok; all valid points; original eight labels; no refitting or BH recomputation | -log10(q); effect ticks show x10^-3 | 39570 (records) |
| Figure2B | source_data/archived/forest.tsv.gz | feature_id;display_id;estimate;lower;upper;singular; key: feature_id | Exact archived 12 features | effect ticks show x10^-3 | 12 (records) |
| Figure2C | source_data/archived/effect_heatmap.tsv.gz;source_data/archived/selected_interactions.tsv.gz | feature_id;display_id;standardized_effect;glm_standardized_effect;q_value;glm_q;support_fraction;icc;singular; key: feature_id | Exact archived 28 rows; LMM/OLS share color range; dot q < .05; support denominator 45 | zero-centered effect fill; no new clipping | 28 (records) |
| Figure2D_slices | source_data/archived/slice_scores.tsv.gz | feature_id;sample_id;patient_id;condition;log1p_score; key: feature_id;sample_id | All 45 slices per original feature | no y transform; fixed x jitter seed 13 | 135 (records) |
| Figure2D_patients | source_data/archived/patient_scores.tsv.gz | feature_id;patient_id;condition;log1p_score; key: feature_id;patient_id | All 35 patients per original feature | fixed x jitter seed 19 | 105 (records) |
| Figure2D_feature_identities | source_data/archived/patient_examples.tsv.gz | feature_id;display_id;ligand;receptor;sender_group_name;receiver_group_name; key: feature_id | Three archived examples; exact saved identities | none | 3 (records) |
| Figure2E | source_data/archived/pathway_pair_representatives.tsv.gz | feature_id;pathway_name;sender_group_name;receiver_group_name;standardized_effect;q_value; key: feature_id | Archived 52 pathway/cell-pair representatives; no pathway-level enrichment inference | area = min(-log10(q),10) | 52 (records) |
| Figure1B_marginals | source_data/current/cosmx_50000_paired.tsv | probability_official;probability_ess; key: log10_bin_left;log10_bin_right | All 18,484 strictly positive paired probabilities | 24 equal log10 probability bins | 24 (display histogram bins) |
| Figure2A_labels | source_data/archived/volcano_labels.tsv.gz;source_data/archived/lmm.tsv.gz | feature_id;display_id;ligand;receptor;estimate;q_value; key: feature_id | Same eight frozen labels; one explicit vector arrow per complete feature key | exact target coordinate; fixed display label positions | 8 (records) |
| FigureS1A | provenance/current_docs/RESULTS.md;source_data/archived/sparkle_fidelity.tsv.gz | documented input bridge;method; key: step | Conceptual summary of existing compatibility workflow | none | 3 (steps) |
| FigureS1B | source_data/archived/bridge_regression.tsv.gz | case;comparison;max_abs_probability_difference;pvalue_exact_fraction;active_jaccard;significance_disagreements; key: case;comparison | All 12 historical bridge/engine checks; keep comparison labels | none | 12 (records) |
| FigureS1C | source_data/archived/sparkle_probability_records.tsv.gz | official_probability;spatialess_probability;official_pvalue;spatialess_pvalue; key: sender_group_name;receiver_group_name;interaction_name | All supplied rows | log10 positive probabilities; count repeated p coordinates | 37808 (records) |
| FigureS1D | source_data/archived/sparkle_probability_records.tsv.gz | official_pvalue;spatialess_pvalue; key: official_pvalue;spatialess_pvalue | Group identical saved p-value pairs; all records retained | area proportional to grouped record count | 101 (distinct p-value pairs) |
| FigureS1F | source_data/archived/raw_corrected_transitions.tsv.gz | category;records; key: category | RAW/SPARKLE input sensitivity, not engine equivalence | none | 3 (records) |
| FigureS1C_marginals | source_data/archived/sparkle_probability_records.tsv.gz | official_probability;spatialess_probability; key: log10_bin_left | All paired records | 24 equal log10 probability bins | 24 (display bins) |
| FigureS1E | source_data/archived/sparkle_fidelity.tsv.gz;source_data/archived/sparkle_probability_records.tsv.gz | active_jaccard;pvalue_exact_fraction;official_pvalue;spatialess_pvalue; key: metric | Three existing agreement metrics | saved values; calls use raw p < 0.05 | 3 (records) |
| FigureS2A | source_data/archived/cervical_spatial_metadata.tsv.gz | cell_id;x_centroid;y_centroid;group; key: cell_id | All 5k stratified cells; same five most frequent clusters + Other | equal aspect; original coordinate orientation | 5000 (records) |
| FigureS2B | source_data/current/official_fidelity.tsv;source_data/archived/cervical_fidelity.tsv.gz | case;filter;cells;max_abs_probability_difference;pvalue_exact_fraction;active_jaccard;significance_disagreements; key: case;filter | All eight matched configs; dimensions-only official comparator | none | 8 (records) |
| FigureS2C | source_data/current/cervical_*_default_paired.tsv | case;probability_official;probability_ess; key: case;sender_group_name;receiver_group_name;interaction_name | Default filters; four full-key paired fixtures | log10 positive probabilities | 204 (records) |
| FigureS2C_center_1000_marginals | source_data/current/cervical_center_1000_default_paired.tsv | probability_official;probability_ess; key: log10_bin_left | All positive paired records | 24 equal log10 probability bins | 24 (display bins) |
| FigureS2C_center_5000_marginals | source_data/current/cervical_center_5000_default_paired.tsv | probability_official;probability_ess; key: log10_bin_left | All positive paired records | 24 equal log10 probability bins | 24 (display bins) |
| FigureS2C_stratified_1000_marginals | source_data/current/cervical_stratified_1000_default_paired.tsv | probability_official;probability_ess; key: log10_bin_left | All positive paired records | 24 equal log10 probability bins | 24 (display bins) |
| FigureS2C_stratified_5000_marginals | source_data/current/cervical_stratified_5000_default_paired.tsv | probability_official;probability_ess; key: log10_bin_left | All positive paired records | 24 equal log10 probability bins | 24 (display bins) |
| FigureS2D | source_data/archived/cervical_fidelity.tsv.gz | case;filter;cells;ess_active;official_active;ess_significant_lt_005;official_significant_lt_005; key: fixture;metric;method | All eight saved configurations; significance is raw permutation p < 0.05 | grouped bars, paired methods | 32 (fixture by metric by method) |
| FigureS3A | source_data/current/spatial_sparsity.tsv | fixture;N;csr_total_bytes;dense_float64_bytes; key: fixture | Same five archived representative fixtures; CSR numeric-array payload; dense theoretical float64 | log10 y | 5 (records) |
| FigureS3B | source_data/current/spatial_sparsity.tsv | fixture;N;mean_out_degree;sampling_design; key: fixture | Four Xenium fixtures; existing graph metrics only | none | 4 (records) |
| FigureS3C | source_data/current/spatial_sparsity.tsv | fixture;N;isolated_cells; key: fixture | All 14 graph fixtures; no graph recomputation | 100 * isolated_cells / N | 14 (records) |
| FigureS3D | source_data/current/spatial_sparsity.tsv | fixture;N;csr_total_bytes;all_numeric_array_bytes; key: fixture | Nine CosMx fixtures; measured numeric arrays only | bytes / 1e6; log10 axes | 9 (records) |
| FigureS4A | source_data/archived/model_summary.tsv.gz | lmm_* counts; key: status | Mutually exclusive documented fit states; no refitting | none | 5 (records) |
| FigureS4B | source_data/archived/lmm_glm_summary.tsv.gz | glm_significant;lmm_significant;overlap_significant; key: summary | Existing BH q < .05 overlap | none | 1 (records) |
| FigureS4C | source_data/archived/pathway_pair_representatives.tsv.gz | feature_id;pair;pathway_name;standardized_effect;q_value; key: feature_id | Same 52 main-panel features; expanded labels only | area min(-log10(q),10) | 52 (records) |
| FigureS4D | source_data/archived/patient_examples.tsv.gz | feature_id;estimate;lower;upper; key: feature_id | Same three main examples, saved Wald intervals | effect ticks x10^-3 | 3 (records) |

## Interpretation

- Figure 1A is a readable editable summary of the tested workflow. The original supplied PDF remains byte-identical in assets for provenance.
- Figure 1D/E include completed measurements only. ESS values are saved v0.1.4 measurements; official measurements are saved comparator runs, not newly timed here.
- F1B marginal histograms count every paired record in equal log10 probability bins. They do not introduce simulated points or change fidelity calculations.
- Light fills in Figure 1D-F mark gaps between existing paired curves or CosMx N-squared and E; they are comparison guides, not confidence intervals.
- Figure 1G retains all 14 mean-degree values, now arranged as fixture-wise horizontal bars with endpoint symbols and printed values; input N remains explicit in each label. No descriptive CI is invented.
- Fidelity uses the original permutation-p call thresholds recorded per panel. These are not BH-FDR discoveries.
- Figure 2 uses the frozen 2026-09-07 patient-random-intercept analysis, verified against the current released model table by full key.
- Figure 2A has one explicit vector arrow per saved labeled full key. Label and arrow-target coordinates are provided in Figure2A_labels.csv.
- Supplementary record-count, fit-state and degree panels share the endpoint-bar style. S1/S2 probability marginals count the existing records; S3 isolated percentages use all 14 fixtures.
- Figure 2C prints effect / response SD; LMM and OLS have the same zero-centered color range. Dots use each saved feature-level BH q < 0.05.
- Patient means are means of the saved slice log1p scores. No group-mean CI is inferred from a condition-effect CI.
- Figure 2E/S4C use saved LR-feature BH q, not pathway-enrichment FDR. Blank locations have no saved representative.
- S1 bridge/official SPARKLE comparisons are the identified archived series. They are not relabeled as a newly rerun comparator.
- S2 official comparator has the documented dimensions-only repair. The tissue map uses public graph-cluster labels, not new cell-type annotation.
- N squared is the potential ordered-pair matrix; E counts stored directed non-self radius edges. All tested graphs use effective radius 35 micrometres.
- S3 dense memory is theoretical float64 matrix storage; CSR is measured numeric-array payload, neither equals peak RSS.
- The CSPG4 display uses the audited receptor ITGA3+ITGB1 while preserving the original interaction key and statistics.
- The full 66-comparison table retains the one non-exact permutation-p result; the main CosMx ladder is not a claim of exact agreement for every dataset.
