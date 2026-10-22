# Figure 1 QA report

## Outcome

The new Figure 1 is a 180 × 250 mm vector composition. The author-provided workflow was placed across a full-width 176 × 78.52 mm panel at the top, preserving aspect ratio and all internal content. The input asset hash is unchanged and matches `SpatialESS_workflow_comparison.pdf`. The estimated smallest source lettering in the enlarged A panel is 4.38 pt; this is retained original lettering, not newly typeset text.

Panels B–I were regenerated from the frozen project source tables. Atera was added to runtime, peak-memory, and completion panels. Panel C includes the 10 Atera paired comparisons, all of which passed the stored strict criterion; the largest probability difference was 2.9837 × 10−16, permutation p values were identical, and there were no significance-call disagreements.

## Atera scope represented

- Inputs used 18,028 measured genes, 3,081 fixed ligand–receptor records and 100 permutations.
- SpatialESS completed the breast input at 169,940 cells and cervical input at 715,512 cells.
- SpatialCellChat completed paired comparisons through 50,000 cells in both tissues.
- Both 100,000-cell SpatialCellChat attempts failed during allocation under the configured memory budget. At 715,512 cervical cells, the full-input reference attempt failed requesting a 96.1 GiB dense intermediate under the 96 GiB virtual-memory limit (31.89 s; peak RSS 12.96 GB); no reference result was generated. The full-breast reference was not attempted. These outcomes do not establish a universal engine limit.
- Atera memory is decimal GB. CosMx RSS values originally recorded in GiB were converted to decimal GB for the shared axis.

## Validation

Automated results are in `composition_checks.tsv`, `export_checks.tsv`, `source_data_checks.tsv` and `data_checks.tsv`. Checks cover page dimensions, in-bounds text, embedded fonts, absence of Type 3 fonts, 600 dpi previews, vector SVG content, complete Atera source rows, all ten strict paired comparisons, and byte identity of frozen raw inputs.

The composed PDF contains selectable text and vector graphics. The SVG retains editable text for panels B–I and vector paths for the supplied A workflow; it contains no raster image replacement. All panel-level D/E source rows retain completion state, including failures and unattempted scales, while failed attempts are excluded from successful runtime and memory curves.

Panels D and E use blue circles for SpatialESS and orange squares for SpatialCellChat. Solid, dashed and dotted lines distinguish CosMx, Atera breast and Atera cervical; shared legends identify both encodings. Neutral shading marks only input sizes with completed measurements from both methods. The restored D runtime axis starts at 8 s; smaller CosMx measurements remain in source data but fall below the displayed range. No measured values were moved to separate nearby curves. Panels F–H have a shared legend for the four platform/sampling-design series and the N² reference curve. Panel I uses equally spaced size ticks at 1k, 5k, 10k, 25k, 50k, 100k, 200k and 700k; 200k is followed directly by 700k and the explanatory x-axis title is omitted. The 400k cervical and nominal 450k (443,515-cell CosMx) categories are hidden in I, with all run states retained in source data; the 443,515-cell measurements remain in D/E. The 169,940-cell breast category and two 250-cell CosMx states are also hidden only in I. At 700k, the same reference-failed cross marks the attempted run, which failed requesting a 96.1 GiB dense intermediate under the 96 GiB virtual-memory limit; no reference result was generated. The D, E and I panel previews were refreshed from their current vector PDFs.

The Atera runs use a fixed-LR comparison setup. They should not be described as running SpatialCellChat's complete overexpression-based LR selection workflow. The failed 700k attempt produced no probabilities or p values and is not a paired comparison.

The frozen baseline inputs and author workflow retain byte identity in the provenance checks; the cervical extension table was synchronized to terminal run records. All 35 numerical checks, 27 export checks and 56 baseline source checks passed. A color-vision-deficiency simulation was not performed.
