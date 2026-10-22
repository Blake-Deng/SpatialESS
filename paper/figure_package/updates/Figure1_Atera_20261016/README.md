# SpatialESS Figure 1 with Atera scaling

This revision enlarges the author-supplied Figure 1A to a full-width top row and redraws panels D, E and I with the saved CosMx and Atera benchmark measurements. Panel C also reports the 10 paired Atera consistency checks. The original workflow PDF is placed as an unchanged vector asset; its source SHA-256 matches `SpatialESS_workflow_comparison.pdf`.

Panels D/E use the original SpatialESS-blue / SpatialCellChat-orange encoding, with solid, dashed and dotted lines for CosMx, Atera breast and Atera cervical. Their shared legends identify method and cohort; paired-size shading is retained. Panel I reads the cervical 200k and 400k completion states from the validated extension run table. Its size ticks stand alone without an explanatory axis title. The 169,940-cell breast category is omitted from I only; its measurements remain in D/E and source data. Cervical SpatialESS completed both 200k and 400k in the run table; the 400k category is now hidden in I, along with the nominal 450k (443,515-cell CosMx) category. Their run states remain in source data; the CosMx 443,515-cell measurements remain in D/E. The 700k reference attempt failed on a 96.1 GiB dense intermediate request; the figure uses the same reference-failed cross as earlier failed attempts, and the exact failure is retained in source data. The visible 200k tick is followed directly by 700k.

## Main deliverables

- `figures/Figure1.pdf`, `figures/Figure1.svg`, `figures/Figure1.png`: composed full figure, editable vector master and 600 dpi preview.
- `figures/Figure1A.pdf`, `figures/Figure1A.svg`, `figures/Figure1A.png`: standalone enlarged workflow panel.
- `figures/Figure1_data.pdf`, `.svg`, `.png`: data-only master; the A slot is blank by design.
- `panels/Figure1B.pdf` through `Figure1I.pdf` and matching SVG files: editable panel exports.
- `source_data/Figure1_Source_Data.xlsx` and `source_data/panels/`: panel-level data, including the 46 CosMx/Atera runtime and memory states, 10 Atera paired comparisons, and completion status.
- `Figure1_caption.md`: panel-specific caption with the comparison and failure scope.
- `reports/QA_REPORT.md` and `reports/VISUAL_REVIEW.md`: numerical, provenance, export and visual checks.

## Reproduce

From this directory, use the Python environment recorded in `logs/environment.json`:

```sh
python scripts/redraw.py --figure1-only
python scripts/compose_figure1.py
```

The figure uses saved analysis outputs only; the scripts do not rerun inference or recompute statistical tests. Atera uses 18,028 measured genes, 3,081 fixed ligand–receptor records and 100 permutations. The Atera runs are fixed-LR benchmark results and should not be described as a full overexpression-selection workflow.
