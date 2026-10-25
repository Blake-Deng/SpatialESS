# Figure 1: numerical agreement, scaling and spatial sparsity

The workflow panel is supplied as an unchanged vector asset, placed across the
top row. Panels D/E show saved CosMx and Atera measurements; panel C includes
ten paired Atera consistency checks. The Atera inputs contain 18,028 measured
genes and 3,081 fixed LR records, with 100 permutations.

The blue/orange encoding distinguishes SpatialESS and SpatialCellChat; line
styles distinguish CosMx, Atera breast and Atera cervical. Completion panel I
shows the 200k and 700k categories. The cervical 400k and CosMx 443,515-cell
measurements remain in source data; CosMx measurements also appear in D/E.
The 700k reference attempt failed on a 96.1 GiB dense intermediate request.
Failures retain their recorded resource budgets and are not successful runtimes.

## Files

Exports and source data are supplied through the companion data package;
follow the [reproduction guide](../../../../docs/REPRODUCE.md) before using
these paths:

- `figures/Figure1.pdf`, `.svg`, `.png`: composed figure and preview.
- `figures/Figure1A.pdf`, `.svg`, `.png`: workflow panel.
- `figures/Figure1_data.pdf`, `.svg`, `.png`: data panels without panel A.
- `panels/Figure1B.pdf` through `Figure1I.pdf`: individual panel exports.
- `source_data/Figure1_Source_Data.xlsx` and `source_data/panels/`: source tables.
- `reports/QA_REPORT.md` and `reports/VISUAL_REVIEW.md`: verification records.

## Reproduce

From the repository root, install `requirements-reproduce.txt`, then run:

```sh
python3 tools/redraw.py Figure1 --materials ../SpatialESS_reproduction_20261023
```

Figure composition requires Poppler. The script verifies fixed input hashes
and plots saved values without inference or statistical testing. Fixed-LR
Atera benchmarks do not test overexpression-based LR selection.
