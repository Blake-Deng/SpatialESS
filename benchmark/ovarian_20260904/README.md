# Ovarian SPARKLE integration benchmark

This directory contains compact, manuscript-ready source tables for the
2026-09-04 ovarian validation. It does not contain raw 10x data, the large H5AD
files, RCTD reference objects or full intermediate RDS files.

## Comparisons

- `SPARKLE -> official SpatialCellChat V3` versus
  `SPARKLE -> SpatialESS`: implementation fidelity, runtime and peak RSS.
- `RAW -> SpatialESS` versus `SPARKLE -> SpatialESS`: sensitivity of inferred
  communication to leakage correction.

The two comparisons answer different questions and should not be merged into a
single accuracy claim.

## Reproduce figures

From the repository root:

```bash
python scripts/make_publication_figures.py
```

Outputs are written to `figures/` as PDF, SVG and 450-dpi PNG.

See `docs/VALIDATION.md` for the complete design, values and interpretation
boundary.
