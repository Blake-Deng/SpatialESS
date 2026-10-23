# Checks

- Primary scores, coefficients and Figure 1A unchanged.
- Full feature keys align across all models.
- R and Python HC3 results agree within 1e-10; independent statsmodels checks passed.
- BH excludes invalid/warning fits and is recomputed over each full eligible universe.
- Nonsingular subset re-BH is labelled separately.
- Tissue-affect confounding and run/TMA redundancy are recorded, not silently resolved.
- Balanced subset ranges are not confidence intervals.
- Font: /home/dzf/.local/share/fonts/arial/Arial.TTF; figure width 180 mm, height 192 mm, PNG 600 dpi.
- Text bounding-box check: 0 outside-canvas text elements.
- PDF/SVG retain vector text and axes; the dense hexbin layer alone is rasterized.
- PDF rendering and visual inspection must be recorded separately.
- Limitations: only nine controls; no independently identifiable anatomical-region effect; approximate robust inference, no causal or experimental validation.

## Completed inspection

The final PDF was rendered with pdftoppm and visually inspected. All a–f letters, legends, points, intervals, labels and zero-centred heatmap are visible; no clipped labels or missing glyphs remain. Arial regular and bold are embedded TrueType. The original 45 profile, metadata and model-source hashes match the primary archive. A full independent rerun from the bundled prepared matrix reproduced both adjusted LMM tables exactly. Word XML/media checks are complete; rendered Word pagination remains uninspected.
