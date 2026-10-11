# SpatialESS data figures: compact editorial style, revision 2

This package contains two main-figure data masters and four supplementary figures.
The data and statistical selections are frozen; only presentation has changed.

## Files to open

- `figures/Figure1_data.pdf` and `.svg`: editable panels B-I, with the user-owned A slot empty.
- `figures/Figure2.pdf` and `.svg`: complete redrawn data figure, including panel A.
- `figures/FigureS1` through `FigureS4`: editable supplementary masters.
- `figures/*.png`: 600 dpi images at the declared physical size.
- `preview/*_from_pdf.png`: previews rendered from the actual PDF masters.
- `panels/`: separate editable PDF/SVG panels.
- `source_data/Source_Data.xlsx`: panel data, provenance, stable feature identities and all 66 official comparisons.
- `FIGURE_PROVENANCE.md`, `reports/QA_REPORT.md` and `FIGURE_CAPTIONS.md`: definitions, checks and captions.

Figure 1 retains the neighbor-graph analysis in panels F-H. Figure 2 retains all
28 selected heatmap rows, 12 forest features, eight original volcano labels and
three original patient examples. The full feature-key mapping is supplied.

## Protected Figure 1A

`assets/Figure1A_user_original.pdf` is byte-identical to the supplied workflow.
No elements, fonts, colors or internal geometry were edited. The slot is
`[x=2, y=2, width=104, height=46.404]` mm from the top-left of a 180 x 190 mm page.
`config/Figure1A_composition.json` records the proportional placement.

`Figure1_with_original_A_reference.pdf` is a reference composition only. The
original A's text scales below the target minimum in this slot; the user will
complete A. The editable data master deliberately leaves that slot blank.

## Reproduce from the packaged data

```sh
python3 -m pip install -r requirements.txt
python3 scripts/run_all.py
```

Arial is used when installed. Set `SPATIALESS_FONT_DIR` to an available Arial font
directory if required; the default is `~/.local/share/fonts/arial`. Font files are
not redistributed. A missing Arial installation is reported and uses Liberation
Sans. Python package versions and the actual font path are saved in `logs/`.

`run_all.py` verifies the original input hashes, redraws from staged tables,
renders and audits the exports, and writes the reports. It does not fit models,
recompute BH values or rerun communication inference. `prepare_sources.py` and
`extract_pairs.R` document the original source staging and full-key binary
extraction; these optional steps require the original project and saved RDS files
and are not needed to reproduce the figures from this package.

## Editing

Physical panel rectangles are in `config/layout.json`. Font, colors and line
settings are in `scripts/theme.py`; plot-specific settings are in
`scripts/redraw.py`. PDF/SVG text, axes, curves, hexagons, heatmap cells and legends
are vectors. The dense volcano point layers are intentionally rasterized at
600 dpi; their text and axes remain editable. SVG editing requires the selected
font on the receiving computer. PDF fonts are embedded.

The new Figure 1 uses a full-width fidelity matrix, wider performance plots,
comparison fills between existing paired curves, actual probability marginals,
and fixture-wise mean-degree bars alongside the pair-space and graph-density
curves. Fills are comparison guides, not uncertainty intervals. No new inference
was performed; marginal bars count the existing paired records.

Figure 2A has a separate full-LR label column with eight explicit editable arrows.
Each arrow target is verified against the full feature key, saved effect and BH q;
the physical coordinates and feature identity are in `Figure2A_labels.csv`.
Figure 2B retains the effect-bar / saved-CI design. Blue/rose effect colors,
pathway-grouped heatmap rows with saved values, and larger data marks are kept.
No numerical results or major feature selections were changed. The degree panel
now uses explicit fixture rows, endpoint marks and numeric labels, retaining all
14 inputs. All four supplementary figures use the same fidelity-tile, marginal
distribution, endpoint-bar and comparison-fill visual language. Descriptive
bars have no invented error intervals. Decimal MB graph-array
payloads and GiB measured peak RSS
remain distinct. Supporting plots do not assert an unperformed spatial-null test.
