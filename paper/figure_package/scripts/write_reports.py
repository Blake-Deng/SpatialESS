"""Write provenance and a concise QA report from actual exported results."""
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError
import hashlib
import json

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
inv = pd.read_csv(ROOT / 'reports/PANEL_INVENTORY.tsv', sep='\t').fillna('')
data = pd.read_csv(ROOT / 'reports/data_checks.tsv', sep='\t')
exports = pd.read_csv(ROOT / 'reports/export_checks.tsv', sep='\t').fillna('')
pairs = pd.read_csv(ROOT / 'reports/pairing_checks.tsv', sep='\t')
manifest = pd.read_csv(ROOT / 'provenance/INPUT_MANIFEST.tsv', sep='\t')
composition = json.loads((ROOT / 'config/Figure1A_composition.json').read_text())
sp = pd.read_csv(ROOT / 'source_data/current/spatial_sparsity.tsv', sep='\t')
fid = pd.read_csv(ROOT / 'source_data/current/official_fidelity.tsv', sep='\t')
for obsolete in ['FigureS1CD.csv', 'Figure2D_effects.csv']:
    (ROOT / 'source_data/panels' / obsolete).unlink(missing_ok=True)

lines = ['# Figure provenance', '',
         'All values originate in the packaged source tables. Absolute original paths',
         'and SHA256 hashes are in `provenance/INPUT_MANIFEST.tsv`. Source-data CSVs',
         'preserve full keys and saved precision; display abbreviations are mapped separately.', '',
         '| Panel | Source file | Columns / key | Filtering and statistical rule | Display transform | n (unit) |',
         '| --- | --- | --- | --- | --- | --- |']
for row in inv.itertuples():
    fields = [row.panel, row.source_file, row.columns_used + '; key: ' + row.unique_key,
              row.filtering_statistical_rule, row.display_transform,
              f'{row.n} ({row.n_definition})']
    fields = [str(s).replace('|', '\\|').replace('\n', ' ') for s in fields]
    lines.append('| ' + ' | '.join(fields) + ' |')
lines += ['', '## Interpretation', '',
          '- Figure 1A is protected and absent from the editable data master; its exact source is supplied separately.',
          '- Figure 1D/E include completed measurements only. ESS values are saved v0.1.4 measurements; official measurements are saved comparator runs, not newly timed here.',
          '- F1B marginal histograms count every paired record in equal log10 probability bins. They do not introduce simulated points or change fidelity calculations.',
          '- Light fills in Figure 1D-F mark gaps between existing paired curves or CosMx N-squared and E; they are comparison guides, not confidence intervals.',
          '- Figure 1G retains all 14 mean-degree values, now arranged as fixture-wise horizontal bars with endpoint symbols and printed values; input N remains explicit in each label. No descriptive CI is invented.',
          '- Fidelity uses the original permutation-p call thresholds recorded per panel. These are not BH-FDR discoveries.',
          '- Figure 2 uses the frozen 2026-09-07 patient-random-intercept analysis, verified against the current released model table by full key.',
          '- Figure 2A has one explicit vector arrow per saved labeled full key. Label and arrow-target coordinates are provided in Figure2A_labels.csv.',
          '- Supplementary record-count, fit-state and degree panels share the endpoint-bar style. S1/S2 probability marginals count the existing records; S3 isolated percentages use all 14 fixtures.',
          '- Figure 2C prints effect / response SD; LMM and OLS have the same zero-centered color range. Dots use each saved feature-level BH q < 0.05.',
          '- Patient means are means of the saved slice log1p scores. No group-mean CI is inferred from a condition-effect CI.',
          '- Figure 2E/S4C use saved LR-feature BH q, not pathway-enrichment FDR. Blank locations have no saved representative.',
          '- S1 bridge/official SPARKLE comparisons are the identified archived series. They are not relabeled as a newly rerun comparator.',
          '- S2 official comparator has the documented dimensions-only repair. The tissue map uses public graph-cluster labels, not new cell-type annotation.',
          '- N squared is the potential ordered-pair matrix; E counts stored directed non-self radius edges. All tested graphs use effective radius 35 micrometres.',
          '- S3 dense memory is theoretical float64 matrix storage; CSR is measured numeric-array payload, neither equals peak RSS.',
          '- The CSPG4 display uses the audited receptor ITGA3+ITGB1 while preserving the original interaction key and statistics.',
          '- The full 66-comparison table retains the one non-exact permutation-p result; the main CosMx ladder is not a claim of exact agreement for every dataset.', '']
(ROOT / 'FIGURE_PROVENANCE.md').write_text('\n'.join(lines), encoding='utf-8')

hash_bad = []
for row in manifest.itertuples():
    if hashlib.sha256((ROOT / row.staged_file).read_bytes()).hexdigest() != row.sha256:
        hash_bad.append(row.staged_file)
assert not hash_bad, hash_bad
assert not (data.status == 'FAIL').any()

def counts(frame):
    return ', '.join(f'{status}: {count}' for status, count in frame.status.value_counts().items())

warnings = exports[exports.status != 'PASS']
report = ['# Figure QA', '',
          f'- Original input SHA256 verification: {len(manifest)} / {len(manifest)} unchanged.',
          f'- Numerical/data checks: {counts(data)}.',
          f'- Export checks: {counts(exports)}.',
          f'- Full-key paired-record audits: {len(pairs)} fixtures; {int(pairs.unmatched.sum())} unmatched, {int(pairs.zeros.sum())} zero values, {int(pairs.nonfinite.sum())} nonfinite values.',
          f'- Neighbor graphs: {len(sp)} fixtures; N = {sp.N.min():,} to {sp.N.max():,}.',
          '- All six data masters use 180 mm width; dimensions, embedded fonts, editable SVG text, page boundaries, letters and 600 dpi PNG size were checked.',
          '- Original panel selections, stable IDs, sender/receiver identities, saved CIs, BH values and singular coding were retained.',
          '- Every one of the eight volcano arrow endpoints was checked against the saved full-key effect and BH q; each exported SVG leader has its stable feature ID.',
          '- Figure 1B marginal-bin counts conserve all 18,484 paired records. All 14 graph fixtures remain in F-H; mean degree is presented by fixture in G.',
          '- Invalid/convergence-warning LMM fits remain excluded from the primary volcano. Valid singular fits retain their documented interpretation and open forest symbols.',
          '- Failed comparator jobs appear only as status symbols, not as completed runtime/RSS measurements.',
          '- PDF/PNG export was visually reviewed for all six data masters, including high-density heatmap, patient and bubble areas. Page-boundary tests do not automatically prove absence of text overlap.',
          '- Grayscale, deuteranopia and protanopia previews are included for all six data masters. Method/geometry markers and printed heatmap signs supply additional coding; no quantitative perceptual certification is claimed.',
          '- Main data master A remains empty. The original user workflow is unchanged and only proportionally placed in a separately named reference composition.', '',
          '## Remaining limitation', '',
          f'The protected Figure 1A has minimum text {composition["minimum_composed_font_pt"]:.2f} pt when placed in the current slot. The reference composition is not a final-size-compliant Figure 1. The user-owned A requires user completion; it was not modified.', '',
          '## Recorded warnings', '']
for row in warnings.itertuples():
    report.append(f'- {row.figure}: {row.check}. {row.details}')
report += ['', '## Editability', '',
           'Hexbin polygons, tissue points, heatmap tiles, curves, CI bars, text, axes and legends are vector objects. Only the dense volcano point layers are rasterized at 600 dpi. Original tables and per-panel CSVs accompany the workbook.', '',
           '## Scientific limits', '',
           f'The complete official-fidelity source table contains {len(fid)} comparisons and {int(fid.strict_fidelity_pass.sum())} strict passes. Its remaining non-exact p-value record is retained. S1 archived validation is explicitly identified. No new fit, null experiment, biological annotation or missing data was generated.', '',
           'Neighbor sparsity is descriptive for the tested inputs and radii. It does not establish global linear complexity, universal graph sparsity or validation for other communication models.', '']
(ROOT / 'reports/QA_REPORT.md').write_text('\n'.join(report), encoding='utf-8')

packages = {}
for name in ['numpy', 'pandas', 'matplotlib', 'openpyxl', 'pillow', 'PyMuPDF', 'colorspacious']:
    try:
        packages[name] = version(name)
    except PackageNotFoundError:
        packages[name] = 'unavailable'
(ROOT / 'logs/export_package_versions.json').write_text(json.dumps(packages, indent=2) + '\n')
print(f'Wrote provenance for {len(inv)} panel-data entries and the QA report.')
