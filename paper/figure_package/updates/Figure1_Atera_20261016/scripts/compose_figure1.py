"""Place the unchanged workflow PDF into the enlarged Figure 1A slot."""
import hashlib
import json
import subprocess
import copy
import xml.etree.ElementTree as ET
from pathlib import Path

import pymupdf as fitz
import pandas as pd
from PIL import Image, ImageOps
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MM = 72 / 25.4
asset = ROOT / 'assets/Figure1A_user_original.pdf'
data_path = ROOT / 'figures/Figure1_data.pdf'
figure_path = ROOT / 'figures/Figure1.pdf'
composition = json.loads((ROOT / 'config/layout.json').read_text())['Figure1']
x, y, width, height = composition['panels']['A']
page_w, page_h = composition['size_mm']
source = fitz.open(asset)
data = fitz.open(data_path)
full = fitz.open()
page = full.new_page(width=page_w * MM, height=page_h * MM)
page.show_pdf_page(page.rect, data, 0)
slot = fitz.Rect(x * MM, y * MM, (x + width) * MM, (y + height) * MM)
page.show_pdf_page(slot, source, 0, keep_proportion=True)
full.set_metadata({
    'title': 'SpatialESS Figure 1: CosMx and Atera scaling',
    'subject': 'Figure 1A is placed from the author-supplied workflow PDF without internal edits.'
})
full.save(figure_path, garbage=4, deflate=True)

# Standalone A uses the same displayed physical dimensions and unchanged vector source.
a_doc = fitz.open()
a_page = a_doc.new_page(width=180 * MM, height=(height + 4) * MM)
a_slot = fitz.Rect(2 * MM, 2 * MM, 178 * MM, (height + 2) * MM)
a_page.show_pdf_page(a_slot, source, 0, keep_proportion=True)
a_doc.set_metadata({'title': 'Figure 1A: author-supplied SpatialESS workflow'})
a_doc.save(ROOT / 'figures/Figure1A.pdf', garbage=4, deflate=True)

for pdf_name, svg_name, png_name in [
    ('Figure1.pdf', 'Figure1.svg', 'Figure1.png'),
    ('Figure1A.pdf', 'Figure1A.svg', 'Figure1A.png')]:
    pdf = ROOT / 'figures' / pdf_name
    subprocess.run(['pdftocairo', '-svg', str(pdf), str(ROOT / 'figures' / svg_name)], check=True)
    doc = fitz.open(pdf)
    doc[0].get_pixmap(dpi=600, alpha=False).save(ROOT / 'figures' / png_name)
    doc[0].get_pixmap(dpi=180, alpha=False).save(ROOT / 'preview' / png_name)

# Keep data-panel labels as editable SVG text; place the workflow as vector paths.
svg_ns = 'http://www.w3.org/2000/svg'
ET.register_namespace('', svg_ns)
data_svg = ET.parse(ROOT/'figures/Figure1_data.svg')
a_svg = ET.parse(ROOT/'figures/Figure1A.svg')
data_root = data_svg.getroot()
a_root = a_svg.getroot()
workflow_layer = ET.Element(f'{{{svg_ns}}}svg', {
    'x': '0', 'y': '0', 'width': data_root.attrib['viewBox'].split()[2],
    'height': a_root.attrib['viewBox'].split()[3],
    'viewBox': a_root.attrib['viewBox'], 'preserveAspectRatio': 'none',
    'aria-label': 'Original author-supplied Figure 1A workflow, scaled without distortion'})
for child in list(a_root):
    workflow_layer.append(copy.deepcopy(child))
data_root.append(workflow_layer)
data_svg.write(ROOT/'figures/Figure1.svg', encoding='utf-8', xml_declaration=True)

# Higher-resolution, print-size A and full-figure previews.
for name in ['Figure1.pdf', 'Figure1A.pdf']:
    doc = fitz.open(ROOT / 'figures' / name)
    target = ROOT / 'preview' / name.replace('.pdf', '_600dpi.png')
    doc[0].get_pixmap(dpi=600, alpha=False).save(target)

status_panel = fitz.open(ROOT / 'panels/Figure1I.pdf')
status_panel[0].get_pixmap(dpi=300, alpha=False).save(
    ROOT / 'preview/Figure1I_detail.png')
for panel in ['D', 'E']:
    panel_doc = fitz.open(ROOT / f'panels/Figure1{panel}.pdf')
    panel_page = panel_doc[0]
    detail_clip = fitz.Rect(0, 0, panel_page.rect.width,
                            panel_page.rect.height - 1 * MM)
    panel_page.get_pixmap(dpi=300, alpha=False, clip=detail_clip).save(
        ROOT / f'preview/Figure1{panel}_detail.png')
full_preview = fitz.open(figure_path)
full_preview[0].get_pixmap(dpi=300, alpha=False,
    clip=fitz.Rect(0, 124*MM, 180*MM, 171*MM)).save(
        ROOT/'preview/Figure1_DE_detail.png')
ImageOps.grayscale(Image.open(ROOT/'preview/Figure1_DE_detail.png')).save(
    ROOT/'preview/Figure1_DE_grayscale.png')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

asset_hash = sha(asset)
input_pdf = fitz.open(asset)
source_span_sizes = [span['size'] for block in input_pdf[0].get_text('dict')['blocks']
                     if 'lines' in block for line in block['lines'] for span in line['spans']]
minimum_display_font = min(source_span_sizes) * width * MM / input_pdf[0].rect.width
comp = {
    'user_source': 'assets/Figure1A_user_original.pdf',
    'sha256': asset_hash,
    'source_pdf_sha256_matches_workflow_comparison': asset_hash == hashlib.sha256(
        Path('/home/dzf/cellchat_acceleration/1/SpatialESS_workflow_comparison.pdf').read_bytes()
    ).hexdigest(),
    'destination_mm': [x, y, width, height],
    'keep_proportion': True,
    'internal_modification': False,
    'minimum_displayed_source_font_pt_approx': minimum_display_font,
    'figure1_canvas_mm': [page_w, page_h],
    'standalone_figure1A_canvas_mm': [180, height + 4],
    'deliverables': ['figures/Figure1.pdf', 'figures/Figure1.svg', 'figures/Figure1.png',
                     'figures/Figure1A.pdf', 'figures/Figure1A.svg', 'figures/Figure1A.png']
}
(ROOT / 'config/Figure1A_composition.json').write_text(json.dumps(comp, indent=2) + '\n')

data_doc = fitz.open(data_path)
text = data_doc[0].get_text()
checks = [
    ('Author A asset unchanged', sha(asset) == asset_hash),
    ('Atera performance embedded in Figure 1', 'Atera' in text),
    ('Atera 50k comparison notes present', 'Atera 50k' in text),
    ('D/E use the original method and cohort legends',
     all(label in text for label in ['SpatialESS','SpatialCellChat','CosMx CancerousLiver',
                                    'Atera breast','Atera cervical'])),
    ('D/E contain the requested original blue/orange method identity',
     'SpatialESS' in fitz.open(ROOT/'panels/Figure1D.pdf')[0].get_text() and
     'SpatialCellChat' in fitz.open(ROOT/'panels/Figure1E.pdf')[0].get_text()),
    ('F-H shared spatial-series and all-pair legend labels present',
     all(label in text for label in ['CosMx', 'Xenium central', 'Xenium stratified',
                                    'Visium HD ovarian', 'N²: all possible ordered pairs'])),
    ('Figure 1 has one page', len(fitz.open(figure_path)) == 1),
    ('Figure 1A has one page', len(fitz.open(ROOT / 'figures/Figure1A.pdf')) == 1),
    ('Panel I uses requested equally spaced labels with 200k followed by 700k',
     all(label in status_panel[0].get_text().split() for label in
         ['1k', '5k', '10k', '25k', '50k', '100k', '200k', '700k']) and
     not any(label in status_panel[0].get_text().split() for label in ['400k','450k']) and
     status_panel[0].get_text().find('200k') < status_panel[0].get_text().find('700k')),
    ('Panel I retains source rows and hides four requested categories',
     len(pd.read_csv(ROOT/'source_data/panels/Figure1I_completion.csv')) == 50 and
     (~pd.read_csv(ROOT/'source_data/panels/Figure1I_completion.csv').displayed_in_panel).sum() == 8 and
     set(pd.read_csv(ROOT/'source_data/panels/Figure1I_completion.csv').query(
         '~displayed_in_panel').N.unique()) == {250,169940,400000,443515}),
    ('Panel I hides 169,940 category while D/E retain the completed measurement',
     not pd.read_csv(ROOT/'source_data/panels/Figure1I_completion.csv').query(
         'N == 169940').displayed_in_panel.any() and
     pd.read_csv(ROOT/'source_data/panels/Figure1D.csv').query(
         'N == 169940').status.eq('completed').any()),
    ('Panel I cervical 200k and 400k show terminal SpatialESS run states',
     len(pd.read_csv(ROOT/'source_data/panels/Figure1I_completion.csv').query(
         "tissue == 'cervical' and method == 'SpatialESS' and N in [200000,400000]")) == 2 and
     not pd.read_csv(ROOT/'source_data/panels/Figure1I_completion.csv').query(
         "tissue == 'cervical' and method == 'SpatialESS' and N in [200000,400000]")
         .status.isin(['not_attempted','running','in_progress']).any()),
    ('Panel I marks 700k reference failed while retaining the allocation failure',
     len(pd.read_csv(ROOT/'source_data/panels/Figure1I_completion.csv').query(
         "tissue == 'cervical' and method == 'SpatialCellChat' and N == 715512 and "
         "display_label == 'reference failed' and status == 'allocation_failure_under_memory_budget'")) == 1),
    ('Panel I has no explanatory input-size axis label',
     'Input size' not in status_panel[0].get_text()),
    ('Panel I uses the reference-failed legend without obsolete status labels',
     'reference failed' in status_panel[0].get_text() and
     'not reference' not in status_panel[0].get_text() and
     'not attempted' not in status_panel[0].get_text()),
    ('Panel I nominal categories have equal x spacing',
     np.array_equal(np.sort(pd.read_csv(
         ROOT/'source_data/panels/Figure1I_completion.csv').query(
         'displayed_in_panel').axis_position.unique()),
         np.arange(pd.read_csv(ROOT/'source_data/panels/Figure1I_completion.csv').query(
             'displayed_in_panel').axis_position.nunique()))),
    ('Atera full-scale N retained in I source data',
     set(pd.read_csv(ROOT/'source_data/panels/Figure1I_completion.csv').query(
         "dataset == 'Atera' and method == 'SpatialESS'").N) >= {169940,200000,400000,715512}),
]
main_svg = ET.parse(ROOT/'figures/Figure1.svg').getroot()
main_text = main_svg.findall(f'.//{{{svg_ns}}}text')
main_paths = main_svg.findall(f'.//{{{svg_ns}}}path')
main_images = main_svg.findall(f'.//{{{svg_ns}}}image')
checks.extend([
    ('Composed SVG retains editable data-panel text', len(main_text) >= 100),
    ('Composed SVG keeps A and data as vector paths', len(main_paths) >= 100),
    ('Composed SVG contains no raster replacement', len(main_images) == 0),
    ('Atera D/E panel source has all 46 observations/states',
     len(pd.read_csv(ROOT/'source_data/panels/Figure1D.csv')) == 46 and
     len(pd.read_csv(ROOT/'source_data/panels/Figure1E.csv')) == 46),
])
full_doc = fitz.open(figure_path)
full_page = full_doc[0]
figure_spans = [span for block in full_page.get_text('dict')['blocks'] if 'lines' in block
                for line in block['lines'] for span in line['spans']]
full_fonts = full_page.get_fonts(full=True)
main_png = Image.open(ROOT/'figures/Figure1.png')
a_png = Image.open(ROOT/'figures/Figure1A.png')
checks.extend([
    ('Figure 1 canvas is 180 x 250 mm',
     np.allclose([full_page.rect.width/MM, full_page.rect.height/MM], [180,250], atol=.01)),
    ('Figure 1 text stays inside page',
     all(full_page.rect.contains(fitz.Rect(span['bbox'])) for span in figure_spans)),
    ('Figure 1 PDF fonts are embedded and non-Type-3',
     all(font[1] != 'n/a' and font[2] != 'Type3' for font in full_fonts)),
    ('Figure 1 full-size PNG is 600 dpi', main_png.size == (4252,5906)),
    ('Figure 1A full-width PNG is 600 dpi', a_png.size == (4252,1950)),
    ('Atera paired table retains all ten strict passes',
     len(pd.read_csv(ROOT/'source_data/panels/Figure1C_Atera_agreement.csv')) == 10 and
     pd.read_csv(ROOT/'source_data/panels/Figure1C_Atera_agreement.csv').strict_pass.astype(bool).all()),
])
pd.DataFrame(checks, columns=['check', 'passed']).to_csv(
    ROOT / 'reports/composition_checks.tsv', sep='\t', index=False)
pd.DataFrame(checks, columns=['check','status']).assign(
    status=lambda d: d.status.map({True:'PASS',False:'FAIL'})).to_csv(
        ROOT/'reports/export_checks.tsv',sep='\t',index=False)

# Confirm frozen inputs were copied without modification; the added Atera sources are new.
baseline = ROOT.parent/'SpatialESS_Figure1_revised_20261012'
source_checks=[]
for folder in ['source_data/current','source_data/archived']:
    for path in sorted((baseline/folder).rglob('*')):
        if not path.is_file():
            continue
        rel=path.relative_to(baseline)
        staged=ROOT/rel
        same=staged.is_file() and hashlib.sha256(path.read_bytes()).digest()==hashlib.sha256(staged.read_bytes()).digest()
        source_checks.append(dict(file=str(rel),status='PASS' if same else 'FAIL'))
source_checks.append(dict(file='assets/Figure1A_user_original.pdf',
    status='PASS' if asset_hash==hashlib.sha256(
        (baseline/'assets/Figure1A_user_original.pdf').read_bytes()).hexdigest() else 'FAIL'))
pd.DataFrame(source_checks).to_csv(ROOT/'reports/source_data_checks.tsv',sep='\t',index=False)
if not all(ok for _, ok in checks):
    raise SystemExit('Composition validation failed; see reports/composition_checks.tsv')
if not all(row['status']=='PASS' for row in source_checks):
    raise SystemExit('A frozen source input changed; see reports/source_data_checks.tsv')
print(pd.DataFrame(checks, columns=['check', 'passed']).to_string(index=False))
