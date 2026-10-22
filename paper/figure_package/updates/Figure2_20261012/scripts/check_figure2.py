"""Check source preservation, support counts, editable exports and PDF text bounds."""
from pathlib import Path
import argparse
import json
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--original-package',type=Path,required=True)
args=parser.parse_args();original=args.original_package;checks=[]
def check(name,passed,detail=''):
    checks.append(dict(check=name,passed=bool(passed),detail=detail))
    if not passed:raise AssertionError(name+': '+detail)
for panel in ['Figure2A','Figure2A_labels','Figure2B','Figure2C','Figure2E']:
    a=pd.read_csv(original/'source_data/panels'/f'{panel}.csv',float_precision='round_trip')
    b=pd.read_csv(ROOT/'source_data/panels'/f'{panel}.csv',float_precision='round_trip')
    if panel=='Figure2A_labels':
        # Arrow positions in mm are measured from the lower edge of the taller page.
        core=[c for c in a if c not in ['arrow_start_y_mm','arrow_tip_y_mm']]
        check(panel+' numerical labels and targets unchanged',a[core].equals(b[core]))
        check('Arrow physical position retained on taller page',
              np.allclose(b.arrow_start_y_mm-a.arrow_start_y_mm,12) and
              np.allclose(b.arrow_tip_y_mm-a.arrow_tip_y_mm,12))
    else:
        check(panel+' numerical data unchanged',a.equals(b))
slices=pd.read_csv(ROOT/'source_data/panels/Figure2D_slices.csv',float_precision='round_trip')
patients=pd.read_csv(ROOT/'source_data/panels/Figure2D_patients.csv',float_precision='round_trip')
identities=pd.read_csv(ROOT/'source_data/panels/Figure2D_feature_identities.csv')
check('Example order and IDs',identities.display_id.tolist()==['I15','I06','I19'])
check('All slices and patients retained',len(slices)==135 and len(patients)==105 and
      slices.groupby('feature_id').size().eq(45).all() and patients.groupby('feature_id').size().eq(35).all())
check('Zero-inclusive response transformed once',np.allclose(np.log1p(slices.score),slices.log1p_score,rtol=1e-14,atol=0))
for _,r in identities.iterrows():
    s=slices[slices.feature_id==r.feature_id]
    check(r.display_id+' support equals model table',int((s.score>0).sum())==int(r.nonzero_units))
for (feature,patient),g in slices.groupby(['feature_id','patient_id']):
    mean=patients[(patients.feature_id==feature)&(patients.patient_id==patient)].log1p_score.iloc[0]
    check('Patient summary '+feature+' '+patient,np.isclose(mean,g.log1p_score.mean(),rtol=1e-14,atol=1e-17))
prior_slices=pd.read_csv(original/'source_data/panels/Figure2D_slices.csv',float_precision='round_trip')
key=identities.loc[identities.display_id=='I19','feature_id'].iloc[0]
a=prior_slices[prior_slices.feature_id==key].set_index('sample_id').sort_index()
b=slices[slices.feature_id==key].set_index('sample_id').sort_index()
check('I19 retained with original slice scores',np.allclose(a.log1p_score,b.log1p_score,rtol=1e-14,atol=0))
xml=ROOT/'reports/Figure2_text_bboxes.html'
subprocess.run(['pdftotext','-bbox-layout',str(ROOT/'figures/Figure2.pdf'),str(xml)],check=True)
tree=ET.parse(xml);ns={'h':'http://www.w3.org/1999/xhtml'};page=tree.find('.//h:page',ns)
width=float(page.get('width'));height=float(page.get('height'))
check('Fixed final dimensions 180x208 mm',np.allclose([width,height],np.array([180,208])*72/25.4,atol=.01))
words=tree.findall('.//h:word',ns);outside=[]
for word in words:
    box=[float(word.get(v)) for v in ['xMin','yMin','xMax','yMax']]
    if box[0]<-.5 or box[1]<-.5 or box[2]>width+.5 or box[3]>height+.5:outside.append(word.text)
check('All PDF text inside page',not outside,str(outside))
contents=' '.join(w.text or '' for w in words)
check('A-E panel letters present',all(any(w.text==letter for w in words) for letter in 'ABCDE'))
check('Singular spelled out', 'Singular' in contents and 'Sing.' not in contents)
check('All six support count labels visible',all(label in contents for label in ['2/10','30/35','5/10','34/35','1/10']))
# The separate C header and first data row must not intersect vertically.
headers=[word for word in words if word.text in ['Nonzero-slice','fraction','Patient','ICC','Singular']
         and 53*72/25.4<float(word.get('yMin'))<66*72/25.4]
check('All C annotation headers identified',len(headers)==5)
rows=[word for word in words if word.text=='I01' and 63*72/25.4<float(word.get('yMin'))<80*72/25.4]
check('C headers above data with clearance',bool(rows) and max(float(w.get('yMax')) for w in headers)<min(float(w.get('yMin')) for w in rows)-2)
font_info=subprocess.check_output(['pdffonts',str(ROOT/'figures/Figure2.pdf')],text=True)
check('Arial embedded as editable TrueType','ArialMT' in font_info and 'Arial-BoldMT' in font_info and 'yes yes yes' in font_info and 'Type 3' not in font_info)
(ROOT/'reports/Figure2_fonts.txt').write_text(font_info)
svg=ET.parse(ROOT/'figures/Figure2.svg');svg_ns={'s':'http://www.w3.org/2000/svg'}
check('SVG has editable text',len(svg.findall('.//s:text',svg_ns))>150)
leaders={n.get('id') for n in svg.iter()}
check('Eight original volcano leaders retained',all('leader_'+i in leaders for i in pd.read_csv(ROOT/'source_data/panels/Figure2A_labels.csv').display_id))
im=Image.open(ROOT/'figures/Figure2.png')
check('PNG final size at 600 dpi',np.allclose(im.size,np.array([180,208])/25.4*600,atol=2),str(im.size))
check('All earlier numerical checks pass',pd.read_csv(ROOT/'reports/data_checks.tsv',sep='\t').status.eq('PASS').all())
pd.DataFrame(checks).to_csv(ROOT/'reports/Figure2_QA.tsv',sep='\t',index=False)
(ROOT/'reports/QA_REPORT.md').write_text(
    f'# Figure 2 revision\n\n{len(checks)} checks passed. Figure 2 alone was regenerated; original manuscripts, package and other figures were not overwritten.\n\n'
    '- C: Nonzero-slice fraction, Patient ICC and Singular are placed in an independent aligned header band above the data. All 28 rows, coefficients, q values, support fractions and ICC values are unchanged. The singular marker is explained.\n'
    '- D: first two examples explicitly changed at the supervisor request to I15 FN1-(ITGAV+ITGB1) in Subpleural FBs and I06 COL1A1-(ITGA3+ITGB1) in Myofibroblasts; I19 HLA-DRA-CD4 in pDCs retained. All slices and patients, including zeros, are plotted. Patient means are mean(log1p(slice score)), not log1p(mean(score)).\n'
    '- D support: I15 Control 2/10, PF 30/35; I06 Control 5/10, PF 34/35; I19 Control 1/10, PF 34/35. Counts are of slices, not patients.\n'
    '- Existing A, B, C and E numeric plotting tables match the previous exports exactly. No inference, model, BH or molecular annotation was changed. The arrow target identities are retained.\n'
    '- Final masters: 180 x 208 mm, Arial regular/bold embedded, editable SVG text and vector axes, 600 dpi PNG. Only dense volcano data points are rasterized. No PDF text extends outside the page.\n'
    '- The rendered PDF was inspected at whole-figure and C/D detail sizes. C headings are separated from one another and data rows; all D group labels fit.\n'
    '- The existing S4d interval panel still describes the previous three examples and was not revised in this Figure-2-only request; the new standalone caption does not refer its intervals to the new examples. The model tables retain original intervals for the new features.\n'
    '- Biological interpretation is unchanged: the I06 patient HC3 sensitivity q is >0.05; replacing the example does not establish robust biological significance.\n')
print(f'PASS: {len(checks)} data/export checks')
