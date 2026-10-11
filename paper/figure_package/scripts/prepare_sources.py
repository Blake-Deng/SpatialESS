"""Stage immutable existing result tables; no inference or model fitting."""
from pathlib import Path
import argparse
import hashlib
import shutil
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--project', type=Path, default=ROOT.parents[1])
args = parser.parse_args()
project = args.project
old = project / '1/SpatialESS_figures_bluered_editable_20260923'
release = project / 'releases/SpatialESS_v0.1.4_20261001/repositories/SpatialESS'
manifest = []

def copy(source, target):
    dest = ROOT / target
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, dest)
    manifest.append(dict(source=str(source), staged_file=str(dest.relative_to(ROOT)),
                         sha256=hashlib.sha256(dest.read_bytes()).hexdigest()))

for path in (old / 'source_data').glob('*.tsv.gz'):
    copy(path, 'source_data/archived/' + path.name)
for path in (release / 'results/tables').glob('*.tsv'):
    copy(path, 'source_data/current/' + path.name)
for name in ['CSPG4_DISPLAY_MAPPING.tsv', 'CSPG4_SOURCE_AUDIT.md', 'FIGURE2D_IDENTITIES.tsv']:
    copy(old / 'audit' / name, 'source_data/annotation/' + name)
for name in ['redraw_figures.R', 'figure2_rendered_source.R', 'supplement_rendered_source.R']:
    copy(old / 'scripts' / name, 'provenance/original_scripts/' + name)
for name in ['README.md', 'docs/RESULTS.md', 'docs/MULTISAMPLE_LMM.md']:
    copy(release / name, 'provenance/current_docs/' + Path(name).name)
for name in ['MAIN_VALIDATION.md', 'LMM_RESULTS_20260907.md', 'REVIEWER_GUIDE.md']:
    copy(project / 'SpatialESS/docs' / name, 'provenance/historical_docs/' + name)
copy(project / '1/SpatialESS_workflow_comparison.pdf', 'assets/Figure1A_user_original.pdf')
copy(project / 'releases/SpatialESS_v0.1.4_20261001/build_inputs.tsv',
     'provenance/build_inputs.tsv')
copy(release / 'results/lmm/patient_random_intercept_lmm.tsv.gz',
     'source_data/current/patient_random_intercept_lmm.tsv.gz')
for sub in ['figures', 'panels', 'source_data/panels', 'reports', 'logs', 'preview', 'config']:
    (ROOT / sub).mkdir(parents=True, exist_ok=True)
pd.DataFrame(manifest).to_csv(ROOT / 'provenance/INPUT_MANIFEST.tsv', sep='\t', index=False)
print(f'Staged {len(manifest)} original files with SHA256 checksums.')
