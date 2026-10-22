"""Archive the revised Figure 2, exact sources and reproducible plotting scripts."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
import pandas as pd

root = Path(__file__).resolve().parents[1]
assert pd.read_csv(root/'reports/Figure2_QA.tsv', sep='\t').passed.all()
files = []
for folder in ['scripts', 'config', 'figures', 'panels', 'reports', 'logs']:
    files.extend(p for p in (root/folder).rglob('*') if p.is_file()
                 and '__pycache__' not in p.parts)
files.extend([root/'README.md', root/'FIGURE2_CAPTION.md'])
files.extend(p for p in (root/'source_data').glob('*') if p.is_file())
files.extend((root/'source_data/panels').glob('Figure2*.csv'))
for name in ['cosmx_scaling.tsv', 'official_fidelity.tsv', 'spatial_sparsity.tsv',
             'cosmx_50000_paired.tsv', 'patient_random_intercept_lmm.tsv.gz']:
    files.append(root/'source_data/current'/name)
for name in ['lmm', 'selected_interactions', 'effect_heatmap', 'forest',
             'patient_examples', 'slice_scores', 'patient_scores', 'volcano_labels',
             'pathway_pair_representatives', 'sample_metadata',
             'figure2_examples_slice_scores', 'figure2_examples_patient_scores']:
    files.append(root/'source_data/archived'/f'{name}.tsv.gz')
files = sorted(set(files))
manifest = pd.DataFrame([dict(file=str(p.relative_to(root)), bytes=p.stat().st_size,
                              sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files])
manifest.to_csv(root/'FILE_SHA256.tsv', sep='\t', index=False)
files.append(root/'FILE_SHA256.tsv')
archive = root/'Figure2_revised_20261012.zip'
with ZipFile(archive, 'w', ZIP_DEFLATED, compresslevel=6) as z:
    for p in files:
        z.write(p, 'Figure2_revised/'+str(p.relative_to(root)))
with ZipFile(archive) as z:
    assert z.testzip() is None
print(archive, f'{archive.stat().st_size/1e6:.2f} MB; {len(files)} files; ZIP verified')
