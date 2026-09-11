from pathlib import Path
import anndata as ad
import pandas as pd
import json

root = Path(__file__).resolve().parent
source = Path('/home/dzf/cellchat_acceleration/SpatialESS_RNA_leakage_20260902')
h5root = source / 'software/SPARKLE/evaluation/reports/h5ad_ovarian_annotated'
meta = pd.read_csv(source/'inputs/spatialess_rctd/cell_metadata.tsv', sep='\t', index_col='cell')
out = root/'inputs'; out.mkdir(exist_ok=True)
rows = []
for condition, tag in [('raw','raw'), ('sparkle','SPARKLE')]:
    src = h5root/f'ovarian_x1000-1800_y300-1100_{tag}.h5ad'
    a = ad.read_h5ad(src, backed='r')
    assert list(a.obs_names) == list(meta.index)
    assert list(a.obs['annotation'].astype(str)) == list(meta['group'])
    for n in [1000,5000,16247]:
        path = out/f'{condition}_{n}.h5ad'
        if n < a.n_obs:
            v = a[:n,:].to_memory()
            # Store the original values/dtypes and complete genes; no new normalization.
            v.write_h5ad(path, compression='gzip')
        else:
            path = src
        m = meta.iloc[:n].copy(); m.index.name='cell'
        mp = out/f'{condition}_{n}_metadata.tsv'; m.to_csv(mp, sep='\t')
        rows.append(dict(condition=condition,cells=n,h5ad=str(path),metadata=str(mp),source=str(src)))
    a.file.close()
(out/'cases.json').write_text(json.dumps(rows,indent=2)+'\n')
print('Prepared six nested RAW/corrected H5AD cases, all genes and fixed shared RCTD labels retained.')
