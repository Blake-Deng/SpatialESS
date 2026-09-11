"""Freeze nested center and submitted stratified cervical subsets."""
import hashlib
import json
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import io, sparse
from scipy.spatial import cKDTree

root = Path(__file__).resolve().parent
raw = Path('/home/dzf/cellchat_acceleration/external/cervical_xenium_20260910/raw')
submitted = Path('/home/dzf/cellchat_acceleration/validation/collaborator_cervical_20260911/inputs')
prefix = 'Xenium_Prime_Cervical_Cancer_FFPE'
out = root / 'validation/dual_track'
out.mkdir(parents=True, exist_ok=True)
with h5py.File(raw / f'{prefix}_cell_feature_matrix.h5') as h:
    m = h['matrix']
    barcodes = pd.Index(m['barcodes'].asstr()[:])
    counts = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]),
                              shape=tuple(m['shape'][:]))
    mask = m['features/feature_type'].asstr()[:] == 'Gene Expression'
    genes = m['features/name'].asstr()[:][mask]
    counts = counts[mask, :].tocsc()
meta = pd.read_csv(raw / f'{prefix}_cells.csv.gz').set_index('cell_id').loc[barcodes]
with tarfile.open(raw / f'{prefix}_analysis.tar.gz') as t:
    with t.extractfile('analysis/clustering/gene_expression_graphclust/clusters.csv') as stream:
        cluster = pd.read_csv(stream).set_index('Barcode').Cluster.reindex(barcodes)
xy = meta[['x_centroid', 'y_centroid']].to_numpy()
eligible = ((meta.transcript_counts.to_numpy() >= 20)
            & (np.asarray(counts.sum(axis=0)).ravel() > 0)
            & np.isfinite(xy).all(axis=1) & cluster.notna().to_numpy())
eligible_idx = np.flatnonzero(eligible)
median = np.median(xy[eligible_idx], axis=0)
center_idx = eligible_idx[np.argmin(np.sum((xy[eligible_idx]-median)**2, axis=1))]
ranked = eligible_idx[np.argsort(np.sum((xy[eligible_idx]-xy[center_idx])**2, axis=1), kind='stable')]
pd.DataFrame({'gene': genes}).to_csv(out / 'genes.tsv', sep='\t', index=False)
summaries = []
for track in ('center', 'stratified'):
    previous = set()
    for n in (1000, 5000):
        if track == 'center':
            idx = np.sort(ranked[:n])
        else:
            ids = pd.read_csv(submitted / f'n{n}.tsv', sep='\t').cell_id
            idx = barcodes.get_indexer(ids)
            assert (idx >= 0).all() and ids.is_unique
        assert eligible[idx].all() and previous.issubset(set(idx))
        previous = set(idx)
        d = out / f'{track}_{n}'
        d.mkdir(exist_ok=True)
        cells = meta.iloc[idx].copy()
        cells['group'] = [f'Cluster_{int(k):02d}' for k in cluster.iloc[idx]]
        cells['source_column_0based'] = idx
        cells.index.name = 'cell_id'
        cells.reset_index().to_csv(d / 'metadata.tsv', sep='\t', index=False)
        io.mmwrite(d / 'counts.mtx', counts[:, idx], field='integer')
        edges = cKDTree(xy[idx]).query_pairs(35, output_type='ndarray')
        lengths = np.linalg.norm(xy[idx][edges[:, 0]]-xy[idx][edges[:, 1]], axis=1)
        edges, lengths = edges[lengths > 0], lengths[lengths > 0]
        summary = dict(case=f'{track}_{n}', cells=n, groups=int(cells.group.nunique()),
            edges_directed=2*len(edges), isolated_cells=n-len(np.unique(edges)),
            original_inferred_dimension=int(edges.max()+1),
            min_distance_um=float(lengths.min()), nested=True,
            order='original matrix index' if track == 'center' else 'submitted barcode order')
        summaries.append(summary)
        print(summary, flush=True)
minimum = min(x['min_distance_um'] for x in summaries)
scale = float(max(1, np.ceil(1.1/minimum)))
manifest = dict(dataset='10x Xenium Prime FFPE Human Cervical Cancer',
    dataset_url='https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-cervical-cancer',
    barcode_source='/home/dzf/cellchat_acceleration/benchmark_cells.rds',
    barcode_sha256=hashlib.sha256(Path('/home/dzf/cellchat_acceleration/benchmark_cells.rds').read_bytes()).hexdigest(),
    annotation_boundary='Public 23 graph clusters; collaborator biological labels not supplied',
    stratification_boundary='Submitted stratified barcode sets; biological proportionality not independently audited',
    center_cell=str(barcodes[center_idx]), center_um=xy[center_idx].tolist(),
    scale_distance=scale, scale_policy='One shared scale, ceil(1.1/minimum nonzero distance across all cases)',
    parameters=dict(ratio=1, tol=5, interaction_range=30, contact_range=10, nboot=100, seed=1),
    cases=summaries)
(out / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
pd.DataFrame(summaries).to_csv(out / 'sampling.tsv', sep='\t', index=False)
