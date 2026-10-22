"""Read official H5 in blocks; export exact measured rows as sparse binary arrays."""
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / 'config/design.json').read_text())


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        while chunk := handle.read(8 * 1024**2):
            digest.update(chunk)
    return digest.hexdigest()


def main(tissue):
    raw = Path(CONFIG['raw_root']) / tissue
    out = ROOT / 'inputs' / tissue
    out.mkdir(exist_ok=True)
    if (out / 'matrix_complete.json').exists():
        print('RESUME matrix', tissue, flush=True)
        return
    validation = json.loads((Path(CONFIG['raw_root']) / 'input_validation.json').read_text())
    previous = next(x for x in validation if x['dataset'] == tissue)
    if previous['status'] == 'FAIL':
        raise ValueError('Raw input validation failed')
    cells = pd.read_csv(raw / 'cells.csv.gz', dtype={'cell_id': str})
    labels = pd.read_csv(next(raw.glob('*_cell_groups.csv')), dtype={'cell_id': str})
    if labels.cell_id.duplicated().any():
        raise ValueError('Ambiguous multiple official annotations')
    labels = labels.set_index('cell_id')
    with h5py.File(raw / 'cell_feature_matrix.h5', 'r') as handle:
        matrix = handle['matrix']
        barcodes = matrix['barcodes'].asstr()[:]
        names = matrix['features/name'].asstr()[:]
        types = matrix['features/feature_type'].asstr()[:]
        feature_ids = matrix['features/id'].asstr()[:]
        genes = np.flatnonzero(types == 'Gene Expression')
        gene_names = names[genes]
        if len(genes) != 18028 or len(set(gene_names)) != len(genes):
            raise ValueError('Unexpected measured gene rows')
        # Every exported row has one explicit gene-expression feature identity.
        gene_mapping = np.full(len(names), -1, dtype=np.int32)
        gene_mapping[genes] = np.arange(len(genes), dtype=np.int32)
        cells = cells.set_index('cell_id').loc[barcodes].reset_index()
        cells['group'] = labels['group'].reindex(barcodes).to_numpy()
        pointer = matrix['indptr'][:]
        total_counts = np.zeros(len(barcodes), dtype=np.float64)
        features = np.zeros(len(barcodes), dtype=np.int32)
        print('AUDIT', tissue, len(barcodes), 'cells', len(genes), 'genes', flush=True)
        for lo in range(0, len(barcodes), 10000):
            hi = min(len(barcodes), lo + 10000)
            begin, end = int(pointer[lo]), int(pointer[hi])
            indices = matrix['indices'][begin:end]
            values = matrix['data'][begin:end]
            lens = np.diff(pointer[lo:hi+1])
            columns = np.repeat(np.arange(hi-lo), lens)
            keep = (gene_mapping[indices] >= 0) & (values > 0)
            total_counts[lo:hi] = np.bincount(columns[keep], weights=values[keep], minlength=hi-lo)
            features[lo:hi] = np.bincount(columns[keep], minlength=hi-lo)
        cells['measured_gene_counts'] = total_counts.astype(np.int64)
        cells['measured_gene_features'] = features
        reasons = np.where(cells['group'].isna(), 'missing_official_annotation', '')
        reasons = np.where(total_counts < CONFIG['qc_min_counts'],
                           np.char.add(reasons, ';counts_below_threshold'), reasons)
        reasons = np.where(features < CONFIG['qc_min_features'],
                           np.char.add(reasons, ';features_below_threshold'), reasons)
        reasons = np.char.lstrip(reasons, ';')
        eligible = reasons == ''
        if eligible.sum() < max(CONFIG['central_sizes']):
            raise ValueError('Too few eligible cells')
        cells['exclusion_reason'] = reasons
        cells.loc[~eligible].to_csv(out / 'excluded_cells.tsv', sep='\t', index=False)
        metadata = cells.loc[eligible].copy()
        metadata['original_matrix_index_0based'] = np.flatnonzero(eligible)
        metadata.to_csv(out / 'cells.tsv', sep='\t', index=False)
        pd.DataFrame({'feature_id': feature_ids[genes], 'gene': gene_names}).to_csv(
            out / 'genes.tsv', sep='\t', index=False)
        coordinates = metadata[['x_centroid', 'y_centroid']].to_numpy()
        tree = cKDTree(coordinates)
        nearest, _ = tree.query(coordinates, k=2, workers=1)
        minimum_distance = float(nearest[:, 1].min())
        if not np.isfinite(minimum_distance) or minimum_distance <= 0:
            raise ValueError('Duplicate/invalid centroids; resolve rather than perturb coordinates')
        graph_degrees = tree.query_ball_point(coordinates,
                         CONFIG['parameters']['interaction.range'] + CONFIG['parameters']['tol'],
                         return_length=True, workers=1) - 1
        graph_edges = int(graph_degrees.sum())
        if graph_edges > CONFIG['max_edges']:
            raise ValueError(f'Graph exceeds configured CSR bound: E={graph_edges}')
        # A fixed per-tissue scale >= 2/minimum distance also satisfies the reference guard.
        distance_scale = float(2 / minimum_distance)
        print('QC', tissue, 'retained', int(eligible.sum()), 'excluded', int((~eligible).sum()),
              'E', graph_edges, 'distance_scale', distance_scale, flush=True)
        del tree, nearest
        output_pointer = [0]
        global_max = 0.0
        with (out / 'x.bin').open('wb') as xfile, (out / 'i.bin').open('wb') as ifile:
            for lo in range(0, len(barcodes), 10000):
                hi = min(len(barcodes), lo + 10000)
                begin, end = int(pointer[lo]), int(pointer[hi])
                indices = matrix['indices'][begin:end]
                values = matrix['data'][begin:end]
                columns = np.repeat(np.arange(lo, hi), np.diff(pointer[lo:hi+1]))
                keep = (gene_mapping[indices] >= 0) & (values > 0) & eligible[columns]
                columns = columns[keep]
                values = np.log1p(values[keep].astype(np.float64) *
                         CONFIG['normalization_scale'] / total_counts[columns])
                normalized_indices = gene_mapping[indices[keep]]
                global_max = max(global_max, float(values.max(initial=0)))
                values.astype('<f8').tofile(xfile)
                normalized_indices.astype('<i4').tofile(ifile)
                counts = np.bincount(columns-lo, minlength=hi-lo)[eligible[lo:hi]]
                cumulative = np.cumsum(counts, dtype=np.int64) + output_pointer[-1]
                output_pointer.extend(cumulative.tolist())
                print('EXPORT', tissue, hi, '/', len(barcodes), 'nnz', output_pointer[-1], flush=True)
        nnz = output_pointer[-1]
        if nnz > np.iinfo(np.int32).max:
            raise ValueError('Expression itself exceeds dgCMatrix limit; do not cast indices')
        np.asarray(output_pointer, dtype='<i4').tofile(out / 'p.bin')
        metadata['group'].value_counts().rename_axis('group').reset_index(name='cells').to_csv(
            out / 'group_counts.tsv', sep='\t', index=False)
        contract = {
            'tissue': tissue, 'raw_cells': len(barcodes), 'eligible_cells': int(eligible.sum()),
            'genes': len(genes), 'nnz': nnz, 'normalization': 'log1p(10000 * count / measured-GEX library size)',
            'normalized_expression_max': global_max, 'qc': {
                'min_counts': CONFIG['qc_min_counts'], 'min_features': CONFIG['qc_min_features'],
                'missing_annotation_excluded': int(cells['group'].isna().sum())},
            'minimum_nonzero_distance_um': minimum_distance,
            'scale_distance': distance_scale, 'full_graph_edges_ckdtree_diagnostic': graph_edges,
            'self_edges_in_csr': False, 'self_scores': 'preserved by SpatialESS reference convention',
            'binary_endianness': 'little', 'x_dtype': 'float64', 'i_p_dtype': 'int32',
            'source': str(raw), 'official_experiment': json.loads((raw/'experiment.xenium').read_text()),
            'hashes': {p.name: sha256(p) for p in out.iterdir() if p.is_file()},
            'environment': {'python': sys.version, 'numpy': np.__version__, 'h5py': h5py.__version__,
                            'pandas': pd.__version__}}
        (out / 'matrix_complete.json').write_text(json.dumps(contract, indent=2) + '\n')
        print('MATRIX COMPLETE', tissue, json.dumps({k:contract[k] for k in ['eligible_cells','genes','nnz']}), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
