"""Select the requested saved features; no inference, model fitting or FDR changes."""
from pathlib import Path
import argparse
import hashlib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--score-table', type=Path, required=True)
args = parser.parse_args()
selection_path = ROOT/'source_data/archived/selected_interactions.tsv.gz'
metadata_path = ROOT/'source_data/archived/sample_metadata.tsv.gz'
selected = pd.read_csv(selection_path, sep='\t', float_precision='round_trip').set_index('display_id').loc[
    ['I15', 'I06', 'I19']].reset_index()
matrix = pd.read_csv(args.score_table, sep='\t', usecols=['sample_id']+selected.feature_id.tolist(),
                     float_precision='round_trip').set_index('sample_id')
metadata = pd.read_csv(metadata_path, sep='\t')
assert matrix.shape == (45, 3) and matrix.index.is_unique
assert set(matrix.index) == set(metadata.sample_id) and metadata.patient_id.nunique() == 35
long = matrix.reset_index().melt(id_vars='sample_id', var_name='feature_id', value_name='score').merge(
    metadata, on='sample_id', validate='many_to_one')
long['log1p_score'] = np.log1p(long.score)
assert len(long) == 135 and not long.isna().any().any() and (long.score >= 0).all()
long.to_csv(ROOT/'source_data/archived/figure2_examples_slice_scores.tsv.gz', sep='\t', index=False,
            float_format='%.17g')
patients = long.groupby(['feature_id', 'patient_id', 'condition'], sort=False).agg(
    log1p_score=('log1p_score', 'mean'), slices=('score', 'size')).reset_index()
assert len(patients) == 105
patients.to_csv(ROOT/'source_data/archived/figure2_examples_patient_scores.tsv.gz', sep='\t',
                index=False, float_format='%.17g')
counts = long.groupby(['feature_id', 'condition'], sort=False).agg(
    nonzero_slices=('score', lambda s: int((s > 0).sum())), total_slices=('score', 'size')).reset_index().merge(
    selected[['feature_id', 'display_id', 'ligand', 'receptor', 'sender_group_name', 'receiver_group_name']],
    on='feature_id', validate='many_to_one')
counts.to_csv(ROOT/'source_data/Figure2D_nonzero_counts.csv', index=False)
pd.DataFrame([dict(file=str(path), bytes=path.stat().st_size,
                   sha256=hashlib.sha256(path.read_bytes()).hexdigest())
              for path in [args.score_table, selection_path, metadata_path]]).to_csv(
    ROOT/'reports/EXAMPLE_INPUT_PROVENANCE.tsv', sep='\t', index=False)
print(counts[['display_id', 'condition', 'nonzero_slices', 'total_slices']].to_string(index=False))
