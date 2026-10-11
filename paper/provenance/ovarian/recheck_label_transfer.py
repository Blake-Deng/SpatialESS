from pathlib import Path
import pandas as pd
root = Path(__file__).resolve().parent / 'rctd'
labels = pd.read_csv(root / 'RAW_first_type.csv')
meta = pd.read_csv(root / 'cell_metadata.tsv', sep='\t')
assert not labels.cell_id.duplicated().any()
mapped = meta.cell_id.map(labels.set_index('cell_id').first_type)
assert (meta.group == mapped.fillna('Unknown')).all()
assert mapped.notna().sum() == 16213 and mapped.isna().sum() == 34
print('PASS: 16,213 mapped labels; 34 Unknown; actual saved cell_id transfer')
