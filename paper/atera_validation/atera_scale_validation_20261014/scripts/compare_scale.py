"""Compare saved full keys and write new comparisons outside reused source results."""
import sys
import numpy as np
import pandas as pd
from scale_common import CFG, ROOT, write_json, read_row

KEY = ['sender_group_name', 'receiver_group_name', 'interaction_name']


def compare(ess_tag, ref_tag):
    ess_path, ref_path = ROOT / 'results' / ess_tag, ROOT / 'results' / ref_tag
    ess = pd.read_csv(ess_path / 'records.tsv', sep='\t', float_precision='round_trip')
    ref = pd.read_csv(ref_path / 'records.tsv', sep='\t', float_precision='round_trip')
    if ess.duplicated(KEY).any() or ref.duplicated(KEY).any():
        raise ValueError('Nonunique complete sender/receiver/LR key')
    finite = bool(np.isfinite(ess[['probability', 'pvalue']]).all().all() and
                  np.isfinite(ref[['probability', 'pvalue']]).all().all())
    ess, ref = ess[ess.probability > 0], ref[ref.probability > 0]
    joined = ess.merge(ref, on=KEY, how='outer', suffixes=('_ess', '_reference'), indicator=True)
    matched = joined[joined['_merge'] == 'both']
    dp = (matched.probability_ess - matched.probability_reference).abs()
    dq = (matched.pvalue_ess - matched.pvalue_reference).abs()
    left = read_row(ess_path / 'summary.tsv')
    right = read_row(ref_path / 'summary.tsv')
    contract = all(left[x] == right[x] for x in ['input_md5', 'N', 'genes', 'lr', 'nboot', 'seed_use'])
    n, union = len(matched), len(joined)
    decisions = int(((matched.pvalue_ess < .05) != (matched.pvalue_reference < .05)).sum())
    inclusive = int(((matched.pvalue_ess <= .05) != (matched.pvalue_reference <= .05)).sum())
    strict = bool(contract and finite and n > 0 and n == union and
                  dp.max() <= CFG['agreement_probability_tolerance'] and dq.max() == 0 and
                  decisions == 0 and inclusive == 0)
    row = dict(ess_tag=ess_tag, reference_tag=ref_tag, tissue=left['tissue'], fixture=left['fixture'],
               N=int(left['N']), genes=int(left['genes']), lr=int(left['lr']), nboot=int(left['nboot']),
               ess_active=len(ess), reference_active=len(ref), matched=n,
               ess_only=int((joined['_merge'] == 'left_only').sum()),
               reference_only=int((joined['_merge'] == 'right_only').sum()),
               active_jaccard=n / union if union else np.nan,
               max_absolute_probability_error=float(dp.max()) if n else np.nan,
               max_absolute_pvalue_error=float(dq.max()) if n else np.nan,
               exact_pvalue_fraction=float((dq == 0).mean()) if n else np.nan,
               significance_disagreements=decisions, inclusive_significance_disagreements=inclusive,
               same_input_contract=contract, finite=finite, strict_pass=strict)
    out = ROOT / 'comparisons' / ess_tag
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(out / 'agreement.tsv', sep='\t', index=False)
    mask = ((joined['_merge'] != 'both') |
            ((joined.probability_ess - joined.probability_reference).abs() > CFG['agreement_probability_tolerance']) |
            (joined.pvalue_ess != joined.pvalue_reference))
    joined[mask].to_csv(out / 'agreement_exceptions.tsv', sep='\t', index=False)
    write_json(out / 'complete.json', row)
    print(row, flush=True)
    return strict


if __name__ == '__main__':
    sys.exit(0 if compare(sys.argv[1], sys.argv[2]) else 2)
