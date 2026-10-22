"""Pair exact sender/receiver/LR keys without epsilon or silent missing records."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

KEY = ['sender_group_name', 'receiver_group_name', 'interaction_name']


def compare(root, ess_tag, reference_tag):
    cfg = json.loads((root/'config/design.json').read_text())
    ess_dir, ref_dir = root/'results'/ess_tag, root/'results'/reference_tag
    ess = pd.read_csv(ess_dir/'records.tsv', sep='\t')
    ref = pd.read_csv(ref_dir/'records.tsv', sep='\t')
    if ess.duplicated(KEY).any() or ref.duplicated(KEY).any():
        raise ValueError('Nonunique full feature key')
    ess = ess[ess.probability>0].copy()
    ref = ref[ref.probability>0].copy()
    joined = ess.merge(ref,on=KEY,how='outer',suffixes=('_ess','_reference'),indicator=True)
    paired = joined[joined['_merge']=='both']
    union, matched = len(joined),len(paired)
    dprob = np.abs(paired.probability_ess-paired.probability_reference)
    dp = np.abs(paired.pvalue_ess-paired.pvalue_reference)
    finite = bool(np.isfinite(paired[['probability_ess','probability_reference','pvalue_ess','pvalue_reference']]).all().all())
    summary = pd.read_csv(ess_dir/'summary.tsv',sep='\t').iloc[0].to_dict()
    ref_summary = pd.read_csv(ref_dir/'summary.tsv',sep='\t').iloc[0]
    same_contract = (summary['input_md5']==ref_summary['input_md5'] and
                     summary['nboot']==ref_summary['nboot'] and summary['lr']==ref_summary['lr'] and
                     summary['N']==ref_summary['N'] and summary['genes']==ref_summary['genes'])
    jaccard = matched/union if union else np.nan
    maxprob = float(dprob.max()) if matched else np.nan
    maxp = float(dp.max()) if matched else np.nan
    p_fraction = float((dp==0).mean()) if matched else np.nan
    decisions = int(((paired.pvalue_ess<.05)!=(paired.pvalue_reference<.05)).sum())
    # Empty comparisons are not evidence of numerical agreement.
    strict = bool(same_contract and finite and matched>0 and matched==union and
                  maxprob<=cfg['agreement_probability_tolerance'] and maxp==0 and decisions==0)
    row = dict(ess_tag=ess_tag,reference_tag=reference_tag,tissue=summary['tissue'],
      fixture=summary['fixture'],phase=summary['phase'],N=summary['N'],genes=summary['genes'],
      lr=summary['lr'],nboot=summary['nboot'],ess_active=len(ess),reference_active=len(ref),
      matched=matched,ess_only=int((joined['_merge']=='left_only').sum()),
      reference_only=int((joined['_merge']=='right_only').sum()),active_jaccard=jaccard,
      max_absolute_probability_error=maxprob,max_absolute_pvalue_error=maxp,
      exact_pvalue_fraction=p_fraction,significance_disagreements=decisions,
      probability_pearson_raw=float(paired.probability_ess.corr(paired.probability_reference)) if matched>1 else np.nan,
      same_input_contract=same_contract,finite=finite,strict_pass=strict)
    pd.DataFrame([row]).to_csv(ess_dir/'agreement.tsv',sep='\t',index=False)
    mismatch = joined[(joined['_merge']!='both') |
       ((joined.probability_ess-joined.probability_reference).abs()>cfg['agreement_probability_tolerance']) |
       (joined.pvalue_ess!=joined.pvalue_reference)]
    mismatch.to_csv(ess_dir/'agreement_exceptions.tsv',sep='\t',index=False)
    print(json.dumps(row,default=str),flush=True)
    return strict


if __name__ == '__main__':
    passed = compare(Path(sys.argv[1]),sys.argv[2],sys.argv[3])
    sys.exit(0 if passed else 2)
