"""Refresh descriptive outputs from successful checkpoints and measured resources."""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT/'config/design.json').read_text())
KEY=['sender_group_name','receiver_group_name','interaction_name']


def table(name, frames):
    value=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()
    value.to_csv(ROOT/'tables'/f'{name}.tsv',sep='\t',index=False)
    value.to_csv(ROOT/'tables'/f'{name}.csv',index=False)
    return value


def get_records(tag):
    return pd.read_csv(ROOT/'results'/tag/'records.tsv',sep='\t')


def stability(performance):
    rows=[]
    if performance.empty:
        return rows
    for tissue in CFG['dataset_order']:
        full_tag=f'{tissue}_full_SpatialESS_formal_b{CFG["nboot"]}'
        if not (ROOT/'results'/full_tag/'complete.json').exists():
            continue
        candidates=performance[(performance.tissue==tissue)&(performance.engine=='SpatialESS')&
            (performance.phase=='formal')&(performance.fixture!='full')]
        todo=[]
        for row in candidates.to_dict('records'):
            cached=ROOT/'results'/row['tag']/'sampling_stability.json'
            if cached.exists():
                rows.append(json.loads(cached.read_text()))
            else:
                todo.append(row)
        if not todo:
            continue
        full=get_records(full_tag)
        full=full[full.probability>0][KEY+['probability']]
        full=full.sort_values(['probability',*KEY],ascending=[False,True,True,True])
        top=full.head(100)
        for row in todo:
            sub=get_records(row['tag'])
            sub=sub[sub.probability>0][KEY+['probability']]
            groups=pd.read_csv(ROOT/'inputs'/tissue/f'{row["fixture"]}_groups.tsv',sep='\t')['group']
            eligible=full[full.sender_group_name.isin(groups)&full.receiver_group_name.isin(groups)]
            joined=eligible.merge(sub,on=KEY,how='outer',suffixes=('_full','_subset'),indicator=True)
            paired=joined[joined['_merge']=='both']
            matched=len(paired); union=len(joined)
            probs=joined[['probability_full','probability_subset']].fillna(0)
            top_matched=len(top.merge(sub,on=KEY,how='inner'))
            result=dict(tissue=tissue,fixture=row['fixture'],N=row['N'],groups=row['groups'],
              full_N=len(pd.read_csv(ROOT/'inputs'/tissue/'full_barcodes.tsv',sep='\t',usecols=['cell_id'])),
              full_active=len(full),full_active_with_groups_present=len(eligible),subset_active=len(sub),
              matched_active=matched,active_jaccard=matched/union if union else None,
              full_active_retained_fraction=matched/len(eligible) if len(eligible) else None,
              rank_spearman_matched_positive=float(paired.probability_full.corr(paired.probability_subset,method='spearman')) if matched>1 else None,
              rank_spearman_active_union_zero_inclusive=float(probs.probability_full.corr(probs.probability_subset,method='spearman')) if union>1 else None,
              top100_full_active_retained=top_matched,top_full_count=len(top),
              note='Descriptive full-tissue vs subset comparison; sender/receiver groups must be present in subset for eligible domain. Missing inactive probabilities are 0 only within this explicit domain. Geometry/composition change the estimand; full tissue is not biological ground truth.')
            (ROOT/'results'/row['tag']/'sampling_stability.json').write_text(json.dumps(result,indent=2,allow_nan=True)+'\n')
            rows.append(result)
    return rows


def main():
    performances=[]; failures=[]; agreements=[]; graphs=[]; manifests=[]; coverage=[]
    for path in sorted((ROOT/'logs').glob('*.exit.json')):
        record=json.loads(path.read_text())
        if 'engine' not in record:
            continue
        if record['returncode']!=0:
            failures.append({k:record.get(k) for k in ['tag','tissue','fixture','engine','phase','nboot','status','returncode','elapsed_seconds']})
            continue
        summary=ROOT/'results'/record['tag']/'summary.tsv'
        if not summary.exists():
            continue
        frame=pd.read_csv(summary,sep='\t')
        time_path=ROOT/'logs'/f'{record["tag"]}.time.txt'
        text=time_path.read_text() if time_path.exists() else ''
        match=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',text)
        frame['peak_RSS_KiB']=int(match[1]) if match else np.nan
        frame['peak_RSS_GB_decimal']=frame.peak_RSS_KiB*1024/1e9
        frame['whole_process_seconds']=record['elapsed_seconds']
        performances.append(frame)
        directory=summary.parent
        if (directory/'agreement.tsv').exists(): agreements.append(pd.read_csv(directory/'agreement.tsv',sep='\t'))
        if (directory/'graph_metrics.tsv').exists(): graphs.append(pd.read_csv(directory/'graph_metrics.tsv',sep='\t'))
    for tissue in CFG['dataset_order']:
        directory=ROOT/'inputs'/tissue
        if (directory/'fixture_manifest.tsv').exists(): manifests.append(pd.read_csv(directory/'fixture_manifest.tsv',sep='\t'))
        if (directory/'lr_coverage.tsv').exists():
            frame=pd.read_csv(directory/'lr_coverage.tsv',sep='\t'); frame.insert(0,'tissue',tissue)
            coverage.append(frame)
    performance=table('atera_performance',performances)
    failure=table('atera_failed_runs',[pd.DataFrame(failures)] if failures else [])
    agreement=table('atera_official_agreement',agreements)
    graph=table('atera_spatial_graph_metrics',graphs)
    manifest=table('atera_fixture_manifest',manifests)
    lr=table('atera_lr_coverage',coverage)
    sample=table('atera_sampling_stability',[pd.DataFrame(stability(performance))])
    with pd.ExcelWriter(ROOT/'tables/Atera_supplementary_results.xlsx',engine='openpyxl') as writer:
        for name,frame in [('performance',performance),('agreement',agreement),('failed_runs',failure),
                          ('graph_metrics',graph),('fixtures',manifest),('LR_coverage',lr),('sampling',sample)]:
            frame.to_excel(writer,sheet_name=name,index=False)
    formal=performance[performance.phase=='formal'] if not performance.empty else pd.DataFrame()
    completed=0 if formal.empty else len(formal)
    report=['# Atera supplementary experiment — live summary','',
      f'Completed inference runs: {len(performance)}; completed formal runs: {completed}; failed inference attempts: {len(failures)}.',
      '', 'Two independent Atera WTA preview tissues are evaluated. There are 18,028 measured Gene Expression rows per input; control/probe features are excluded. Missing official annotations and QC exclusions are listed, and exact barcode files are retained.',
      '', 'Both engines receive the same normalized full panel, spatial coordinates, official groups, mandatory-component-complete LR database rows, spatial settings, seeds and permutation counts. This is a fixed-LR computational benchmark that bypasses overexpression-based LR selection; it is not the complete standard LR-selection workflow.',
      '', 'Performance: whole-process wall time includes R startup, input loading and output writing. Inference time is separately measured inside R. Peak RSS is measured by GNU time; CSR array bytes and theoretical dense matrix bytes are separate quantities. Failed/timeout runs contribute no completed performance values or numerical-agreement claims.',
      '', 'Sampling comparisons describe stability under changed geometry/composition. The full tissue is not biological ground truth. Replicate subsets are not independent patients; these two specimens do not establish cohort-level, causal or general biological validation.',
      '', 'Tables are refreshed after each task. Raw permutation p < 0.05 is descriptive and is not pathway-enrichment FDR. No new model, significance rule or biochemical validation is introduced.']
    if not agreement.empty:
        report += ['',f'Strictly passing reference comparisons: {int(agreement.strict_pass.sum())}/{len(agreement)}. Strict pass requires full active-key agreement, maximum absolute probability difference <= 1e-15, identical p values and no decision disagreements. Exceptions are retained per job.']
    (ROOT/'reports/Atera_results_summary.md').write_text('\n'.join(report)+'\n')
    (ROOT/'reports/table_status.json').write_text(json.dumps(dict(completed_inference=len(performance),
          completed_formal=completed,failed_inference=len(failures),agreement_pairs=len(agreement)),indent=2)+'\n')
    if not formal.empty:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none','font.size':7,'axes.linewidth':.5})
        fig,axes=plt.subplots(1,2,figsize=(7.09,2.4))
        colors={'SpatialESS':'#0072B2','SpatialCellChat':'#6B7280'}
        markers={'breast':'o','cervical':'s'}
        ladder=formal[formal.fixture.str.startswith('central_')|(formal.fixture=='full')]
        for (tissue,engine),frame in ladder.groupby(['tissue','engine']):
            frame=frame.sort_values('N')
            for ax,y in zip(axes,['whole_process_seconds','peak_RSS_GB_decimal']):
                ax.plot(frame.N,frame[y],marker=markers[tissue],ms=3,lw=.9,color=colors[engine],label=f'{tissue}: {engine}')
                ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Cells'); ax.spines[['top','right']].set_visible(False)
        axes[0].set_ylabel('Whole-process wall time (s)'); axes[1].set_ylabel('Measured peak RSS (GB)')
        axes[0].legend(fontsize=6,frameon=False)
        fig.tight_layout()
        for extension in ['pdf','svg','png']:
            fig.savefig(ROOT/'figures'/f'Atera_performance_live.{extension}',dpi=400)
        plt.close(fig)
    print('TABLES refreshed',len(performance),'inference runs',completed,'formal',flush=True)


if __name__=='__main__':
    main()
