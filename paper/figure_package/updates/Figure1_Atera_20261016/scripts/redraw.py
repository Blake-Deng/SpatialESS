"""Redraw frozen data. This script does not fit models or recompute FDR."""
import argparse
import hashlib
import json
import platform
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch
from matplotlib.path import Path as MplPath
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, MaxNLocator

from theme import (ROOT, COLORS as C, LAYOUT, FONT, FONT_PATH, DIVERGING, DENSITY, effect_norm,
                   figure, axes_mm, text_mm, heading, data_axis, cell_label,
                   log_cells, save, lr_label, abbrev, ABBREVIATIONS)

parser = argparse.ArgumentParser()
parser.add_argument('--theme-only', action='store_true')
parser.add_argument('--main-only', action='store_true')
parser.add_argument('--figure1-only', action='store_true')
args = parser.parse_args()
PROV = []
CHECKS = []

def read(name, current=False):
    path = ROOT/'source_data'/('current' if current else 'archived')/name
    return pd.read_csv(path, sep='\t', float_precision='round_trip')

def old(name):
    return read(name+'.tsv.gz')

def check(name, condition, note=''):
    CHECKS.append(dict(check=name, status='PASS' if condition else 'FAIL', note=note))
    if not condition:
        pd.DataFrame(CHECKS).to_csv(ROOT/'reports/data_checks.tsv', sep='\t', index=False)
        raise AssertionError(name + ': ' + note)

def source(panel, data, files, columns, key, units, transform='none',
           rule='All supplied rows', n_definition='records', interval='none', version=''):
    path = ROOT/'source_data/panels'/f'{panel}.csv'
    data.to_csv(path, index=False, float_format='%.17g')
    PROV.append(dict(panel=panel, source_file=files, panel_csv=str(path.relative_to(ROOT)),
                     columns_used=columns, unique_key=key, units=units,
                     display_transform=transform, filtering_statistical_rule=rule,
                     n=len(data), n_definition=n_definition, interval_source=interval,
                     result_version=version))
    print(f'{panel}\t{files}\t{columns}\tn={len(data)}', flush=True)

SC = read('cosmx_scaling.tsv', True)
ATERA = pd.read_csv(ROOT/'source_data/current/atera_scale_performance.tsv', sep='\t')
ATERA = ATERA[ATERA.engine.isin(['SpatialESS', 'SpatialCellChat'])].copy()
ATERA['cells'] = pd.to_numeric(ATERA.N, errors='coerce')
ATERA['status'] = ATERA.status.astype(str)
FULL_REFERENCE_PATH = ROOT/'source_data/current/cervical_700k_reference_attempt.tsv'
FULL_REFERENCE = (pd.read_csv(FULL_REFERENCE_PATH, sep='\t')
                  if FULL_REFERENCE_PATH.exists() else pd.DataFrame())
if len(FULL_REFERENCE):
    full_result = FULL_REFERENCE.query("N == 715512 and engine == 'SpatialCellChat'")
    check('Full cervical reference attempt is unique and failed without output',
          len(full_result) == 1 and
          full_result.iloc[0].status == 'allocation_failure_under_96_GiB_virtual_memory_budget' and
          not bool(full_result.iloc[0].reference_result_generated))
    full_result = full_result.iloc[0]
    full_mask = ((ATERA.tissue == 'cervical') & (ATERA.fixture == 'full') &
                 (ATERA.engine == 'SpatialCellChat'))
    check('Full cervical reference attempt matches one saved 700k row', int(full_mask.sum()) == 1)
    ATERA.loc[full_mask, 'status'] = 'allocation_failure_under_memory_budget'
    ATERA.loc[full_mask, 'whole_process_seconds'] = full_result.whole_process_seconds
    ATERA.loc[full_mask, 'peak_RSS_GB_decimal'] = full_result.peak_RSS_GB_decimal
    ATERA.loc[full_mask, 'failure_reason'] = full_result.failure_reason
EXTENSION_PATH = ROOT/'source_data/current/cervical_extension_200k_400k.tsv'
EXTENSION = (pd.read_csv(EXTENSION_PATH, sep='\t') if EXTENSION_PATH.exists()
             else pd.DataFrame())
ATERA_AGREE = pd.read_csv(ROOT/'source_data/current/atera_scale_agreement.tsv', sep='\t')
performance_rows = []
for _, r in SC.iterrows():
    for engine, status_col, time_col, rss_col in [
        ('SpatialESS','ESS_status','ESS_wall_seconds','ESS_peak_rss_GiB'),
        ('SpatialCellChat','official_status','official_wall_seconds','official_peak_rss_gib')]:
        rss = r[rss_col]
        performance_rows.append(dict(dataset='CosMx', tissue='CancerousLiver',
            fixture=f'cosmx_{int(r.cells)}', N=int(r.cells), engine=engine,
            status=r[status_col], whole_process_seconds=r[time_col],
            peak_RSS_GB_decimal=rss*(2**30/1e9) if pd.notna(rss) else np.nan,
            inference_seconds=np.nan))
for _, r in ATERA.iterrows():
    performance_rows.append(dict(dataset='Atera', tissue=r.tissue,
        fixture=r.fixture, N=int(r.cells), engine=r.engine, status=r.status,
        whole_process_seconds=r.whole_process_seconds,
        peak_RSS_GB_decimal=r.peak_RSS_GB_decimal,
        failure_reason=getattr(r, 'failure_reason', ''),
        inference_seconds=r.inference_seconds))
D_PERF = pd.DataFrame(performance_rows)
PERFORMANCE_STYLE = json.loads((ROOT/'config/performance_styles.json').read_text())
FID = read('official_fidelity.tsv', True)
SP = read('spatial_sparsity.tsv', True)
PAIR = read('cosmx_50000_paired.tsv', True)
COS = FID[FID['case'].str.startswith('cosmx_')].copy()
COS['N'] = COS['case'].str.split('_').str[-1].astype(int)
COS = COS.sort_values('N')
LMM = old('lmm')
SEL = old('selected_interactions')
HM = old('effect_heatmap').set_index('feature_id').loc[SEL.feature_id].reset_index()
FOREST = old('forest')
EX = old('patient_examples').sort_values(['ligand','receptor','feature_id'])
SS, PS = old('slice_scores'), old('patient_scores')
LAB = old('volcano_labels').head(8).copy()
BUBBLES = old('pathway_pair_representatives')
EFFECT_LIMIT=max(np.nanmax(np.abs(HM[['standardized_effect','glm_standardized_effect']].values)),
                 np.nanmax(np.abs(BUBBLES.standardized_effect)))
BUBBLE_AREA_PER_UNIT = 4.4
META = old('sample_metadata')
VALID = LMM[LMM.status == 'ok'].copy()
IDS = SEL.set_index('feature_id')['display_id']
LAB['display_id'] = LAB.feature_id.map(IDS)
EX['display_id'] = EX.feature_id.map(IDS)

# Compare frozen plotting coefficients to the released model table by full key.
current_lmm = read('patient_random_intercept_lmm.tsv.gz', True).set_index('feature_id')
plotted_lmm = LMM.set_index('feature_id')
check('LMM full keys retained', set(current_lmm.index)==set(plotted_lmm.index))
for col in ['estimate','q_value','std_error','nonzero_units','total_units','patient_count']:
    check('LMM saved field '+col,
          np.allclose(current_lmm.loc[plotted_lmm.index,col],plotted_lmm[col],
                      rtol=1e-12,atol=0,equal_nan=True), 'Only text-serialization precision allowed')
check('Warning fits excluded from primary FDR', LMM.loc[LMM.status!='ok','q_value'].isna().all())
check('No zero or invalid q in displayed volcano', ((VALID.q_value>0)&(VALID.q_value<=1)).all())
check('28 selected heatmap rows preserved',len(HM)==28 and HM.feature_id.tolist()==SEL.feature_id.tolist())
check('12 forest features retained',len(FOREST)==12 and FOREST.feature_id.is_unique)
check('Stable display IDs', HM.display_id.tolist()==SEL.display_id.tolist())
check('Volcano eight archived labels', LAB.display_id.notna().all())
check('Feature CI from saved columns', np.isfinite(FOREST[['lower','upper']]).all().all())
check('Support denominator', (SEL.total_units==len(META)).all() and
      np.allclose(SEL.support_fraction,SEL.nonzero_units/SEL.total_units))
check('Patient identity and slice count', META.patient_id.nunique()==35 and len(META)==45)
check('Saved patient means', all(np.isclose(g.log1p_score.mean(),
      PS.loc[(PS.feature_id==key[0])&(PS.patient_id==key[1]),'log1p_score'].iloc[0])
      for key,g in SS.groupby(['feature_id','patient_id'])))
check('No double log transform', np.allclose(np.log1p(SS.score),SS.log1p_score))
check('CSR definitions retained', (SP.validation_status=='PASS').all() and
      np.allclose(SP.N**2,SP.N_squared) and np.allclose(SP.E/SP.N,SP.mean_out_degree) and
      np.allclose(SP.E/SP.N_squared,SP.graph_density) and
      (SP.csr_total_bytes==SP.offset_bytes+SP.neighbor_bytes+SP.distance_bytes).all())
check('No self-edges in radius graph counts', (SP.self_edges==0).all())
check('Fidelity exception preserved', len(FID)==66 and FID.strict_fidelity_pass.sum()==65)
mapping = SEL[['feature_id','display_id','interaction_name','sender_group_name',
               'receiver_group_name','ligand','receptor','pathway_name']].copy()
mapping['LR_display'] = mapping.apply(lr_label,axis=1)
mapping['sender_display'] = mapping.sender_group_name.map(abbrev)
mapping['receiver_display'] = mapping.receiver_group_name.map(abbrev)
mapping.to_csv(ROOT/'source_data/feature_mapping.csv',index=False)
pd.DataFrame(list(ABBREVIATIONS.items()),columns=['original','display']).to_csv(
    ROOT/'source_data/abbreviations.csv',index=False)

for panel,data,files,cols,key,units,transform,rule in [
 ('Figure1B',PAIR,'current/cosmx_50000_paired.tsv','probability_official;probability_ess',
  'sender_group_name;receiver_group_name;interaction_name','probability','log10 for strictly positive values',
  'Full-key pairing; count nonpositive, NA and unpaired records; Pearson on original probability'),
 ('Figure1C',COS,'current/official_fidelity.tsv','N;probability_pearson;max_probability_difference;exact_pvalue_fraction;active_jaccard;significant_disagreements_le005',
  'case','metric-specific','none','Six completed CosMx pairs; calls use raw p <= 0.05'),
 ('Figure1D',D_PERF,'current/cosmx_scaling.tsv;current/atera_scale_performance.tsv',
  'dataset;tissue;fixture;N;engine;status;whole_process_seconds',
  'dataset;fixture;engine','seconds','log10 x and y','Successful measurements plotted; all failed and unattempted states retained in source data'),
 ('Figure1E',D_PERF,'current/cosmx_scaling.tsv;current/atera_scale_performance.tsv',
  'dataset;tissue;fixture;N;engine;status;peak_RSS_GB_decimal',
  'dataset;fixture;engine','decimal GB','log10 x and y','GiB source measurements converted to decimal GB; failed resources are not plotted'),
 ('Figure1F',SP,'current/spatial_sparsity.tsv','N;N_squared;E;fixture;sampling_design',
  'fixture','ordered pairs / directed non-self edges','log10 x and y','All 14 PASS spatial fixtures'),
 ('Figure1G',SP,'current/spatial_sparsity.tsv','N;E;mean_out_degree;fixture',
  'fixture','neighbors per cell','horizontal bars + endpoint + numeric label; fixture labels include N',
  'All 14 fixtures; original fixture order within platform; no duplication by LR filter'),
 ('Figure1H',SP,'current/spatial_sparsity.tsv','N;graph_density;fixture',
  'fixture','E / N squared','log10 x and y','All fixtures; only nested CosMx values connected'),
 ('Figure1I',SC,'current/cosmx_scaling.tsv','cells;ESS_status;official_status;official_failure_reason',
  'cells','completion status','none','Paired filled; ESS-only open; official failed cross'),
 ('Figure2A',VALID,'archived/lmm.tsv.gz;archived/volcano_labels.tsv.gz','feature_id;estimate;q_value;status',
  'feature_id','log1p-score contrast PF minus Control; BH q','-log10(q); effect ticks show x10^-3',
  'status == ok; all valid points; original eight labels; no refitting or BH recomputation'),
 ('Figure2B',FOREST,'archived/forest.tsv.gz','feature_id;display_id;estimate;lower;upper;singular',
  'feature_id','log1p-score contrast','effect ticks show x10^-3','Exact archived 12 features'),
 ('Figure2C',HM,'archived/effect_heatmap.tsv.gz;archived/selected_interactions.tsv.gz',
  'feature_id;display_id;standardized_effect;glm_standardized_effect;q_value;glm_q;support_fraction;icc;singular',
  'feature_id','effect / response SD; nonzero-slice fraction; ICC','zero-centered effect fill; no new clipping',
  'Exact archived 28 rows; LMM/OLS share color range; dot q < .05; support denominator 45'),
 ('Figure2D_slices',SS,'archived/slice_scores.tsv.gz','feature_id;sample_id;patient_id;condition;log1p_score',
  'feature_id;sample_id','saved log1p score','no y transform; fixed x jitter seed 13','All 45 slices per original feature'),
 ('Figure2D_patients',PS,'archived/patient_scores.tsv.gz','feature_id;patient_id;condition;log1p_score',
  'feature_id;patient_id','mean saved log1p score','fixed x jitter seed 19','All 35 patients per original feature'),
 ('Figure2D_feature_identities',EX,'archived/patient_examples.tsv.gz',
  'feature_id;display_id;ligand;receptor;sender_group_name;receiver_group_name',
  'feature_id','feature identity','none','Three archived examples; exact saved identities'),
 ('Figure2E',BUBBLES,'archived/pathway_pair_representatives.tsv.gz',
  'feature_id;pathway_name;sender_group_name;receiver_group_name;standardized_effect;q_value',
  'feature_id','effect / SD; feature-level BH q','area = min(-log10(q),10)',
  'Archived 52 pathway/cell-pair representatives; no pathway-level enrichment inference')]:
    if args.figure1_only and not panel.startswith('Figure1'):
        continue
    source(panel,data,'source_data/'+files.replace(';',';source_data/'),cols,key,units,transform,rule,
           interval='saved lower/upper Wald columns' if panel in ['Figure2B','Figure2D_effects'] else 'none',
           version='ESS 0.1.4 / identified saved V3' if panel.startswith('Figure1') else '2026-09-07 frozen LMM')

pd.DataFrame(PROV).to_csv(ROOT/'reports/PANEL_INVENTORY.tsv',sep='\t',index=False)

def concise(x):
    if pd.isna(x): return 'NA'
    if x == 0: return '0'
    if abs(x) < .001: return f'{x:.2e}'
    if x == int(x): return str(int(x))
    return f'{x:.4g}'

def metric_table(ax, rows, headings, values, widths=None, fontsize=6, rowshade=True,
                 fidelity_tiles=False):
    ax.set_axis_off()
    nr, nc = len(rows), len(headings)
    widths = np.array(widths if widths is not None else [1.1]+[1]*nc,dtype=float)
    widths /= widths.sum()
    edges = np.r_[0,np.cumsum(widths)]
    centers = (edges[:-1]+edges[1:])/2
    ax.set_xlim(0,1); ax.set_ylim(nr+.4,-1.7)
    for j, label in enumerate(['']+headings):
        ax.text(centers[j],-1.03,label,ha='center',va='center',size=fontsize,linespacing=1.05)
    ax.hlines(-.47,0,1,color=C['secondary'],lw=.45)
    for i,row in enumerate(rows):
        if rowshade and not fidelity_tiles and i%2==0:
            ax.add_patch(Rectangle((0,i-.43),1,.86,fc='#F1F4F6',ec='none'))
        ax.text(edges[0]+.005,i,row,ha='left',va='center',size=fontsize)
        for j,v in enumerate(values[i]):
            if fidelity_tiles:
                label=headings[j]
                perfect_one=label in {'r','Exact p','Jaccard'} and str(v)=='1'
                perfect_zero=label=='Δ calls' and str(v)=='0'
                color='#E2EEF7' if perfect_one else '#EDF1F4' if perfect_zero else '#F5F7F9'
                ax.add_patch(Rectangle((edges[j+1]+.004,i-.44),widths[j+1]-.008,.88,
                                       fc=color,ec='none'))
            ax.text(centers[j+1],i,str(v),ha='center',va='center',size=fontsize,
                    weight='bold' if fidelity_tiles and j==nc-1 else 'normal')

def agreement(ax, d, xcol='probability_official', ycol='probability_ess',
              summary=None, colorbar=True, xlabel='V3 probability (log10)',
              ylabel='ESS probability (log10)'):
    good = np.isfinite(d[xcol]) & np.isfinite(d[ycol]) & (d[xcol]>0) & (d[ycol]>0)
    check('Positive finite agreement '+ax.figure._name+str(len(CHECKS)),good.all(),
          f'n={len(d)}; invalid={(~good).sum()}')
    x,y = np.log10(d.loc[good,xcol]),np.log10(d.loc[good,ycol])
    lo, hi = min(x.min(),y.min()),max(x.max(),y.max())
    pad=(hi-lo)*.05
    hb=ax.hexbin(x,y,gridsize=34,mincnt=1,cmap=DENSITY,norm=LogNorm(),linewidths=0)
    # Explicit polygons preserve the same hexbin counts across PDF/SVG backends.
    hb.update_scalarmappable()
    hex_vertices=hb.get_paths()[0].vertices.copy()
    counts=np.asarray(hb.get_array())
    centers=np.asarray(hb.get_offsets())
    check('Hexbin count conservation '+ax.figure._name,counts.sum()==len(x))
    for center,color in zip(centers,hb.get_facecolors()):
        ax.add_patch(Polygon(hex_vertices+center,facecolor=color,edgecolor='none'))
    mappable=mpl.cm.ScalarMappable(norm=hb.norm,cmap=hb.cmap)
    hb.remove();hb=mappable
    pd.DataFrame({'log10_x_center':centers[:,0], 'log10_y_center':centers[:,1],
                  'records':counts.astype(int)}).to_csv(
        ROOT/'source_data/panels'/f'{ax.figure._name}_probability_bins.csv',index=False)
    ax.plot([lo-pad,hi+pad],[lo-pad,hi+pad],color=C['secondary'],lw=.45,ls='--')
    ax.set_xlim(lo-pad,hi+pad);ax.set_ylim(lo-pad,hi+pad);ax.set_aspect('equal',adjustable='box')
    ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
    if summary is not None:
        label=f'n = {len(d):,}\nr = {summary[0]:.6g}\nmax |ΔP| = {summary[1]:.2e}'
        ax.text(.04,.96,label,transform=ax.transAxes,va='top',size=5.6,
                bbox=dict(fc='white',ec='none',pad=1,alpha=.9))
    if colorbar:
        pos=ax.get_position(); fw,fh=ax.figure.get_size_inches()*25.4
        cax=ax.figure.add_axes([pos.x1+1.0/fw,pos.y0,1.5/fw,pos.height*.53])
        cb=ax.figure.colorbar(hb,cax=cax)
        cb.solids.set_rasterized(False)
        cb.solids.set_edgecolor('face')
        cb.outline.set_visible(False);cb.ax.tick_params(labelsize=5.3,width=.4,length=1.5)
        cb.ax.set_title('Records',fontsize=5.5,pad=3)
    return hb

def probability_marginals(ax,d,xcol,ycol,marginal_height=2.8,marginal_width=6):
    fig=ax.figure;fig.canvas.draw()
    pos=ax.get_position();fw,fh=fig.get_size_inches()*25.4
    bins=np.linspace(*ax.get_xlim(),25)
    top=fig.add_axes([pos.x0,pos.y1+.6/fh,pos.width,marginal_height/fh])
    right=fig.add_axes([pos.x1+.8/fw,pos.y0,marginal_width/fw,pos.height])
    cx,_=np.histogram(np.log10(d[xcol]),bins=bins)
    cy,_=np.histogram(np.log10(d[ycol]),bins=bins)
    top.bar(bins[:-1],cx,width=np.diff(bins),align='edge',color=C['reference'],lw=0)
    right.barh(bins[:-1],cy,height=np.diff(bins),align='edge',color=C['ESS'],lw=0)
    top.set_xlim(ax.get_xlim());right.set_ylim(ax.get_ylim())
    top.set_axis_off();right.set_axis_off()
    check('Marginal conservation '+fig._name,cx.sum()==len(d) and cy.sum()==len(d))
    return pd.DataFrame({'log10_bin_left':bins[:-1],'log10_bin_right':bins[1:],
                         'reference_records':cx,'ESS_records':cy})

def scaling(ax, what):
    value = 'whole_process_seconds' if what == 'time' else 'peak_RSS_GB_decimal'
    cohorts = [
        ('CosMx', 'CancerousLiver', '-'),
        ('Atera', 'breast', '--'),
        ('Atera', 'cervical', ':'),
    ]
    plotted = 0
    for dataset, tissue, linestyle in cohorts:
        cohort = D_PERF[(D_PERF.dataset == dataset) & (D_PERF.tissue == tissue)]
        completed = cohort[cohort.status == 'completed'].dropna(subset=[value])
        paired = completed.pivot_table(index='N', columns='engine', values=value, aggfunc='first')
        paired = paired.dropna(subset=['SpatialESS', 'SpatialCellChat']).sort_index()
        if len(paired) > 1:
            ax.fill_between(paired.index.to_numpy(dtype=float),
                            paired['SpatialESS'].to_numpy(dtype=float),
                            paired['SpatialCellChat'].to_numpy(dtype=float),
                            color='#DCE4EB', alpha=.48, linewidth=0, zorder=1)
        for method, color, marker in [
            ('SpatialESS', C['ESS'], 'o'),
            ('SpatialCellChat', '#D55E00', 's')]:
            d = completed[completed.engine == method].sort_values('N')
            if d.empty:
                continue
            ax.plot(d.N, d[value], color=color, marker=marker, ms=4.2,
                    mfc=color, mec='white', mew=.55, lw=1.2,
                    ls=linestyle, zorder=3 if method == 'SpatialESS' else 2)
            plotted += len(d)
    check('D/E successful point count '+what, plotted ==
          len(D_PERF[D_PERF.status == 'completed'].dropna(subset=[value])))
    log_cells(ax, limits=(200, 900000));ax.set_yscale('log');ax.minorticks_off()
    ax.set_ylabel('Wall-clock time (s)' if what=='time' else 'Peak RSS (GB, decimal)')
    if what=='time':
        ax.set_yticks([10,100,1000,10000],['10','100','1,000','10,000']);ax.set_ylim(8,8e4)
    else:
        ax.set_yticks([.1,1,10,100],['0.1','1','10','100']);ax.set_ylim(.08,150)
    if what == 'time':
        b = ATERA[(ATERA.tissue == 'breast') & (ATERA.cells == 50000)]
        c = ATERA[(ATERA.tissue == 'cervical') & (ATERA.cells == 50000)]
        ratios = []
        for d in [b, c]:
            e = d[d.engine == 'SpatialESS'].whole_process_seconds
            r = d[d.engine == 'SpatialCellChat'].whole_process_seconds
            if len(e) and len(r): ratios.append(float(r.iloc[0] / e.iloc[0]))
        note = f'Atera 50k: {min(ratios):.0f}–{max(ratios):.0f}× faster' if ratios else ''
    else:
        b = ATERA[(ATERA.tissue == 'breast') & (ATERA.cells == 50000)]
        c = ATERA[(ATERA.tissue == 'cervical') & (ATERA.cells == 50000)]
        ratios = []
        for d in [b, c]:
            e = d[d.engine == 'SpatialESS'].peak_RSS_GB_decimal
            r = d[d.engine == 'SpatialCellChat'].peak_RSS_GB_decimal
            if len(e) and len(r): ratios.append(float(r.iloc[0] / e.iloc[0]))
        note = f'Atera 50k: {min(ratios):.1f}–{max(ratios):.1f}× lower RSS' if ratios else ''
    ax.text(.035,.97,note,va='top',size=5.0,weight='bold',color=C['ESS'],transform=ax.transAxes)

SERIES=['CosMx central','Xenium central','Xenium stratified','Visium HD ovarian']
SCOL=[C['ESS'],C['positive'],C['positive'],C['green']]
SMARK=['o','^','v','D']
SP['series']=np.select([SP.fixture.str.startswith('cosmx'),SP.fixture.str.startswith('cervical_center'),
                       SP.fixture.str.startswith('cervical_stratified')],SERIES[:3],default=SERIES[3])

def spatial(ax, field):
    for series,color,marker in zip(SERIES,SCOL,SMARK):
        d=SP[SP.series==series].sort_values('N')
        ax.plot(d.N,d[field],ls='-' if series==SERIES[0] else 'none',color=color,
                marker=marker,ms=4.6,mew=.85,mfc='white' if series==SERIES[2] else color,lw=1.05)
    log_cells(ax)
    if field!='mean_out_degree':
        ax.set_yscale('log');ax.minorticks_off()
    if field=='E':
        d=SP[SP.series==SERIES[0]].sort_values('N')
        ax.fill_between(d.N,d.E,d.N_squared,color='#E4EEF7',lw=0,zorder=0)
        r=SP.sort_values('N').drop_duplicates('N')
        ax.plot(r.N,r.N_squared,color=C['reference'],ls='--',lw=.7)
        ax.text(.96,.94,'N²',ha='right',va='top',size=6,transform=ax.transAxes)
        ax.set_ylabel('Pairs / directed edges');ax.set_yticks([1e2,1e5,1e8,1e11]);ax.set_ylim(20,8e11)
        largest=SP[SP.series==SERIES[0]].sort_values('N').iloc[-1]
        ax.text(.03,.86,f'N²/E = {largest.N_squared/largest.E:,.0f}\nat {cell_label(largest.N)} cells',
                transform=ax.transAxes,size=5.2,color=C['secondary'],va='top')
    elif field=='mean_out_degree':
        ax.set_ylabel('Mean neighbors, E/N');ax.set_ylim(-1.2,SP.mean_out_degree.max()*1.13)
        ax.set_yticks([0,10,20,30])
    else:
        ax.set_ylabel('Graph density, E/N²');ax.set_yticks([1e-4,1e-3,1e-2,1e-1]);ax.set_ylim(2e-5,.15)

def spatial_legend(fig,y):
    labels=['CosMx','Xenium central','Xenium stratified','Visium HD ovarian']
    handles=[Line2D([],[],color=c,marker=m,ls='none',ms=3.5,mew=.7,
                    mfc='white' if i==2 else c,label=labels[i])
             for i,(c,m) in enumerate(zip(SCOL,SMARK))]
    handles.append(Line2D([],[],color=C['reference'],ls='--',lw=.8,
                          label='N²: all possible ordered pairs'))
    h=fig.get_size_inches()[1]*25.4
    fig.legend(handles=handles,ncol=5,loc='center',bbox_to_anchor=(.52,1-y/h),fontsize=5.2,
               handletextpad=.35,columnspacing=1.05)

def volcano(ax, compact=False):
    x=VALID.estimate;y=-np.log10(VALID.q_value)
    sig=VALID.q_value<.05
    for mask,color,size,alpha in [(~sig,C['nonsignificant'],2.5,.66),
                                  (sig&(x<0),C['negative'],4.1,.79),
                                  (sig&(x>=0),C['positive'],4.1,.79)]:
        ax.scatter(x[mask],y[mask],s=size,alpha=alpha,c=color,edgecolors='none',rasterized=True)
    ax.axhline(-np.log10(.05),ls='--',lw=.5,c=C['secondary']);ax.axvline(0,lw=.4,c=C['secondary'])
    ax.set_ylim(-.6,max(y.max()*1.08,25));ax.set_xlim(x.min()*1.12,x.max()*1.17)
    ax.set_ylabel('−log10(BH q)');ax.set_xlabel('PF − Control effect (× 0.001)')
    ax.set_xticks([-.001,0,.001,.002,.003]);ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x*1000:g}'))
    if compact:
        return
    # Labels occupy a separate column; paths are not covered by another label.
    # Every arrow tip uses the saved, full-key-matched effect and BH q coordinate.
    ordered=LAB.sort_values('q_value',kind='stable')
    ax.figure.canvas.draw()
    fw,fh=ax.figure.get_size_inches()*25.4
    pos=ax.get_position()
    point_mm=ax.transData.transform(np.column_stack([ordered.estimate,-np.log10(ordered.q_value)]))
    point_mm=point_mm/ax.figure.dpi*25.4
    rows=[]
    for i,(_,r) in enumerate(ordered.iterrows()):
        label_x=pos.x1*fw+6.0
        label_top=15.5+i*3.0
        anchor=np.array([label_x-.65,fh-label_top-1.05])
        tip=point_mm[i]
        elbow=np.array([pos.x1*fw+.8,tip[1]])
        verts=np.array([anchor,elbow,tip])/[fw,fh]
        path=MplPath(verts,[MplPath.MOVETO,MplPath.LINETO,MplPath.LINETO])
        arrow=FancyArrowPatch(path=path,transform=ax.figure.transFigure,
                              arrowstyle='-|>',mutation_scale=4.6,linewidth=.55,
                              color='#596D80',clip_on=False,zorder=5)
        arrow.set_gid('leader_'+r.display_id)
        ax.figure.add_artist(arrow)
        text_mm(ax.figure,label_x,label_top,r.display_id+'  '+lr_label(r),
                size=5.05,weight='bold')
        ax.scatter(r.estimate,-np.log10(r.q_value),s=13,facecolors='white',
                   edgecolors=C['positive'],linewidths=.75,zorder=6)
        rows.append(dict(feature_id=r.feature_id,display_id=r.display_id,
                         LR_display=lr_label(r),estimate=r.estimate,q_value=r.q_value,
                         target_x=r.estimate,target_y=-np.log10(r.q_value),
                         label_x_mm=label_x,label_top_mm=label_top,
                         arrow_start_x_mm=anchor[0],arrow_start_y_mm=anchor[1],
                         arrow_tip_x_mm=tip[0],arrow_tip_y_mm=tip[1]))
    audit=pd.DataFrame(rows)
    frozen=VALID.set_index('feature_id').loc[audit.feature_id]
    check('Volcano leader full-key identity',np.array_equal(frozen.estimate.values,audit.target_x.values)
          and np.array_equal(-np.log10(frozen.q_value.values),audit.target_y.values))
    source('Figure2A_labels',audit,'source_data/archived/volcano_labels.tsv.gz;source_data/archived/lmm.tsv.gz',
           'feature_id;display_id;ligand;receptor;estimate;q_value','feature_id',
           'effect; -log10 BH q; physical mm','exact target coordinate; fixed display label positions',
           'Same eight frozen labels; one explicit vector arrow per complete feature key')

def effect_forest(ax):
    """A lollipop-style forest plot with a visible zero-to-estimate bar."""
    max_upper=FOREST.upper.max()
    value_x=max_upper*1.30
    ax.set_xlim(-max_upper*.08,value_x*1.08)
    for i in range(len(FOREST)):
        if i % 2 == 0:
            ax.axhspan(i-.44,i+.44,color='#F5F7F9',lw=0,zorder=0)
    for i,(_,r) in enumerate(FOREST.iterrows()):
        color=C['positive'] if r.estimate>0 else C['negative']
        ax.plot([0,r.estimate],[i,i],c='#E8A7B1',lw=4.0,solid_capstyle='butt',zorder=1)
        ax.plot([r.lower,r.upper],[i,i],c=color,lw=1.25,zorder=3)
        ax.plot([r.lower,r.lower],[i-.14,i+.14],c=color,lw=.75,zorder=3)
        ax.plot([r.upper,r.upper],[i-.14,i+.14],c=color,lw=.75,zorder=3)
        ax.plot(r.estimate,i,'o',ms=5.3,mfc='white' if r.singular else color,
                mec=color,mew=1.0,zorder=4)
        ax.text(value_x,i,f'{r.estimate*1000:+.2f}',ha='right',va='center',
                size=5.1,color=C['secondary'],weight='bold')
    ax.axvline(0,lw=.5,c=C['secondary']);ax.set_ylim(len(FOREST)-.5,-.6)
    ax.set_yticks(range(len(FOREST)),[r.display_id+'  '+lr_label(r) for _,r in FOREST.iterrows()],fontsize=5.65)
    ax.tick_params(axis='y',length=0,pad=3);ax.spines['left'].set_visible(False)
    ax.set_xlabel('PF − Control effect (× 0.001; 95% CI)')
    ax.set_xticks([0,.001,.002,.003]);ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x*1000:g}'))

def effect_heatmap(fig):
    x,y,w,h=heading(fig,'C','Communication effects and model comparison')
    ax=axes_mm(fig,x,y+5.5,w,h-10)
    ax.set_axis_off();ax.set_xlim(0,176)
    # Fixed row spacing, additional whitespace only between pathways.
    yvals=[]; pos=0; previous=None
    for _,r in HM.iterrows():
        if previous is not None and previous!=r.pathway_name:pos+=.23
        yvals.append(pos);pos+=1;previous=r.pathway_name
    ax.set_ylim(pos+.1,-2.1)
    xs=dict(path=0,ID=17,LR=25,group=69,LMM=114,OLS=127,support=143,icc=158,singular=172)
    for key,label in [('path','Pathway'),('ID','ID'),('LR','Ligand–receptor'),
                      ('group','Sender = receiver group'),('LMM','LMM'),('OLS','OLS'),
                      ('support','Nonzero-slice\nfraction'),('icc','Patient\nICC'),('singular','Sing.')]:
        ax.text(xs[key],-1.35,label,size=5.7,ha='left' if key in ['path','ID','LR','group'] else 'center',va='center')
    lim=EFFECT_LIMIT
    norm=effect_norm(lim)
    check('Heatmap zero center',norm(0)==.5)
    for block,p in enumerate(HM.pathway_name.unique()):
        loc=[yy for yy,pp in zip(yvals,HM.pathway_name) if p==pp]
        if block % 2 == 0:
            ax.add_patch(Rectangle((0,min(loc)-.48),176,max(loc)-min(loc)+.96,
                                   fc='#F5F7F9',ec='none',zorder=-2))
        ax.plot([-.4,-.4],[min(loc)-.40,max(loc)+.40],color=C['secondary'],
                lw=.6,clip_on=False)
        ax.text(0,np.mean(loc),p,size=5.55,va='center',color=C['secondary'],weight='bold')
    ax.hlines(-.68,0,176,color='#AAB3BC',lw=.4)
    for yy,(_,r) in zip(yvals,HM.iterrows()):
        for xx,txt,col in [(xs['ID'],r.display_id,C['secondary']),
                            (xs['LR'],lr_label(r),C['ink']),
                            (xs['group'],abbrev(r.sender_group_name),C['secondary'])]:
            ax.text(xx,yy,txt,va='center',size=5.7,color=col)
        for xx,ec,q in [(xs['LMM'],'standardized_effect','q_value'),(xs['OLS'],'glm_standardized_effect','glm_q')]:
            color=DIVERGING(norm(r[ec])) if pd.notna(r[ec]) else C['pale']
            ax.add_patch(Rectangle((xx-5.3,yy-.45),10.6,.9,fc=color,ec='none'))
            label=f'{r[ec]:.2f}' if pd.notna(r[ec]) else 'NA'
            rgb=np.asarray(color[:3])
            lettering='white' if np.dot(rgb,[.2126,.7152,.0722])<.52 else C['ink']
            ax.text(xx-.55,yy,label,ha='center',va='center',size=5.05,color=lettering)
            if pd.notna(r[q]) and r[q]<.05:
                ax.plot(xx+4.1,yy,'o',ms=1.2,c=lettering)
        for xx,field,color in [(xs['support'],'support_fraction',C['ESS']),(xs['icc'],'icc',C['green'])]:
            ax.plot([xx-5,xx+5],[yy,yy],lw=2.8,c=C['pale'],solid_capstyle='butt')
            if pd.notna(r[field]):
                ax.plot([xx-5,xx-5+10*r[field]],[yy,yy],lw=2.8,c=color,solid_capstyle='butt')
            else:ax.text(xx,yy,'NA',size=5,va='center',ha='center')
        if r.singular:ax.plot(xs['singular'],yy,'s',ms=2.2,c=C['reference'])
    cax=axes_mm(fig,75,y+h-4.8,24,1.5)
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=DIVERGING),cax=cax,orientation='horizontal',ticks=[-1,0,1])
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face')
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=5.3,length=1.2,pad=.8)
    text_mm(fig,55,y+h-4.9,'Effect / SD',size=5.5)
    text_mm(fig,104,y+h-4.9,'• q < 0.05; bars: 0–1; nonzero / 45 slices',size=5.3)

def patient_examples(fig):
    x,y,w,h=heading(fig,'D','Slices and patient summaries')
    text_mm(fig,91,y+.3,'Small: slices   Open: patient means; 35 patients / 45 slices',size=5.4)
    rngs=np.random.default_rng(13);rngp=np.random.default_rng(19)
    for i,(_,r) in enumerate(EX.iterrows()):
        xx=x+i*59
        text_mm(fig,xx+30,y+4,lr_label(r)+'  |  '+r.display_id,ha='center',size=6.1,weight='bold')
        text_mm(fig,xx+30,y+7.1,abbrev(r.sender_group_name)+' → '+abbrev(r.receiver_group_name),
                ha='center',size=5.5)
        ax=axes_mm(fig,xx+10,y+11.0,45,16.0)
        s=SS[SS.feature_id==r.feature_id];p=PS[PS.feature_id==r.feature_id]
        for j,(condition,color) in enumerate([('Control',C['negative']),('Disease',C['positive'])]):
            a=s[s.condition==condition];b=p[p.condition==condition]
            ax.scatter(j+rngs.uniform(-.16,.16,len(a)),a.log1p_score,s=7,c=color,alpha=.46,edgecolors='none')
            ax.scatter(j+rngp.uniform(-.09,.09,len(b)),b.log1p_score,s=15,facecolors='white',edgecolors=color,lw=.75)
        ax.set_xlim(-.45,1.45);ax.set_ylim(-.00015,max(s.log1p_score)*1.1)
        ax.set_xticks([0,1],['Control','PF']);ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:.3f}'))
        if i==0:ax.set_ylabel('log1p(score)',size=6)
        ax.tick_params(labelsize=5.5)

def bubble(fig):
    x,y,w,h=heading(fig,'E','Pathway / sender–receiver representatives')
    text_mm(fig,98,y+.25,'Effect / SD; area: −log10(feature q), capped at 10',size=5.3)
    pairs=sorted(BUBBLES['pair'].unique());paths=sorted(BUBBLES.pathway_name.unique())
    ax=axes_mm(fig,x+19,y+4.5,128,19.5)
    norm=effect_norm(EFFECT_LIMIT)
    px={v:i for i,v in enumerate(pairs)};py={v:i for i,v in enumerate(paths)}
    ax.scatter(BUBBLES['pair'].map(px),BUBBLES.pathway_name.map(py),
               s=np.minimum(-np.log10(BUBBLES.q_value),10)*BUBBLE_AREA_PER_UNIT,
               c=BUBBLES.standardized_effect,cmap=DIVERGING,norm=norm,edgecolors='white',lw=.25)
    ax.set_xlim(-.5,len(pairs)-.5);ax.set_ylim(-.5,len(paths)-.5)
    short_groups={'Activated Fibrotic FBs':'Act. FBs', 'CD8+ T-cells':'CD8+ T',
                  'Interstitial Macrophages':'Interst. Mph', 'Alveolar FBs':'Alv. FBs',
                  'CD4+ T-cells':'CD4+ T'}
    short_labels=[]
    for pair in pairs:
        row=BUBBLES[BUBBLES['pair']==pair].iloc[0]
        short_labels.append(short_groups.get(row.sender_group_name,row.sender_group_name)+'\n→ '+
                            short_groups.get(row.receiver_group_name,row.receiver_group_name))
    ax.set_yticks(range(len(paths)),paths,size=5.4)
    ax.set_xticks(range(len(pairs)),short_labels,size=5.4)
    ax.tick_params(length=0,pad=1)
    for spine in ax.spines.values():spine.set_visible(False)
    cax=axes_mm(fig,x+155,y+5,18,1.5)
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=DIVERGING),cax=cax,orientation='horizontal',ticks=[-1,0,1]);cb.outline.set_visible(False)
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face')
    cb.ax.tick_params(labelsize=5.2,length=1.0,pad=.5)
    for j,v in enumerate([2,5,10]):
        ax2=axes_mm(fig,x+154+j*7,y+11,7,6);ax2.set_axis_off();ax2.set_xlim(0,1);ax2.set_ylim(0,1)
        ax2.scatter(.5,.7,s=v*BUBBLE_AREA_PER_UNIT,c=C['reference']);ax2.text(.5,.1,str(v),ha='center',size=5.1)
    pairmap=pd.DataFrame([dict(pair_id=f'P{i+1}',pair=p,
      sender=BUBBLES.loc[BUBBLES['pair']==p,'sender_group_name'].iloc[0],
      receiver=BUBBLES.loc[BUBBLES['pair']==p,'receiver_group_name'].iloc[0]) for i,p in enumerate(pairs)])
    pairmap.to_csv(ROOT/'source_data/pair_mapping.csv',index=False)
    # Short pair IDs remain explicitly mapped in the full-width extended panel below.
    return pairmap

def fidelity_ribbon(fig):
    x,y,w,h=heading(fig,'C','Fidelity across paired scales')
    ax=axes_mm(fig,x,y+4.4,w,h-9.0)
    ax.set_axis_off();ax.set_xlim(0,176);ax.set_ylim(5.05,-.5)
    edges=np.linspace(53,176,7)
    fields=[('probability_pearson','Probability r'),
            ('max_probability_difference','Max |Δ probability|'),
            ('exact_pvalue_fraction','Exact p fraction'),
            ('active_jaccard','Active-record Jaccard'),
            ('significant_disagreements_le005','Call disagreements')]
    for j,(_,r) in enumerate(COS.iterrows()):
        center=(edges[j]+edges[j+1])/2
        ax.text(center,-.25,cell_label(r.N),ha='center',va='center',size=5.7,weight='bold')
        for i,(field,label) in enumerate(fields):
            color=['#DFECF6','#F3F5F7','#DFECF6','#DFECF6','#E8EFF3'][i]
            ax.add_patch(Rectangle((edges[j]+.55,i+.17),edges[j+1]-edges[j]-1.1,.72,
                                   fc=color,ec='none'))
            ax.text(center,i+.53,concise(r[field]),ha='center',va='center',
                    size=5.55,weight='bold' if i==4 else 'normal')
    strict = ATERA_AGREE[ATERA_AGREE.strict_pass.astype(bool)]
    max_error = ATERA_AGREE.max_absolute_probability_error.max()
    check('Atera paired tolerance summary',len(ATERA_AGREE)==10 and len(strict)==10 and
          ATERA_AGREE.active_jaccard.eq(1).all() and
          ATERA_AGREE.exact_pvalue_fraction.eq(1).all() and
          ATERA_AGREE.significance_disagreements.eq(0).all(),
          f'10/10 strict; max probability error {max_error:.17g}')
    text_mm(fig,x+1,y+h-1.25,
            f'Atera 1–50k: {len(strict)}/10 strict passes; max |ΔP| = {max_error:.2e}; p identical; 0 call disagreements',
            size=5.0,color=C['secondary'])
    for i,(_,label) in enumerate(fields):
        ax.text(0,i+.53,label,va='center',size=5.5,color=C['secondary'])

def quantified_bars(ax,values,labels,colors,markers=None,open_points=None,
                    formats=None,fontsize=5.3,height=.46):
    """Shared bar/endpoint/value style; no uncertainty is generated."""
    values=np.asarray(values,dtype=float)
    maximum=max(values.max(),1)
    if markers is None:markers=['o']*len(values)
    if open_points is None:open_points=[False]*len(values)
    for i,(v,c,m,hollow) in enumerate(zip(values,colors,markers,open_points)):
        if i%2==0:ax.axhspan(i-.46,i+.46,color='#F4F7F9',lw=0,zorder=0)
        ax.barh(i,v,color=c,alpha=.42,height=height,linewidth=0,zorder=1)
        ax.plot(v,i,marker=m,ls='none',ms=3.5,mfc='white' if hollow else c,
                mec=c,mew=.7,zorder=3)
        label=formats[i] if formats is not None else f'{v:,.0f}'
        ax.text(maximum*1.31,i,label,ha='right',va='center',size=fontsize,
                color=C['secondary'],weight='bold')
    ax.set_yticks(range(len(labels)),labels,fontsize=fontsize)
    ax.set_ylim(len(labels)-.5,-.5);ax.set_xlim(-maximum*.03,maximum*1.38)
    ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0,pad=2)

def local_degree_lollipops(fig):
    x,y,w,h=heading(fig,'G','Local neighbors per cell')
    ax=axes_mm(fig,x+23,y+6,w-25,h-17)
    labels=[];rows=[]
    for s in SERIES:
        d=SP[SP.series==s].sort_values('N')
        prefix={'CosMx central':'CosMx', 'Xenium central':'X cen.',
                'Xenium stratified':'X strat.','Visium HD ovarian':'VHD'}[s]
        for _,r in d.iterrows():
            labels.append(prefix+' '+cell_label(r.N));rows.append(r)
    values=[r.mean_out_degree for r in rows]
    quantified_bars(ax,values,labels,[SCOL[SERIES.index(r.series)] for r in rows],
                    markers=[SMARK[SERIES.index(r.series)] for r in rows],
                    open_points=[r.series==SERIES[2] for r in rows],
                    formats=[f'{v:.3f}' if v<1 else f'{v:.2f}' for v in values],fontsize=5.0)
    ax.set_xticks([0,10,20,30]);ax.set_xlabel('Neighbors, E/N',fontsize=6.1)
    ax.spines['left'].set_visible(False)
    for bound in [8.5,10.5,12.5]:
        ax.axhline(bound,color='#E2E8ED',lw=.4,zorder=0)

def fig1():
    f=figure('Figure1')
    x,y,w,h=heading(f,'B','Probability fidelity · CosMx 50k')
    ax=axes_mm(f,x+11,y+8,35,26)
    z=COS[COS.N==50000].iloc[0]
    mappable=agreement(ax,PAIR,summary=(z.probability_pearson,z.max_probability_difference),
                       colorbar=False,xlabel='SpatialCellChat probability (log10)',
                       ylabel='SpatialESS probability (log10)')
    bins=np.linspace(*ax.get_xlim(),25)
    top=axes_mm(f,x+11,y+4.8,35,2.5)
    right=axes_mm(f,x+47,y+8,5.5,26)
    count_x,_=np.histogram(np.log10(PAIR.probability_official),bins=bins)
    count_y,_=np.histogram(np.log10(PAIR.probability_ess),bins=bins)
    top.bar(bins[:-1],count_x,width=np.diff(bins),align='edge',color=C['reference'],lw=0)
    right.barh(bins[:-1],count_y,height=np.diff(bins),align='edge',color=C['ESS'],lw=0)
    top.set_xlim(ax.get_xlim());right.set_ylim(ax.get_ylim())
    for marginal in [top,right]:
        marginal.set_axis_off()
    check('F1 probability marginal count conservation',count_x.sum()==len(PAIR) and count_y.sum()==len(PAIR))
    source('Figure1B_marginals',pd.DataFrame({'log10_bin_left':bins[:-1],
               'log10_bin_right':bins[1:],'reference_records':count_x,'ESS_records':count_y}),
           'source_data/current/cosmx_50000_paired.tsv','probability_official;probability_ess',
           'log10_bin_left;log10_bin_right','records','24 equal log10 probability bins',
           'All 18,484 strictly positive paired probabilities',n_definition='display histogram bins')
    cax=axes_mm(f,x+57,y+18,1.3,15)
    cb=f.colorbar(mappable,cax=cax);cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face')
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=5.2,width=.4,length=1.3)
    cb.ax.set_title('Records',fontsize=5.3,pad=3)
    fidelity_ribbon(f)
    scaling(data_axis(f,'D','Whole-process runtime',left=11,top=9,bottom=8),'time')
    scaling(data_axis(f,'E','Peak resident memory',left=11,top=9,bottom=8),'memory')
    method_handles=[Line2D([],[],color=C['ESS'],marker='o',ms=4.0,lw=1.1,
                           label='SpatialESS'),
                    Line2D([],[],color='#D55E00',marker='s',ms=4.0,lw=1.1,
                           label='SpatialCellChat')]
    input_handles=[Line2D([],[],color=C['secondary'],ls='-',lw=1.1,label='CosMx CancerousLiver'),
                   Line2D([],[],color=C['secondary'],ls='--',lw=1.1,label='Atera breast'),
                   Line2D([],[],color=C['secondary'],ls=':',lw=1.2,label='Atera cervical'),
                   Rectangle((0,0),1,1,fc='#DCE4EB',ec='none',label='Shading: paired completed sizes')]
    f.legend(handles=method_handles,loc='center',bbox_to_anchor=(.5,1-165.4/250),
             ncol=2,fontsize=5.0,handletextpad=.4,columnspacing=1.2)
    f.legend(handles=input_handles,loc='center',bbox_to_anchor=(.5,1-168.5/250),
             ncol=4,fontsize=5.0,handletextpad=.4,columnspacing=1.0)
    spatial(data_axis(f,'F','Potential pairs and observed edges',left=11,top=8,bottom=11),'E')
    local_degree_lollipops(f)
    spatial(data_axis(f,'H','Graph density',left=10,top=8,bottom=11),'graph_density')
    spatial_legend(f,210.5)
    text_mm(f,4,215,'Directed non-self edges; effective radius = 35 μm. Runtime and memory points are measured values.',size=5.0)
    x,y,w,h=heading(f,'I','Completion and comparison scope')
    text_mm(f,62,y+4.4,'Filled: paired complete   Open: SpatialESS only   ×: reference failed',
            size=5.0,color=C['secondary'])
    rows=[]
    for dataset, data, id_col, ess_col, ref_col in [
        ('CosMx', SC, 'cells', 'ESS_status', 'official_status')]:
        for method in ['SpatialCellChat','SpatialESS']:
            for _,r in data.iterrows():
                status = r[ref_col] if method == 'SpatialCellChat' else r[ess_col]
                rows.append(dict(dataset=dataset,tissue='CancerousLiver',N=int(r[id_col]),
                                 method=method,status=status))
    for tissue in ['breast','cervical']:
        d=ATERA[ATERA.tissue==tissue]
        for method in ['SpatialCellChat','SpatialESS']:
            for _,r in d[d.engine==method].iterrows():
                rows.append(dict(dataset='Atera',tissue=tissue,N=int(r.cells),
                                 method=method,status=r.status))
    status_df=pd.DataFrame(rows)
    status_df['record_type']='saved_run_record'
    status_df['source_note']='Saved benchmark run status'
    status_df['display_label']=''
    if len(FULL_REFERENCE):
        full_reference_row = FULL_REFERENCE.query(
            "N == 715512 and engine == 'SpatialCellChat'").iloc[0]
        full_mask = ((status_df.dataset == 'Atera') & (status_df.tissue == 'cervical') &
                     (status_df.N == 715512) & (status_df.method == 'SpatialCellChat'))
        check('F1 I full-scale reference-failed marker reflects allocation failure',
              int(full_mask.sum()) == 1 and
              full_reference_row.status == 'allocation_failure_under_96_GiB_virtual_memory_budget' and
              not bool(full_reference_row.reference_result_generated))
        status_df.loc[full_mask, 'status'] = 'allocation_failure_under_memory_budget'
        status_df.loc[full_mask, 'record_type'] = 'reference_attempt_result'
        status_df.loc[full_mask, 'source_note'] = (
            'Attempted at 715,512 cells; failed requesting a 96.1 GiB dense intermediate '
            'under the 96 GiB virtual-memory limit; no reference result generated')
        status_df.loc[full_mask, 'display_label'] = 'reference failed'
    # Keep full-breast results in D/E and source data, but omit the unbalanced
    # near-170k column from the categorical completion panel.
    status_df['displayed_in_panel']=(status_df.N >= 1000) & \
        (~status_df.N.isin([169940,400000,443515]))
    planned_rows=[]
    for n in [200000,400000]:
        for method in ['SpatialCellChat','SpatialESS']:
            result = (EXTENSION[(EXTENSION.tissue == 'cervical') &
                                (pd.to_numeric(EXTENSION.N, errors='coerce') == n) &
                                (EXTENSION.engine == method)]
                      if len(EXTENSION) else pd.DataFrame())
            status = str(result.iloc[0].status) if len(result) else 'not_attempted'
            planned_rows.append(dict(dataset='Atera',tissue='cervical',N=n,
                method=method,status=status,
                record_type='extension_run_status' if len(result) else 'display_only_status',
                source_note='Cervical scale-extension run record' if len(result) else
                            'No validated extension run record is available',
                displayed_in_panel=n not in [400000,443515]))
    status_df=pd.concat([status_df,pd.DataFrame(planned_rows)],ignore_index=True)
    check('F1 I hides 250-cell and 169,940-cell categories only',
          int((~status_df.displayed_in_panel).sum()) == 8 and
          set(status_df.loc[~status_df.displayed_in_panel,'N']) ==
              {250,169940,400000,443515})
    row_order=[('CosMx','CancerousLiver','SpatialCellChat'),('CosMx','CancerousLiver','SpatialESS'),
               ('Atera','breast','SpatialCellChat'),('Atera','breast','SpatialESS'),
               ('Atera','cervical','SpatialCellChat'),('Atera','cervical','SpatialESS')]
    row_labels=['CosMx · SpatialCellChat','CosMx · SpatialESS',
                'Atera breast · SpatialCellChat','Atera breast · SpatialESS',
                'Atera cervical · SpatialCellChat','Atera cervical · SpatialESS']
    nominal_label={1000:'1k',5000:'5k',10000:'10k',25000:'25k',50000:'50k',
                   100000:'100k',200000:'200k',400000:'400k',443515:'450k',715512:'700k'}
    visible_sizes=sorted(status_df.loc[status_df.displayed_in_panel,'N'].unique())
    check('F1 I has requested categories and omits 170k',
          set(visible_sizes)=={1000,5000,10000,25000,50000,100000,200000,715512})
    cervical_added=status_df.query(
        "dataset == 'Atera' and tissue == 'cervical' and N in [200000,400000]")
    check('F1 I cervical 200k/400k records retained; 400k category hidden',
          len(cervical_added)==4 and
          not cervical_added.duplicated(['N','method']).any() and
          cervical_added.query('N == 200000').displayed_in_panel.all() and
          not cervical_added.query('N == 400000').displayed_in_panel.any())
    check('F1 I visible states need only completion and reference-failure symbols',
          status_df.loc[status_df.displayed_in_panel, 'status'].isin(
              ['completed', 'failed_signal_11', 'failed_allocation',
               'allocation_failure_under_memory_budget']).all())
    axis_labels=[nominal_label[int(n)] for n in visible_sizes]
    axis_position={int(n):i for i,n in enumerate(visible_sizes)}
    status_df['axis_category']=status_df.N.map(nominal_label)
    status_df['axis_position']=status_df.N.map(axis_position)
    ax=axes_mm(f,51,y+8,127,19);ax.set_xlim(-.5,len(visible_sizes)-.5)
    ax.set_ylim(-.55,len(row_order)-.45);ax.invert_yaxis()
    ax.set_xticks(range(len(visible_sizes)),axis_labels)
    ax.set_yticks(range(len(row_order)),row_labels,fontsize=5.0)
    ax.minorticks_off();ax.tick_params(axis='both',length=0,pad=1.2,labelsize=5.0)
    ax.grid(axis='x',which='major',color='#E4E9ED',lw=.35,zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():spine.set_visible(False)
    for yi,key in enumerate(row_order):
        d=status_df[(status_df.dataset==key[0])&(status_df.tissue==key[1])&
                    (status_df.method==key[2])]
        for _,r in d[d.displayed_in_panel].iterrows():
            actual_n=int(r.N);xval=r.axis_position;status=str(r.status)
            paired= status=='completed' and status_df[(status_df.dataset==key[0])&
                (status_df.tissue==key[1])&(status_df.N==actual_n)&
                (status_df.method!=key[2])].status.eq('completed').any()
            method_color=C['reference'] if key[2]=='SpatialCellChat' else C['ESS']
            if status=='completed':
                ax.plot(xval,yi,'o',ms=3.0,mfc=method_color if paired else 'white',
                        mec=method_color,mew=.7,zorder=3)
            elif status in ['failed_signal_11','allocation_failure_under_memory_budget'] or 'failed' in status:
                ax.plot(xval,yi,'x',ms=3.7,mew=.8,c=method_color,zorder=3)
            elif status=='not_attempted':
                ax.plot(xval,yi,'|',ms=4,mew=.75,c=C['secondary'],zorder=3)
            elif status.lower() in ['running','in_progress']:
                ax.plot(xval,yi,'D',ms=3.2,mfc='white',mec=method_color,mew=.8,zorder=3)
    check('Atera table included in figure 1 status',len(ATERA)>0 and status_df.dataset.eq('Atera').any())
    source('Figure1Atera_performance',ATERA,
           'source_data/current/atera_scale_performance.tsv;source_data/current/cervical_700k_reference_attempt.tsv',
           'tissue;fixture;N;engine;status;whole_process_seconds;peak_RSS_GB_decimal;failure_reason;inference_seconds',
           'tissue;fixture;engine','seconds;decimal GB','log10 x/y; failed and unattempted resources not plotted',
           'Measured completed runs only; Atera inputs use 18,028 genes, 3,081 fixed LR records, 100 permutations')
    source('Figure1I_completion',status_df,'current/cosmx_scaling.tsv;current/atera_scale_performance.tsv;current/cervical_extension_200k_400k.tsv;current/cervical_700k_reference_attempt.tsv',
           'dataset;tissue;N;method;status;record_type;source_note;display_label;displayed_in_panel;axis_category;axis_position',
           'dataset;tissue;N;method','completion state; record type; display flag; nominal input-size category',
           'equally spaced categorical x-axis; N retains the exact cell count',
           'Saved run states plus cervical 200k/400k extensions and the full 715,512-cell reference attempt. Panel I omits 250-cell CosMx, 169,940-cell breast, 400k cervical and 443,515-cell CosMx categories; all states remain in source data, and 443,515 remains in D/E. At 700k, the official run failed allocating a 96.1 GiB dense intermediate under the 96 GiB virtual-memory limit and uses the same reference-failed cross as earlier failed attempts. Filled circles paired complete; open circles SpatialESS-only; crosses failed reference runs.')
    source('Figure1C_Atera_agreement',ATERA_AGREE,'source_data/current/atera_scale_agreement.tsv',
           'tissue;N;active_jaccard;max_absolute_probability_error;exact_pvalue_fraction;significance_disagreements;strict_pass',
           'tissue;fixture','probability error; fractions; count','no transform',
           'All ten paired Atera breast/cervical cases, 1k–50k; thresholds as archived validation table')
    save(f,'Figure1_data');plt.close(f)

def fig2():
    f=figure('Figure2')
    ax=data_axis(f,'A','Primary LMM associations',left=11,right=31,top=7,bottom=10)
    volcano(ax)
    text_mm(f,13,47,f'{len(VALID):,} valid fits; {(VALID.q_value<.05).sum():,} BH q < 0.05',size=5.4)
    handles=[Line2D([],[],ls='none',marker='o',ms=2.8,c=c,label=l) for l,c in
             [('Control higher',C['negative']),('PF higher',C['positive']),('Not significant',C['nonsignificant'])]]
    f.legend(handles=handles,loc='center',bbox_to_anchor=(.254,1-7/196),ncol=3,fontsize=5.2)
    ax=data_axis(f,'B','Primary-cohort effect sizes',left=43,top=6,bottom=10)
    effect_forest(ax)
    text_mm(f,112,47,'● non-singular   ○ singular   |   labels: PF − Control × 0.001',size=5.4)
    effect_heatmap(f);patient_examples(f);bubble(f)
    save(f);plt.close(f)

def theme_preview():
    f=plt.figure(figsize=(180/25.4,67/25.4));f._name='Theme_preview'
    text_mm(f,2,2,'Representative data panels · 180 mm · Arial',size=8,weight='bold')
    a=axes_mm(f,13,12,42,42);scaling(a,'time')
    b=axes_mm(f,72,12,43,42);volcano(b,True)
    c=axes_mm(f,128,12,43,42)
    hm=HM[['standardized_effect','glm_standardized_effect']].values
    lim=np.max(abs(hm));c.imshow(hm[:12],cmap=DIVERGING,norm=effect_norm(lim),aspect='auto')
    c.set_xticks([0,1],['LMM','OLS']);c.set_yticks(range(12),HM.display_id[:12],size=5.6)
    c.tick_params(length=0)
    for s in c.spines.values():s.set_visible(False)
    save(f);plt.close(f)

def fig_s1():
    br=old('bridge_regression');spr=old('sparkle_probability_records');sm=old('sparkle_fidelity').iloc[0]
    tr=old('raw_corrected_transitions')
    workflow=pd.DataFrame({'step':[1,2,3],
      'label':['Corrected H5AD','Sparse matrix + metadata','ESS / V3 on identical input']})
    source('FigureS1A',workflow,'provenance/current_docs/RESULTS.md;source_data/archived/sparkle_fidelity.tsv.gz',
           'documented input bridge;method','step','workflow steps',
           rule='Conceptual summary of existing compatibility workflow',n_definition='steps')
    source('FigureS1B',br,'source_data/archived/bridge_regression.tsv.gz',
           'case;comparison;max_abs_probability_difference;pvalue_exact_fraction;active_jaccard;significance_disagreements',
           'case;comparison','metric-specific',rule='All 12 historical bridge/engine checks; keep comparison labels',version='Archived compatibility series')
    source('FigureS1C',spr,'source_data/archived/sparkle_probability_records.tsv.gz',
           'official_probability;spatialess_probability;official_pvalue;spatialess_pvalue',
           'sender_group_name;receiver_group_name;interaction_name','probability; permutation p',
           'log10 positive probabilities; count repeated p coordinates',version='Archived official paired corrected input')
    d=spr.groupby(['official_pvalue','spatialess_pvalue']).size().reset_index(name='records')
    check('SPARKLE p-value count conservation',d.records.sum()==len(spr))
    source('FigureS1D',d,'source_data/archived/sparkle_probability_records.tsv.gz',
           'official_pvalue;spatialess_pvalue','official_pvalue;spatialess_pvalue',
           'permutation p; records','area proportional to grouped record count',
           rule='Group identical saved p-value pairs; all records retained',n_definition='distinct p-value pairs')
    source('FigureS1F',tr,'source_data/archived/raw_corrected_transitions.tsv.gz','category;records',
           'category','significant records',rule='RAW/SPARKLE input sensitivity, not engine equivalence')
    f=figure('FigureS1');x,y,w,h=heading(f,'A','SPARKLE input bridge and validation')
    for xx,txt in [(4,'Corrected H5AD'),(50,'Sparse matrix + metadata'),(119,'ESS / V3 on identical input')]:
        text_mm(f,xx,y+8,txt,size=6.5)
    for a,b in [(37,47),(101,116)]:
        ax=axes_mm(f,a,y+6,b-a,7);ax.set_axis_off();ax.annotate('',(.95,.4),(.05,.4),
                  arrowprops=dict(arrowstyle='->',lw=.7,color=C['secondary']))
    x,y,w,h=heading(f,'B','Bridge and historical regression')
    ax=axes_mm(f,x,y+5,w,h-10)
    rows=[r['case'].replace('sparkle','SPK').replace('_',' ')+' / '+('E' if r.comparison=='old_vs_new' else 'B') for _,r in br.iterrows()]
    fields=['max_abs_probability_difference','pvalue_exact_fraction','active_jaccard','significance_disagreements']
    metric_table(ax,rows,['Max |ΔP|','Exact p','Jaccard','Δ calls'],
                 [[concise(r[k]) for k in fields] for _,r in br.iterrows()],widths=[2.6,1,1,1,1],
                 fontsize=5.7,fidelity_tiles=True)
    text_mm(f,x+1,y+h-3,'E: historical engine comparison; B: bridge / direct',size=5.5)
    ax=data_axis(f,'C','Corrected-input probabilities',left=12,right=16,top=11,bottom=10)
    mappable=agreement(ax,spr,'official_probability','spatialess_probability',
              (np.corrcoef(spr.official_probability,spr.spatialess_probability)[0,1],sm.max_abs_probability_difference),colorbar=False)
    marginals=probability_marginals(ax,spr,'official_probability','spatialess_probability')
    source('FigureS1C_marginals',marginals,'source_data/archived/sparkle_probability_records.tsv.gz',
           'official_probability;spatialess_probability','log10_bin_left','records',
           '24 equal log10 probability bins','All paired records',n_definition='display bins')
    pos=ax.get_position();fw,fh=f.get_size_inches()*25.4
    cax=f.add_axes([pos.x1+12/fw,pos.y0,1.4/fw,pos.height*.55])
    cb=f.colorbar(mappable,cax=cax);cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face')
    cb.outline.set_visible(False);cb.ax.tick_params(labelsize=5.2,width=.4,length=1.3)
    cb.ax.set_title('Records',fontsize=5.3,pad=3)
    ax=data_axis(f,'D','Permutation p-values',left=12,right=4,top=7,bottom=11)
    ax.plot([0,1],[0,1],ls='--',lw=.5,c=C['secondary'])
    area_per_record=90/d.records.max()
    ax.scatter(d.official_pvalue,d.spatialess_pvalue,s=area_per_record*d.records,c=C['ESS'],alpha=.85,edgecolors='none')
    d.to_csv(ROOT/'source_data/panels/FigureS1D_bins.csv',index=False)
    ax.set_aspect('equal');ax.set_xlim(-.04,1.04);ax.set_ylim(-.04,1.04);ax.set_xticks([0,.5,1]);ax.set_yticks([0,.5,1])
    ax.set_xlabel('SpatialCellChat V3 p-value');ax.set_ylabel('SpatialESS p-value')
    ax.text(.04,.96,'Exact fraction = '+concise(sm.pvalue_exact_fraction),transform=ax.transAxes,va='top',size=6)
    # An area key refers only to measured record counts.
    size_values=[200,1000,3000]
    for n in size_values:ax.scatter([],[],s=area_per_record*n,c=C['ESS'],label=f'{n:,}')
    ax.legend(title='Records',loc='lower right',fontsize=5.5,title_fontsize=5.5)
    x,y,w,h=heading(f,'E','Interaction-call agreement')
    ax=axes_mm(f,x+30,y+5,w-33,h-11)
    sg1=spr.official_pvalue<.05;sg2=spr.spatialess_pvalue<.05
    call_jaccard=(sg1&sg2).sum()/(sg1|sg2).sum()
    metrics=pd.DataFrame({'metric':['Active-record Jaccard','Exact p fraction',
                    'Significant-call Jaccard'],
                    'value':[sm.active_jaccard,sm.pvalue_exact_fraction,call_jaccard]})
    source('FigureS1E',metrics,'source_data/archived/sparkle_fidelity.tsv.gz;source_data/archived/sparkle_probability_records.tsv.gz',
           'active_jaccard;pvalue_exact_fraction;official_pvalue;spatialess_pvalue','metric','fraction',
           'saved values; calls use raw p < 0.05','Three existing agreement metrics')
    quantified_bars(ax,metrics.value,[v.replace('-','–') for v in metrics.metric],
                    [C['ESS'],C['reference'],C['positive']],formats=[f'{v:.4f}' for v in metrics.value],
                    fontsize=5.2,height=.42)
    ax.set_xlim(-.03,1.34);ax.set_xticks([0,.5,1]);ax.set_xlabel('Agreement fraction',fontsize=5.4,labelpad=.4)
    ax.tick_params(axis='x',labelsize=5.2,pad=1)
    ax=data_axis(f,'F','RAW / corrected input sensitivity',left=29,top=7,bottom=9)
    colors=[C['reference'],C['positive'],C['ESS']]
    quantified_bars(ax,tr.records,['Shared','Lost after SPARKLE','Gained after SPARKLE'],colors,
                    fontsize=5.6,height=.44)
    ax.set_xlabel('Significant records (raw p < 0.05)',fontsize=6)
    save(f);plt.close(f)

def cervical_fidelity():
    cv=old('cervical_fidelity').copy()
    cv['current_case']='cervical_'+cv['case']+'_'+cv['filter'].map({'default':'default','legacy':'relaxed','relaxed':'relaxed'})
    # Reuse real current comparison metrics after verifying the archived counts.
    for i,r in cv.iterrows():
        z=FID[FID['case']==r.current_case]
        check('Cervical current comparison '+r.current_case,len(z)==1)
        z=z.iloc[0]
        check('Cervical active identity '+r.current_case,r.ess_active==z.ESS_active and r.official_active==z.reference_active)
        cv.loc[i,'max_abs_probability_difference']=z.max_probability_difference
        cv.loc[i,'pvalue_exact_fraction']=z.exact_pvalue_fraction
        cv.loc[i,'active_jaccard']=z.active_jaccard
        cv.loc[i,'significance_disagreements']=z.significant_disagreements_lt005
    return cv

def fig_s2():
    cm=old('cervical_spatial_metadata');cv=cervical_fidelity()
    cr=pd.concat([read(f'cervical_{design}_{n}_default_paired.tsv',True)
                  for design in ['center','stratified'] for n in [1000,5000]],ignore_index=True)
    source('FigureS2A',cm,'source_data/archived/cervical_spatial_metadata.tsv.gz','cell_id;x_centroid;y_centroid;group',
           'cell_id','micrometres','equal aspect; original coordinate orientation','All 5k stratified cells; same five most frequent clusters + Other')
    source('FigureS2B',cv,'source_data/current/official_fidelity.tsv;source_data/archived/cervical_fidelity.tsv.gz',
           'case;filter;cells;max_abs_probability_difference;pvalue_exact_fraction;active_jaccard;significance_disagreements',
           'case;filter','metric-specific',rule='All eight matched configs; dimensions-only official comparator')
    source('FigureS2C',cr,'source_data/current/cervical_*_default_paired.tsv',
           'case;probability_official;probability_ess','case;sender_group_name;receiver_group_name;interaction_name',
           'probability','log10 positive probabilities',rule='Default filters; four full-key paired fixtures')
    f=figure('FigureS2')
    ax=data_axis(f,'A','Independent cervical Xenium',left=11,right=3,top=7,bottom=20)
    tops=cm.groupby('group').size().sort_values(ascending=False,kind='stable').head(5).index.tolist()
    cols=[C['ESS'],C['positive'],C['green'],C['purple'],'#555555','#BDC4CA']
    for group,color in zip(tops+['Other'],cols):
        d=cm[cm.group==group] if group!='Other' else cm[~cm.group.isin(tops)]
        ax.scatter(d.x_centroid,d.y_centroid,s=1.9,c=color,edgecolors='none',label=group.replace('Cluster_','C'))
    ax.set_aspect('equal',adjustable='box');ax.set_xlabel('x (μm)');ax.set_ylabel('y (μm)')
    ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.27),ncol=3,fontsize=5.5,markerscale=3,
              title='Public graph clusters · 5,000 stratified cells',title_fontsize=5.5)
    x,y,w,h=heading(f,'B','Fidelity across sampling designs')
    ax=axes_mm(f,x,y+6,w,h-15)
    rows=[('Central' if 'center' in r['case'] else 'Stratified')+' '+cell_label(r.cells)+
          ('  0.1 / 5' if r['filter']=='default' else '  0 / 1') for _,r in cv.iterrows()]
    fields=['max_abs_probability_difference','pvalue_exact_fraction','active_jaccard','significance_disagreements']
    metric_table(ax,rows,['Max |ΔP|','Exact p','Jaccard','Δ calls'],[[concise(r[k]) for k in fields] for _,r in cv.iterrows()],
                 widths=[2.65,1.3,1,1,1],fontsize=5.7,fidelity_tiles=True)
    text_mm(f,x,y+h-6,'V3: dimensions-only repair. Row suffix: min.percent / min.cells.sr.',size=5.6)
    x,y,w,h=heading(f,'C','Independent probability agreement')
    for i,(design,n) in enumerate([(d,n) for d in ['center','stratified'] for n in [1000,5000]]):
        xx=x+(i%2)*39;yy=y+6+(i//2)*39
        title=('Central' if design=='center' else 'Stratified')+' '+cell_label(n)
        text_mm(f,xx+21,yy,title,ha='center',size=6,weight='bold')
        ax=axes_mm(f,xx+9,yy+6,26,26)
        d=cr[cr['case']==f'cervical_{design}_{n}_default']
        lx,ly=np.log10(d.probability_official),np.log10(d.probability_ess)
        rg=[min(lx.min(),ly.min())-.1,max(lx.max(),ly.max())+.1]
        ax.plot(rg,rg,lw=.5,c=C['reference'],ls='--');ax.scatter(lx,ly,s=8,c=C['ESS'],alpha=.85,edgecolors='none')
        ax.set_xlim(rg);ax.set_ylim(rg);ax.set_aspect('equal');ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(3))
        marginals=probability_marginals(ax,d,'probability_official','probability_ess',
                                       marginal_height=2,marginal_width=3)
        source('FigureS2C_'+design+'_'+str(n)+'_marginals',marginals,
               'source_data/current/cervical_'+design+'_'+str(n)+'_default_paired.tsv',
               'probability_official;probability_ess','log10_bin_left','records',
               '24 equal log10 probability bins','All positive paired records',n_definition='display bins')
        ax.tick_params(labelsize=5.3);ax.text(.04,.93,f'n = {len(d)}',transform=ax.transAxes,size=5.5,va='top')
        if i//2:ax.set_xlabel('V3 (log10 P)',fontsize=5.6)
        if i%2==0:ax.set_ylabel('ESS (log10 P)',fontsize=5.6)
    x,y,w,h=heading(f,'D','Active and significant records')
    cases=cv['case'].drop_duplicates().tolist()
    long=[]
    for _,r in cv.iterrows():
        design='Central' if 'center' in r['case'] else 'Stratified'
        fixture=design+' '+cell_label(r.cells)+' '+('min.percent/min.cells.sr = 0.1/5' if r['filter']=='default' else 'min.percent/min.cells.sr = 0/1')
        for metric,prefix,suffix in [('Active records','active','_active'),
                                     ('Raw-p significant','significant','_significant_lt_005')]:
            long.extend([dict(fixture=fixture,metric=metric,method='SpatialESS',
                              value=r['ess'+suffix]),
                         dict(fixture=fixture,metric=metric,method='SpatialCellChat V3',
                              value=r['official'+suffix])])
    long=pd.DataFrame(long)
    source('FigureS2D',long,'source_data/archived/cervical_fidelity.tsv.gz',
           'case;filter;cells;ess_active;official_active;ess_significant_lt_005;official_significant_lt_005',
           'fixture;metric;method','records','grouped bars, paired methods',
           'All eight saved configurations; significance is raw permutation p < 0.05',
           n_definition='fixture by metric by method')
    fixtures=cv.apply(lambda r:('Central ' if 'center' in r['case'] else 'Strat. ')+cell_label(r.cells)+
                        (' 0.1/5' if r['filter']=='default' else ' 0/1'),axis=1).tolist()
    methods=['SpatialESS','SpatialCellChat V3']
    for j,metric in enumerate(['Active records','Raw-p significant']):
        ax=axes_mm(f,x+23+j*35,y+11,32,h-24)
        suffix='_active' if j==0 else '_significant_lt_005'
        maximum=cv[['ess'+suffix,'official'+suffix]].max().max()
        for i,(_,r) in enumerate(cv.iterrows()):
            if i%2==0:ax.axhspan(i-.46,i+.46,color='#F1F5F8',lw=0,zorder=0)
            for k,(prefix,c,m) in enumerate([('ess',C['ESS'],'o'),('official',C['reference'],'s')]):
                v=r[prefix+suffix];yy=i+(k-.5)*.32
                ax.barh(yy,v,height=.26,color=c,alpha=.42,lw=0)
                ax.plot(v,yy,marker=m,ms=2.5,c=c,ls='none')
                ax.text(maximum*1.29,yy,f'{v:,.0f}',ha='right',va='center',size=5.0,color=c)
        ax.set_yticks(range(len(fixtures)),fixtures if j==0 else ['']*len(fixtures),fontsize=5.0)
        ax.set_ylim(len(fixtures)-.5,-.5);ax.set_xlim(0,maximum*1.35)
        ax.set_xlabel('Records',fontsize=5.7)
        ax.yaxis.set_tick_params(length=0,pad=2);ax.tick_params(axis='x',labelsize=5.0)
        ax.xaxis.set_major_locator(MaxNLocator(3,integer=True));ax.spines['left'].set_visible(False)
        text_mm(f,x+23+j*35+16,y+6,'Active' if j==0 else 'Significant',ha='center',size=5.7,weight='bold')
    handles=[Line2D([],[],marker='o',ls='none',color=C['ESS'],ms=3,label='SpatialESS'),
             Line2D([],[],marker='o',ls='none',color=C['reference'],ms=3,label='SpatialCellChat V3')]
    f.legend(handles=handles,loc='center',bbox_to_anchor=((x+w*.62)/180,1-(y+h-1)/175),ncol=2,fontsize=5.2)
    text_mm(f,x+2,y+h-5,'Rows: min.percent / min.cells.sr. Significant: raw p < .05.',size=5.1)
    save(f);plt.close(f)

def fig_s3():
    f=figure('FigureS3')
    reps=SP.set_index('fixture').loc[['cosmx_50000','cosmx_443515','cervical_center_5000',
                                    'cervical_stratified_5000','ovarian_16247']].reset_index()
    source('FigureS3A',reps,'source_data/current/spatial_sparsity.tsv',
           'fixture;N;csr_total_bytes;dense_float64_bytes','fixture','decimal MB','log10 y',
           'Same five archived representative fixtures; CSR numeric-array payload; dense theoretical float64')
    source('FigureS3B',SP[SP.technology=='Xenium'],'source_data/current/spatial_sparsity.tsv',
           'fixture;N;mean_out_degree;sampling_design','fixture','neighbors per cell',
           rule='Four Xenium fixtures; existing graph metrics only')
    source('FigureS3C',SP,'source_data/current/spatial_sparsity.tsv',
           'fixture;N;isolated_cells','fixture','percent isolated','100 * isolated_cells / N',
           rule='All 14 graph fixtures; no graph recomputation')
    source('FigureS3D',SP[SP.series==SERIES[0]],'source_data/current/spatial_sparsity.tsv',
           'fixture;N;csr_total_bytes;all_numeric_array_bytes','fixture','decimal MB',
           'bytes / 1e6; log10 axes',rule='Nine CosMx fixtures; measured numeric arrays only')
    ax=data_axis(f,'A','Graph storage: dense versus CSR',left=12,top=8,bottom=14)
    for i,r in reps.iterrows():
        col=SCOL[SERIES.index(r.series)]
        ax.plot([i,i],[r.csr_total_bytes/1e6,r.dense_float64_bytes/1e6],c=col,lw=4.2,alpha=.32,solid_capstyle='butt')
        ax.plot([i,i],[r.csr_total_bytes/1e6,r.dense_float64_bytes/1e6],c=col,lw=.8)
        ax.plot(i,r.csr_total_bytes/1e6,'o',c=col,ms=4.8)
        ax.plot(i,r.dense_float64_bytes/1e6,'o',mfc='white',mec=col,ms=4.8,mew=.9)
    ax.set_yscale('log');ax.minorticks_off();ax.set_yticks([.01,1,100,10000,1e6]);ax.set_ylabel('Storage (MB; decimal)')
    ax.set_xticks(range(5),['CosMx\n50k','CosMx\n443k','Xenium\ncentral 5k','Xenium\nstratified 5k','Visium HD\n16,247'],fontsize=5.6)
    ax.text(.99,.98,'○ theoretical dense float64 (8N²)\n● measured CSR array payload',transform=ax.transAxes,va='top',ha='right',size=5.6)
    ax=data_axis(f,'B','Sampling geometry and degree',left=23,top=8,bottom=12)
    d=SP[SP.technology=='Xenium'].copy().sort_values(['series','N'])
    labels=[('Central ' if 'central' in r.series else 'Stratified ')+cell_label(r.N) for _,r in d.iterrows()]
    quantified_bars(ax,d.mean_out_degree,labels,[C['positive']]*len(d),
                    markers=['^' if 'central' in r.series else 'v' for _,r in d.iterrows()],
                    open_points=['stratified' in r.series for _,r in d.iterrows()],
                    formats=[f'{v:.3f}' if v<1 else f'{v:.2f}' for v in d.mean_out_degree],fontsize=5.5)
    ax.set_xlabel('Mean neighbors, E/N')
    ax=data_axis(f,'C','Isolated cells in tested graphs',left=29,top=8,bottom=10)
    rows=[];labels=[];colors=[];markers=[];hollow=[]
    for s,c,m in zip(SERIES,SCOL,SMARK):
        d=SP[SP.series==s].sort_values('N')
        prefix={'CosMx central':'CosMx', 'Xenium central':'X cen.',
                'Xenium stratified':'X strat.','Visium HD ovarian':'VHD'}[s]
        for _,r in d.iterrows():
            rows.append(100*r.isolated_cells/r.N);labels.append(prefix+' '+cell_label(r.N))
            colors.append(c);markers.append(m);hollow.append(s==SERIES[2])
    quantified_bars(ax,rows,labels,colors,markers,hollow,
                    formats=[f'{v:.3f}%' if 0<v<.01 else f'{v:.2f}%' for v in rows],fontsize=5.0)
    ax.set_xlabel('Isolated cells (%)',fontsize=6.1)
    ax=data_axis(f,'D','Stored arrays and extra graph fields',left=12,top=8,bottom=10)
    d=SP[SP.series==SERIES[0]].sort_values('N')
    ax.fill_between(d.N,d.csr_total_bytes/1e6,d.all_numeric_array_bytes/1e6,
                    color='#E4F0EC',lw=0,zorder=0)
    for field,label,color,marker in [('csr_total_bytes','Core CSR arrays',C['ESS'],'o'),
                                    ('all_numeric_array_bytes','All numeric arrays',C['green'],'s')]:
        d=SP[SP.series==SERIES[0]].sort_values('N')
        ax.plot(d.N,d[field]/1e6,marker=marker,ms=4.2,c=color,lw=1.2,label=label)
    log_cells(ax);ax.set_yscale('log');ax.minorticks_off();ax.set_ylabel('CosMx array payload (MB)')
    ax.legend(loc='upper left',fontsize=5.7)
    text_mm(f,3,122,'Array storage is not peak RSS. Core CSR = offsets + neighbors + distances. Shading shows the measured array-storage gap.',size=5.1)
    save(f);plt.close(f)

def fig_s4():
    f=figure('FigureS4');ms=old('model_summary').iloc[0];summary=old('lmm_glm_summary').iloc[0]
    status_rows=pd.DataFrame({'status':['Valid, non-singular','Valid, singular','Convergence warning','Numerical failure','Insufficient nonzero'],
                              'features':[ms.lmm_non_singular,ms.lmm_singular,ms.lmm_convergence_warning,ms.lmm_numerical_failure,ms.lmm_insufficient_nonzero]})
    source('FigureS4A',status_rows,'source_data/archived/model_summary.tsv.gz','lmm_* counts',
           'status','features',rule='Mutually exclusive documented fit states; no refitting')
    source('FigureS4B',old('lmm_glm_summary'),'source_data/archived/lmm_glm_summary.tsv.gz',
           'glm_significant;lmm_significant;overlap_significant','summary','features',rule='Existing BH q < .05 overlap')
    source('FigureS4C',BUBBLES,'source_data/archived/pathway_pair_representatives.tsv.gz',
           'feature_id;pair;pathway_name;standardized_effect;q_value','feature_id','effect / SD; feature-level q',
           'area min(-log10(q),10)',rule='Same 52 main-panel features; expanded labels only')
    source('FigureS4D',EX,'source_data/archived/patient_examples.tsv.gz','feature_id;estimate;lower;upper',
           'feature_id','condition contrast','effect ticks x10^-3',rule='Same three main examples, saved Wald intervals')
    ax=data_axis(f,'A','Model-fit diagnostics',left=34,top=8,bottom=10)
    quantified_bars(ax,status_rows.features,status_rows.status,
                    [C['ESS'],C['green'],C['positive'],C['purple'],C['nonsignificant']],fontsize=5.6)
    ax.set_xlabel('Features');ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    ax=data_axis(f,'B','LMM / slice-OLS discoveries',left=21,top=8,bottom=10)
    values=[summary.overlap_significant,summary.lmm_significant-summary.overlap_significant,
            summary.glm_significant-summary.overlap_significant]
    quantified_bars(ax,values,['Shared','LMM only','OLS only'],[C['reference'],C['ESS'],C['positive']],fontsize=5.6)
    ax.set_xlabel('Features with BH q < 0.05');ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    x,y,w,h=heading(f,'C','Expanded sender–receiver labels')
    ax=axes_mm(f,x+16,y+6,w-18,24)
    pairs=sorted(BUBBLES['pair'].unique());paths=sorted(BUBBLES.pathway_name.unique())
    px={v:i for i,v in enumerate(pairs)};py={v:i for i,v in enumerate(paths)}
    for i in range(len(paths)):
        if i%2==0:ax.axhspan(i-.45,i+.45,color='#F4F7F9',lw=0,zorder=0)
    ax.scatter(BUBBLES['pair'].map(px),BUBBLES.pathway_name.map(py),
       s=np.minimum(-np.log10(BUBBLES.q_value),10)*BUBBLE_AREA_PER_UNIT,c=BUBBLES.standardized_effect,
       cmap=DIVERGING,norm=effect_norm(EFFECT_LIMIT),edgecolors='none')
    ax.set_yticks(range(len(paths)),paths,fontsize=5.5);ax.set_xticks(range(len(pairs)),[f'P{i+1}' for i in range(len(pairs))],fontsize=5.5)
    ax.set_xlim(-.5,len(pairs)-.5);ax.set_ylim(-.5,len(paths)-.5);ax.tick_params(length=0,pad=1)
    for sp in ax.spines.values():sp.set_visible(False)
    for i,p in enumerate(pairs):
        a=BUBBLES[BUBBLES['pair']==p].iloc[0]
        text_mm(f,x+1,y+34+i*3.3,f'P{i+1}  {abbrev(a.sender_group_name)} → {abbrev(a.receiver_group_name)}',size=5.5)
    ax=data_axis(f,'D','Condition effects for patient examples',left=31,top=9,bottom=16)
    for i,(_,r) in enumerate(EX.iterrows()):
        if i%2==0:ax.axhspan(i-.43,i+.43,color='#F1F5F8',lw=0,zorder=0)
        ax.plot([0,r.estimate],[i,i],color='#E8A7B1',lw=4.2,solid_capstyle='butt',zorder=1)
        ax.plot([r.lower,r.upper],[i,i],color=C['positive'],lw=1)
        ax.plot([r.lower,r.lower],[i-.12,i+.12],color=C['positive'],lw=.7)
        ax.plot([r.upper,r.upper],[i-.12,i+.12],color=C['positive'],lw=.7)
        ax.plot(r.estimate,i,'o',c=C['positive'],ms=4.8)
    ax.axvline(0,lw=.5,c=C['reference']);ax.set_yticks(range(3),[r.display_id+'\n'+lr_label(r) for _,r in EX.iterrows()],fontsize=6)
    ax.set_ylim(2.5,-.5);ax.set_xlabel('PF − Control (× 0.001; 95% Wald CI)',fontsize=6)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x*1000:g}'));ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False)
    text_mm(f,x+94,y+h-8,'Intervals describe the condition contrast,\nnot adjusted group means.',size=5.8)
    save(f);plt.close(f)

if args.figure1_only:
    fig1()
    pd.DataFrame(PROV).to_csv(ROOT/'reports/PANEL_INVENTORY.tsv',sep='\t',index=False)
    pd.DataFrame(CHECKS).to_csv(ROOT/'reports/data_checks.tsv',sep='\t',index=False)
    with pd.ExcelWriter(ROOT/'source_data/Figure1_Source_Data.xlsx',engine='openpyxl') as writer:
        pd.DataFrame(PROV).to_excel(writer,sheet_name='Provenance',index=False)
        for path in sorted((ROOT/'source_data/panels').glob('Figure1*.csv')):
            pd.read_csv(path,float_precision='round_trip').to_excel(writer,sheet_name=path.stem[:31],index=False)
    (ROOT/'logs/environment.json').write_text(json.dumps(dict(
        python=sys.version,platform=platform.platform(),numpy=np.__version__,pandas=pd.__version__,
        matplotlib=mpl.__version__,font=FONT,font_path=FONT_PATH),indent=2)+'\n')
    print('Figure 1 label revision exported; PDF composition and visual QA still required.',flush=True)
    sys.exit(0)

theme_preview()
if args.theme_only:
    print('Theme-only preview completed.',flush=True)
    sys.exit(0)
fig1();fig2()
if args.main_only:
    pd.DataFrame(PROV).to_csv(ROOT/'reports/PANEL_INVENTORY.tsv',sep='\t',index=False)
    pd.DataFrame(CHECKS).to_csv(ROOT/'reports/data_checks.tsv',sep='\t',index=False)
    print('Main-only design preview completed.',flush=True)
    sys.exit(0)
fig_s1();fig_s2();fig_s3();fig_s4()
pd.DataFrame(PROV).to_csv(ROOT/'reports/PANEL_INVENTORY.tsv',sep='\t',index=False)
pd.DataFrame(CHECKS).to_csv(ROOT/'reports/data_checks.tsv',sep='\t',index=False)
pd.DataFrame([dict(semantic=k,hex=v) for k,v in C.items()]).to_csv(ROOT/'config/palette.tsv',sep='\t',index=False)
with pd.ExcelWriter(ROOT/'source_data/Source_Data.xlsx',engine='openpyxl') as writer:
    pd.DataFrame(PROV).to_excel(writer,sheet_name='Provenance',index=False)
    mapping.to_excel(writer,sheet_name='Feature_IDs',index=False)
    for p in sorted((ROOT/'source_data/panels').glob('*.csv')):
        if p.stem not in {'FigureS1CD','Figure2D_effects'}:
            pd.read_csv(p).to_excel(writer,sheet_name=p.stem[:31],index=False)
    FID.to_excel(writer,sheet_name='All_66_official_comparisons',index=False)
    pd.read_csv(ROOT/'source_data/pair_mapping.csv').to_excel(writer,sheet_name='Sender_receiver_IDs',index=False)
(ROOT/'logs/environment.json').write_text(json.dumps(dict(
    python=sys.version,platform=platform.platform(),numpy=np.__version__,pandas=pd.__version__,
    matplotlib=mpl.__version__,font=FONT,font_path=FONT_PATH),indent=2)+'\n')
print('All data panels exported; composition and visual QA still required.',flush=True)
