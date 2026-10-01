"""Render current resource and sparsity figures only from distributed TSVs."""
import csv
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root=Path(sys.argv[1] if len(sys.argv)>1 else '.').resolve()
out=root/'figures';out.mkdir(exist_ok=True)
def read(name):
    with (root/'results/tables'/name).open() as f:
        return list(csv.DictReader(f,delimiter='\t'))
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
                     'font.size':9,'axes.labelsize':9,'axes.titlesize':10,
                     'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
                     'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.6})
blue='#2468A2';red='#C74945';neutral='#696969'
def finish(fig,name):
    for ext in ['pdf','svg','png']:
        fig.savefig(out/(name+'.'+ext),dpi=600 if ext=='png' else None,bbox_inches='tight',facecolor='white')
    plt.close(fig)
rows=read('cosmx_scaling.tsv')
fig,axes=plt.subplots(1,2,figsize=(7.1,2.9),layout='constrained')
for ax,metric,official,ylabel,label in zip(axes,['wall_seconds','peak_rss_GiB'],
        ['official_wall_seconds','official_peak_rss_gib'],['Wall time (s)','Peak RSS (GiB)'],['A','B']):
    n=np.array([int(r['cells']) for r in rows])
    v=np.array([float(r['ESS_'+metric]) for r in rows])
    lo=np.array([float(r['ESS_'+metric+'_min']) for r in rows])
    hi=np.array([float(r['ESS_'+metric+'_max']) for r in rows])
    ax.errorbar(n,v,yerr=np.array([v-lo,hi-v]),fmt='o-',color=blue,ms=3.5,lw=1.1,capsize=2,label='SpatialESS 0.1.4')
    complete=[r for r in rows if r['official_status']=='completed' and r[official] not in ('NA','')]
    ax.plot([int(r['cells']) for r in complete],[float(r[official]) for r in complete],
            's-',color=neutral,ms=3.5,lw=1.1,label='SpatialCellChat V3 reference')
    ax.set(xscale='log',yscale='log',xlabel='Number of cells',ylabel=ylabel)
    ax.text(-.15,1.05,label,transform=ax.transAxes,fontweight='bold',fontsize=12)
    ax.grid(axis='y',alpha=.13,linewidth=.5)
axes[0].legend(frameon=False,fontsize=7,loc='upper left')
finish(fig,'Benchmark_resources')

rows=read('spatial_sparsity.tsv')
fig,axes=plt.subplots(1,3,figsize=(9,2.8),layout='constrained')
styles={'CosMx':(blue,'o'),'Xenium':(red,'s'),'Visium HD':('#547E53','^')}
for technology in sorted({r['technology'] for r in rows}):
    selected=[r for r in rows if r['technology']==technology]
    color,marker=styles.get(technology,(red,'D'))
    selected.sort(key=lambda r:int(r['N']))
    for ax,field in zip(axes,['E','mean_out_degree','graph_density']):
        ax.scatter([int(r['N']) for r in selected],[float(r[field]) for r in selected],
                   c=color,marker=marker,s=25,label=technology,edgecolors='white',linewidth=.3,zorder=3)
        if technology=='CosMx':
            ax.plot([int(r['N']) for r in selected],[float(r[field]) for r in selected],c=color,lw=.8)
n=np.array(sorted({int(r['N']) for r in rows}),dtype=float)
axes[0].plot(n,n*n,'--',color=neutral,lw=.8,label='N squared')
for ax,ylabel,label in zip(axes,['Directed edges / possible pairs','Mean outgoing degree','Graph density (E / N squared)'],['A','B','C']):
    ax.set(xscale='log',yscale='log',xlabel='Number of cells',ylabel=ylabel)
    ax.text(-.15,1.05,label,transform=ax.transAxes,fontweight='bold',fontsize=12)
    ax.grid(axis='y',alpha=.13,linewidth=.5)
axes[0].legend(frameon=False,fontsize=7)
finish(fig,'Spatial_neighborhoods')
(out/'README.md').write_text('# Result figures\n\n'
 'Benchmark_resources: current SpatialESS original-harness measurements; official measurements '
 'are from the identified reference experiment. Only completed observations are plotted. '
 'Failure status and reason remain in cosmx_scaling.tsv; no failed time or RSS is plotted.\n\n'
 'Spatial_neighborhoods: measured directed graph edges, mean outgoing degree and graph density '
 'under the tested settings. The dashed pair-space reference is N squared.\n\n'
 'All plotted values are read from results/tables/cosmx_scaling.tsv and '
 'results/tables/spatial_sparsity.tsv. Text uses Arial with Helvetica/DejaVu Sans fallback.\n')
print('Rendered figures from current tables.')
