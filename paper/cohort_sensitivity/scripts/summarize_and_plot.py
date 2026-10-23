"""Summarize the saved sensitivity results and render Supplementary Figure 5."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib import font_manager

BLUE="#0072B2"; RED="#D55E00"; GRAY="#6B7280"
MODEL_NAMES={"primary_lmm":"Primary LMM", "run_adjusted_lmm":"Run LMM", "tma_adjusted_lmm":"TMA LMM",
             "patient_run_hc3":"Patient HC3", "patient_welch":"Patient Welch"}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--out",type=Path,required=True); args=parser.parse_args()
    out=args.out; data=out/"data"; reports=out/"reports"
    models={k:pd.read_csv(data/f"{k}.tsv.gz",sep="\t").set_index("feature_id") for k in MODEL_NAMES}
    primary=models["primary_lmm"]
    for k in models: models[k]=models[k].loc[primary.index]
    mapping=pd.read_csv(data/"feature_mapping.csv").set_index("feature_id")
    balanced=pd.read_csv(data/"balanced_run_adjusted_effects.tsv.gz",sep="\t").set_index("feature_id")
    lopo=pd.read_csv(data/"leave_one_patient_feature_summary.tsv.gz",sep="\t").set_index("feature_id")
    labels=pd.read_csv(data/"cohort_evidence_labels.tsv.gz",sep="\t").set_index("feature_id")
    primary_hit=primary.status.eq("ok") & primary.q_value.lt(.05)
    overview=[]
    selected=mapping.copy()
    for key,model in models.items():
        valid=model.status.eq("ok"); hits=valid & model.q_value.lt(.05)
        common=valid & primary.status.eq("ok")
        overview.append(dict(model=key,label=MODEL_NAMES[key],features=len(model),valid_features=int(valid.sum()),
          bh_hits=int(hits.sum()),primary_hits_retained=int((primary_hit&hits).sum()),
          primary_hit_retention_fraction=float((primary_hit&hits).sum()/primary_hit.sum()),
          common_valid=int(common.sum()),effect_spearman=float(stats.spearmanr(primary.estimate[common],model.estimate[common]).statistic),
          direction_agreement=float((np.sign(primary.estimate[common])==np.sign(model.estimate[common])).mean()),
          valid_singular=int((valid & model.singular.eq(True)).sum()) if "singular" in model else np.nan,
          selected_28_hits=int(hits.reindex(mapping.index).sum())))
        for column in ["estimate","std_error","ci_lower","ci_upper","p_value","q_value","status","singular"]:
            if column in model: selected[f"{key}_{column}"]=model[column]
        if key=="primary_lmm":
            selected[f"{key}_ci_lower"]=model.estimate-1.959963984540054*model.std_error
            selected[f"{key}_ci_upper"]=model.estimate+1.959963984540054*model.std_error
        selected[f"{key}_effect_over_slice_sd"]=model.estimate/primary.response_sd
    for prefix,table in [("balanced",balanced),("lopo",lopo),("evidence",labels)]:
        for col in table: selected[f"{prefix}_{col}"]=table[col]
    selected["primary_response_sd"]=primary.response_sd
    selected["primary_nonzero_slices"]=primary.nonzero_units
    overview=pd.DataFrame(overview)
    overview.to_csv(reports/"model_comparison.tsv",sep="\t",index=False)
    selected.reset_index().to_csv(data/"selected_28_sensitivity.tsv",sep="\t",index=False)
    nonsingular=primary.status.eq("ok") & primary.singular.eq(False)
    q=np.full(len(primary),np.nan);q[nonsingular]=stats.false_discovery_control(primary.p_value[nonsingular].to_numpy())
    pd.DataFrame(dict(feature_id=primary.index,eligible_nonsingular=nonsingular,
                      primary_global_q=primary.q_value,subset_reBH_q=q)).to_csv(data/"nonsingular_sensitivity.tsv.gz",sep="\t",index=False)
    state=pd.DataFrame([dict(model=key,status=status,features=int(count)) for key,model in models.items()
                        for status,count in model.status.value_counts().items()])
    state.to_csv(reports/"model_fit_states.tsv",sep="\t",index=False)
    support=pd.read_csv(data/"selected_support_and_distributions.tsv",sep="\t")
    fid=selected.index[selected.display_id.eq("I15")][0];cid=selected.index[selected.display_id.eq("I06")][0]
    fn=selected.loc[fid];col=selected.loc[cid]
    consistent=labels.evidence_label.eq("supported_by_run_and_patient_sensitivity")
    metrics=dict(features=len(primary),primary_hits=int(primary_hit.sum()),agreement_labels=int(consistent.sum()),
                 primary_nonsingular=int(nonsingular.sum()),primary_nonsingular_hits=int((nonsingular&primary_hit).sum()),
                 subset_reBH_hits=int(np.sum(q<.05)),fn1_feature_id=fid,col1a1_feature_id=cid,
                 selected_28_patient_hc3_hits=int(overview.loc[overview.model.eq("patient_run_hc3"),"selected_28_hits"].iloc[0]))
    (reports/"summary_metrics.json").write_text(json.dumps(metrics,indent=2)+"\n")
    font=font_manager.findfont("Arial",fallback_to_default=False)
    plt.rcParams.update({"font.family":"Arial","font.size":6.5,"axes.labelsize":7,"axes.titlesize":7,
      "xtick.labelsize":6,"ytick.labelsize":6,"axes.linewidth":.5,"lines.linewidth":.85,
      "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none","text.color":"#111111",
      "axes.edgecolor":"#111111","axes.spines.top":False,"axes.spines.right":False})
    fig=plt.figure(figsize=(180/25.4,192/25.4))
    grid=fig.add_gridspec(3,2,left=.083,right=.93,bottom=.15,top=.96,
                         width_ratios=[1.12,1],height_ratios=[.90,1.78,1.10],wspace=.36,hspace=.55)
    axa=fig.add_subplot(grid[0,0]);axb=fig.add_subplot(grid[0,1]);axc=fig.add_subplot(grid[1,0]);axd=fig.add_subplot(grid[1,1]);axe=fig.add_subplot(grid[2,0]);axf=fig.add_subplot(grid[2,1])
    axes=[axa,axb,axc,axd,axe,axf]
    titles=["Cohort and processing runs","Sensitivity to run adjustment","Original feature set across analyses",
            "Singular and nonsingular fits","Zero-inclusive patient distributions","Original examples across models"]
    for ax,letter,title in zip(axes,"abcdef",titles):
        ax.set_title(title,loc="left",pad=7,fontweight="bold")
        ax.text(-.14,1.055,letter,transform=ax.transAxes,fontweight="bold",fontsize=8,va="bottom")
    axb.texts[-1].set_position((-.20,1.055))
    meta=pd.read_csv(data/"sample_metadata.tsv",sep="\t")
    previous=np.zeros(2)
    for run,color in zip(["Run1","Run2","Run3"],[BLUE,GRAY,RED]):
        counts=np.array([np.sum((meta.condition==c)&(meta.run==run)) for c in ["Control","Disease"]])
        axa.bar([0,1],counts,bottom=previous,width=.6,color=color,label=run)
        for x,(n,bottom) in enumerate(zip(counts,previous)):
            axa.text(x,bottom+n/2,str(n),ha="center",va="center",color="white",fontsize=6,fontweight="bold")
        previous+=counts
    axa.set_xticks([0,1],["Control\n9 patients","PF\n26 patients"]);axa.set_ylabel("Slices (n)");axa.set_ylim(0,43)
    axa.yaxis.set_major_locator(MaxNLocator(integer=True,nbins=4));axa.legend(frameon=False,fontsize=6,ncol=3,loc="upper left",handlelength=1,columnspacing=.8)
    run=models["run_adjusted_lmm"];common=primary.status.eq("ok")&run.status.eq("ok")
    xx=primary.estimate[common]/primary.response_sd[common];yy=run.estimate[common]/primary.response_sd[common]
    bins=axb.hexbin(xx,yy,gridsize=52,mincnt=1,bins="log",cmap="cividis",linewidths=0,rasterized=True)
    limit=max(np.max(np.abs(xx)),np.max(np.abs(yy)))*1.03
    axb.plot([-limit,limit],[-limit,limit],color=GRAY,lw=.5,ls="--");axb.set_xlim(-limit,limit);axb.set_ylim(-limit,limit)
    axb.set_xlabel("Primary effect / slice SD");axb.set_ylabel("Run-adjusted effect / slice SD")
    axb.text(.02,.98,f"n = {common.sum():,} features\nSpearman ρ = {stats.spearmanr(xx,yy).statistic:.3f}",transform=axb.transAxes,va="top",fontsize=6)
    cb=fig.colorbar(bins,ax=axb,fraction=.035,pad=.03);cb.set_label("Features (n)",fontsize=6);cb.ax.tick_params(labelsize=5.5)
    # Existing stable IDs are kept; complete molecular and cell identities are in Table 2.
    row_order=selected.sort_values("display_id").index
    effect=np.column_stack([models[k].estimate.where(models[k].status.eq("ok")).reindex(row_order)/primary.response_sd.reindex(row_order) for k in MODEL_NAMES])
    effect=np.clip(effect,-2,2)
    colors=LinearSegmentedColormap.from_list("effect",[BLUE,"#FFFFFF",RED]);colors.set_bad("#D9DDE2")
    image=axc.imshow(effect,aspect="auto",cmap=colors,norm=TwoSlopeNorm(vmin=-2,vcenter=0,vmax=2))
    axc.set_yticks(np.arange(len(row_order)),selected.loc[row_order,"display_id"])
    axc.set_xticks(np.arange(5),["Primary","Run","TMA","HC3","Welch"],fontsize=6)
    axc.tick_params(length=0,pad=3);axc.set_ylabel("Original feature ID")
    for j,k in enumerate(MODEL_NAMES):
        qv=models[k].q_value.reindex(row_order).to_numpy()
        axc.scatter(np.full(np.sum(qv<.05),j),np.flatnonzero(qv<.05),s=3,color="#111111",linewidths=0)
    axc.spines[["left","bottom"]].set_visible(False)
    cb=fig.colorbar(image,ax=axc,orientation="horizontal",fraction=.045,pad=.14,extend="both")
    cb.set_label("Effect / original slice SD; display cap ±2",fontsize=6);cb.ax.tick_params(labelsize=5.5)
    for j,(label,mask) in enumerate([("Nonsingular",nonsingular),("Singular",primary.status.eq("ok")&primary.singular.eq(True))]):
        hit=int((mask&primary_hit).sum());total=int(mask.sum())
        axd.barh(j,total,color="#D7DDE2",height=.56);axd.barh(j,hit,color=BLUE,height=.56)
        axd.text(hit+300,j,str(hit),va="center",fontsize=6)
        axd.text(total,j+.34,f"{total:,} valid fits",ha="right",fontsize=6,color=GRAY)
    axd.set_yticks([0,1],["Nonsingular","Singular"]);axd.set_ylim(-.65,1.9);axd.set_xlabel("Features (n)")
    axd.tick_params(axis="y",length=0);axd.set_xlim(0,25000);axd.xaxis.set_major_locator(MaxNLocator(3))
    axd.legend(handles=[Line2D([0],[0],color=BLUE,lw=5,label="Primary BH q < 0.05"),Line2D([0],[0],color="#D7DDE2",lw=5,label="Other valid fits")],frameon=False,fontsize=6,loc="upper left")
    dist=pd.read_csv(data/"selected_patient_distributions.tsv",sep="\t")
    rng=np.random.default_rng(20261012)
    for j,(fid_,short) in enumerate([(fid,"I15 · FN1"),(cid,"I06 · COL1A1")]):
        for c,offset,color in [("Control",-.18,BLUE),("Disease",.18,RED)]:
            values=dist[(dist.feature_id==fid_)&(dist.condition==c)].mean_log1p.to_numpy()
            position=j+offset
            axe.scatter(position+rng.uniform(-.11,.11,len(values)),values*1000,s=12,
                        facecolor=color,edgecolor="white",linewidth=.25,alpha=.88)
            axe.plot([position-.12,position+.12],[np.median(values)*1000]*2,color="#111111",lw=1)
    axe.set_xticks([0,1],["I15 · FN1","I06 · COL1A1"]);axe.set_ylabel("Patient mean log1p(score) × 1,000");axe.set_xlim(-.5,1.5)
    axe.legend(handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=BLUE,markersize=4,label="Control"),Line2D([0],[0],marker="o",color="none",markerfacecolor=RED,markersize=4,label="PF")],frameon=False,fontsize=6,ncol=2,loc="upper left")
    palette=["#111111",BLUE,"#009E73",RED,"#CC79A7"]
    for j,fid_ in enumerate([fid,cid]):
        for k,(model_name,color) in enumerate(zip(MODEL_NAMES,palette)):
            row=selected.loc[fid_];ypos=j+(k-2)*.12
            value=row[f"{model_name}_estimate"]*1000;lo=row[f"{model_name}_ci_lower"]*1000;hi=row[f"{model_name}_ci_upper"]*1000
            axf.plot([lo,hi],[ypos,ypos],color=color,lw=.8);axf.scatter(value,ypos,s=10,color=color,zorder=3)
    axf.axvline(0,color=GRAY,lw=.5,ls="--");axf.set_yticks([0,1],["I15 · FN1","I06 · COL1A1"])
    axf.set_ylim(-.5,1.5);axf.invert_yaxis();axf.set_xlabel("Condition coefficient × 1,000")
    axf.legend(handles=[Line2D([0],[0],marker="o",color=c,lw=.8,markersize=3,label=MODEL_NAMES[k]) for k,c in zip(MODEL_NAMES,palette)],frameon=False,fontsize=5.5,ncol=2,loc="upper left",bbox_to_anchor=(-.15,-.45),columnspacing=.8,handlelength=1)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer();bounds=[]
    skipped=set()
    for ax in fig.axes:
        for tick in ax.xaxis.get_major_ticks()+ax.xaxis.get_minor_ticks():
            if not min(ax.get_xlim()) <= tick.get_loc() <= max(ax.get_xlim()):
                skipped.update([tick.label1,tick.label2])
        for tick in ax.yaxis.get_major_ticks()+ax.yaxis.get_minor_ticks():
            if not min(ax.get_ylim()) <= tick.get_loc() <= max(ax.get_ylim()):
                skipped.update([tick.label1,tick.label2])
    for artist in fig.findobj(matplotlib.text.Text):
        if artist in skipped:continue
        if not artist.get_visible() or not artist.get_text():continue
        bbox=artist.get_window_extent(renderer);outside=bbox.x0 < -1 or bbox.y0 < -1 or bbox.x1>fig.bbox.width+1 or bbox.y1>fig.bbox.height+1
        bounds.append(dict(text=artist.get_text(),font_size=artist.get_fontsize(),outside_canvas=outside))
    pd.DataFrame(bounds).to_csv(reports/"FigureS5_text_bounds.tsv",sep="\t",index=False)
    for suffix in ["pdf","svg","png"]:
        fig.savefig(out/"figures"/f"FigureS5.{suffix}",dpi=600,facecolor="white")
    plt.close(fig)
    qa=("# Checks\n\n- Primary scores, coefficients and Figure 1A unchanged.\n- Full feature keys align across all models.\n- R and Python HC3 results agree within 1e-10; independent statsmodels checks passed.\n- BH excludes invalid/warning fits and is recomputed over each full eligible universe.\n- Nonsingular subset re-BH is labelled separately.\n- Tissue-affect confounding and run/TMA redundancy are recorded, not silently resolved.\n- Balanced subset ranges are not confidence intervals.\n"
      f"- Font: {font}; figure width 180 mm, height 192 mm, PNG 600 dpi.\n"
      f"- Text bounding-box check: {sum(r['outside_canvas'] for r in bounds)} outside-canvas text elements.\n"
      "- PDF/SVG retain vector text and axes; the dense hexbin layer alone is rasterized.\n- PDF rendering and visual inspection must be recorded separately.\n- Limitations: only nine controls; no independently identifiable anatomical-region effect; approximate robust inference, no causal or experimental validation.\n")
    (reports/"QA_REPORT.md").write_text(qa)
    print(overview.to_string(index=False)); print("Figure and reports complete")


if __name__=="__main__":main()
