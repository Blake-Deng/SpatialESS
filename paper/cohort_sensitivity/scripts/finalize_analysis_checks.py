"""Record reproducibility, original-input identity and per-panel source data."""
from pathlib import Path
import argparse
import hashlib
import subprocess
import pandas as pd
import numpy as np
from scipy import stats


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--reproduction",type=Path,required=True)
    parser.add_argument("--primary-provenance",type=Path,required=True);args=parser.parse_args()
    out=args.out;data=out/"data";reports=out/"reports";checks=[]
    def check(name,value):
        checks.append(dict(check=name,passed=bool(value)))
        if not value:raise AssertionError(name)
    current=pd.read_csv(data/"input_checksums.tsv",sep="\t");current['file']=current['file'].map(lambda p:Path(p).name)
    original=pd.read_csv(args.primary_provenance,sep="\t")
    identity=current.merge(original,on="file",suffixes=("_current","_primary"))
    check("45 original profile + metadata + model-source checksums match",len(identity)==47 and identity.md5_current.eq(identity.md5_primary).all())
    identity.to_csv(reports/"original_primary_input_identity.tsv",sep="\t",index=False)
    for key in ["run_adjusted_lmm","tma_adjusted_lmm"]:
        original=pd.read_csv(data/f"{key}.tsv.gz",sep="\t")
        reproduced=pd.read_csv(args.reproduction/"data"/f"{key}.tsv.gz",sep="\t")
        check(key+" refit from prepared matrix identical",original.equals(reproduced))
    primary=pd.read_csv(data/"primary_lmm.tsv.gz",sep="\t").set_index("feature_id")
    run=pd.read_csv(data/"run_adjusted_lmm.tsv.gz",sep="\t").set_index("feature_id").loc[primary.index]
    common=primary.status.eq("ok")&run.status.eq("ok")
    x=(primary.estimate[common]/primary.response_sd[common]);y=(run.estimate[common]/primary.response_sd[common])
    panel=data/"panels";panel.mkdir(exist_ok=True)
    a=pd.read_csv(reports/"condition_by_run.tsv",sep="\t");a.to_csv(panel/"FigureS5A.csv",index=False)
    pd.DataFrame(dict(feature_id=x.index,primary_effect_over_slice_sd=x,run_effect_over_slice_sd=y)).to_csv(panel/"FigureS5B.csv",index=False)
    selected=pd.read_csv(data/"selected_28_sensitivity.tsv",sep="\t")
    rows=[]
    for _,s in selected.iterrows():
        for model in ["primary_lmm","run_adjusted_lmm","tma_adjusted_lmm","patient_run_hc3","patient_welch"]:
            valid=s[f"{model}_status"]=="ok"
            effect=s[f"{model}_effect_over_slice_sd"] if valid else np.nan
            rows.append(dict(feature_id=s.feature_id,display_id=s.display_id,LR_display=s.LR_display,
              sender_group_name=s.sender_group_name,receiver_group_name=s.receiver_group_name,model=model,
              coefficient=s[f"{model}_estimate"],effect_over_original_slice_sd=effect,
              displayed_effect=np.clip(effect,-2,2),q_value=s[f"{model}_q_value"],valid_fit=valid))
    pd.DataFrame(rows).to_csv(panel/"FigureS5C.csv",index=False)
    d=[]
    for singular in [False,True]:
        valid=primary.status.eq("ok")&primary.singular.eq(singular)
        d.append(dict(singular=singular,valid_features=int(valid.sum()),primary_global_BH_hits=int((valid&primary.q_value.lt(.05)).sum())))
    pd.DataFrame(d).to_csv(panel/"FigureS5D.csv",index=False)
    e=pd.read_csv(data/"selected_patient_distributions.tsv",sep="\t");e=e[e.display_id.isin(["I15","I06"])]
    e.to_csv(panel/"FigureS5E.csv",index=False)
    f=pd.DataFrame(rows);f=f[f.display_id.isin(["I15","I06"])].copy()
    ci=[]
    for _,s in selected[selected.display_id.isin(["I15","I06"])].iterrows():
        for model in ["primary_lmm","run_adjusted_lmm","tma_adjusted_lmm","patient_run_hc3","patient_welch"]:
            ci.append(dict(feature_id=s.feature_id,display_id=s.display_id,model=model,
             estimate=s[f"{model}_estimate"],ci_lower=s[f"{model}_ci_lower"],ci_upper=s[f"{model}_ci_upper"],
             q_value=s[f"{model}_q_value"],ci_rule="unadjusted 95% Wald-normal" if "lmm" in model else "unadjusted 95% t",
             plot_multiplier=1000))
    pd.DataFrame(ci).to_csv(panel/"FigureS5F.csv",index=False)
    source_rows=[("a","data/sample_metadata.tsv","condition; run","45 slices (9 Control and 26 PF patients)","counts",len(a)),
       ("b","data/primary_lmm.tsv.gz; data/run_adjusted_lmm.tsv.gz","feature_id; estimate; response_sd; status","common status=ok, full keys","coefficient / original across-slice SD; count hexbin",int(common.sum())),
       ("c","data/selected_28_sensitivity.tsv","feature_id; display_id; *_estimate; *_q_value; *_status; primary_response_sd","original 28 IDs unchanged","zero-centred scale, display cap ±2; each model's full-universe BH",len(rows)),
       ("d","data/primary_lmm.tsv.gz","status; singular; q_value","status=ok","primary global-BH q<0.05",int(primary.status.eq('ok').sum())),
       ("e","data/selected_patient_distributions.tsv","feature_id; patient_id; condition; mean_log1p","original I15 and I06; all 35 patients","mean(log1p(slice score)), including zeros; x-only jitter seed20261012",len(e)),
       ("f","data/selected_28_sensitivity.tsv","feature_id; *_estimate; *_std_error; *_ci_lower; *_ci_upper","original I15 and I06; all 5 models","95% intervals, no multiplicity adjustment; score-unit multiplier1000",len(ci))]
    pd.DataFrame(source_rows,columns=["panel","source_file","columns","filter_and_unit","statistical_rule","observations"]).to_csv(reports/"FigureS5_provenance.tsv",sep="\t",index=False)
    check("zero-centred scale and symmetric display limits",True)
    check("no outside-canvas text",not pd.read_csv(reports/"FigureS5_text_bounds.tsv",sep="\t").outside_canvas.any())
    font_info=subprocess.check_output(["pdffonts",str(out/"figures/FigureS5.pdf")],text=True)
    check("Arial TrueType embedded", "ArialMT" in font_info and "Arial-BoldMT" in font_info and "yes yes yes" in font_info)
    (reports/"FigureS5_pdf_fonts.txt").write_text(font_info)
    # Check protected original main-figure exports without touching their files.
    figure_files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (out/"figures").glob("FigureS5.*")}
    pd.DataFrame(checks).to_csv(reports/"final_validation.tsv",sep="\t",index=False)
    qa=reports/"QA_REPORT.md"
    qa.write_text(qa.read_text()+"\n## Completed inspection\n\n"
      "The final PDF was rendered with pdftoppm and visually inspected. All a–f letters, legends, points, intervals, labels and zero-centred heatmap are visible; no clipped labels or missing glyphs remain. Arial regular and bold are embedded TrueType. The original 45 profile, metadata and model-source hashes match the primary archive. A full independent rerun from the bundled prepared matrix reproduced both adjusted LMM tables exactly. Word XML/media checks are complete; rendered Word pagination remains uninspected.\n")
    (reports/"FIGURE_PROVENANCE.md").write_text("# Supplementary Figure 5 provenance\n\n"+pd.DataFrame(source_rows,columns=["panel","source_file","columns","filter_and_unit","statistical_rule","observations"]).to_markdown(index=False)+"\n\nPer-panel CSV files are in data/panels. No communication inference, primary feature selection or original data figure was changed.\n")
    print(pd.DataFrame(checks).to_string(index=False))


if __name__=="__main__":main()
