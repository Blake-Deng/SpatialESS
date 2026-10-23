"""Independent-patient, heteroskedasticity and influence sensitivity analyses.

No communication inference or primary model is changed. Patient outcomes are
mean(log1p(slice_score)), including zeros. BH is applied to each full valid set.
"""
from pathlib import Path
import argparse
import json
import platform
import warnings

import numpy as np
import pandas as pd
import scipy
from scipy import stats

SEED = 20261012


def bh(p):
    q = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    q[valid] = stats.false_discovery_control(p[valid], method="bh")
    return q


def hc3(x, y):
    n, k = x.shape
    if np.linalg.matrix_rank(x) != k:
        raise ValueError("Patient sensitivity design is rank-deficient")
    inverse = np.linalg.inv(x.T @ x)
    linear_weights = inverse @ x.T
    coefficients = linear_weights @ y
    residual = y - x @ coefficients
    leverage = np.sum(x * (x @ inverse), axis=1)
    if np.max(leverage) >= 1 - 1e-10:
        raise ValueError("HC3 has unit leverage")
    variance = (linear_weights[1] ** 2 / (1 - leverage) ** 2) @ (residual ** 2)
    se = np.sqrt(np.maximum(variance, 0))
    eligible = (np.sum(y > 0, axis=0) >= 2) & (np.std(y, axis=0, ddof=1) > 1e-12)
    valid = eligible & np.isfinite(se) & (se > 0)
    p = np.full(y.shape[1], np.nan)
    p[valid] = 2 * stats.t.sf(np.abs(coefficients[1, valid] / se[valid]), n - k)
    critical = stats.t.ppf(.975, n - k)
    return pd.DataFrame(dict(estimate=coefficients[1], std_error=se,
                            ci_lower=coefficients[1] - critical * se,
                            ci_upper=coefficients[1] + critical * se,
                            p_value=p, q_value=bh(p), residual_df=n-k,
                            nonzero_patients=np.sum(y > 0, axis=0), total_patients=n,
                            status=np.where(valid, "ok", "insufficient_support_or_variation")))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--figure-source", required=True, type=Path)
    parser.add_argument("--replicates", type=int, default=1000)
    args = parser.parse_args()
    out = args.out
    meta = pd.read_csv(out / "data/sample_metadata.tsv", sep="\t")
    score = pd.read_csv(out / "data/slice_scores.tsv.gz", sep="\t", index_col=0)
    meta = meta.set_index("sample_id").loc[score.index]
    assert not score.isna().any().any() and (score >= 0).all().all()
    transformed = np.log1p(score.to_numpy())
    primary = pd.read_csv(out / "data/primary_lmm.tsv.gz", sep="\t").set_index("feature_id").loc[score.columns]
    patients = meta.patient_id.drop_duplicates().to_list()
    y = np.stack([transformed[meta.patient_id.to_numpy() == p].mean(axis=0) for p in patients])
    rows = []
    for p in patients:
        part = meta[meta.patient_id == p]
        assert part.condition.nunique() == 1
        rows.append(dict(patient_id=p, condition=part.condition.iloc[0], slices=len(part),
                         fraction_run2=(part.run == "Run2").mean(),
                         fraction_run3=(part.run == "Run3").mean()))
    pm = pd.DataFrame(rows)
    disease = (pm.condition == "Disease").to_numpy()
    assert disease.sum() == 26 and (~disease).sum() == 9
    x = np.column_stack([np.ones(len(pm)), disease, pm.fraction_run2, pm.fraction_run3])
    pm.to_csv(out / "data/patient_metadata.tsv", sep="\t", index=False)
    pd.DataFrame(y, index=patients, columns=score.columns).to_csv(out / "data/patient_mean_log1p.tsv.gz", sep="\t", index_label="patient_id")
    print("35 patients: 9 Control, 26 PF; patient response=mean(log1p(slice score)); HC3 df=31", flush=True)
    patient = hc3(x, y)
    patient.insert(0, "feature_id", score.columns)
    patient.to_csv(out / "data/patient_run_hc3.tsv.gz", sep="\t", index=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        welch = stats.ttest_ind(y[disease], y[~disease], equal_var=False, axis=0)
    difference = y[disease].mean(axis=0) - y[~disease].mean(axis=0)
    se = np.sqrt(y[disease].var(axis=0, ddof=1)/26 + y[~disease].var(axis=0, ddof=1)/9)
    eligible = (np.sum(y > 0, axis=0) >= 2) & (np.std(y, axis=0, ddof=1) > 1e-12)
    p = np.where(eligible & np.isfinite(welch.pvalue), welch.pvalue, np.nan)
    pd.DataFrame(dict(feature_id=score.columns, estimate=difference, std_error=se,
                      ci_lower=difference - stats.t.ppf(.975, welch.df)*se,
                      ci_upper=difference + stats.t.ppf(.975, welch.df)*se,
                      p_value=p, q_value=bh(p), residual_df=welch.df,
                      nonzero_patients=np.sum(y > 0, axis=0), total_patients=35,
                      status=np.where(np.isfinite(p), "ok", "insufficient_support_or_variation")))\
        .to_csv(out / "data/patient_welch.tsv.gz", sep="\t", index=False)
    # Only descriptive effect stability is inferred from these repeated subsets.
    rng = np.random.default_rng(SEED)
    weights = np.zeros((args.replicates, len(pm)))
    controls = np.flatnonzero(~disease); cases = np.flatnonzero(disease)
    selections = []
    for b in range(args.replicates):
        selected_cases = rng.choice(cases, size=9, replace=False)
        selected = np.concatenate([controls, selected_cases])
        if np.linalg.matrix_rank(x[selected]) != x.shape[1]:
            raise ValueError("Balanced sampling generated rank-deficient design")
        weights[b, selected] = np.linalg.pinv(x[selected])[1]
        selections.extend(dict(replicate=b+1, patient_id=patients[i], condition=pm.condition.iloc[i]) for i in selected)
    pd.DataFrame(selections).to_csv(out / "data/balanced_patient_selections.tsv.gz", sep="\t", index=False)
    balanced_parts = []
    full_estimate = patient.estimate.to_numpy()
    for start in range(0, y.shape[1], 2000):
        end = min(y.shape[1], start+2000)
        effects = weights @ y[:, start:end]
        quantiles = np.quantile(effects, [.025, .5, .975], axis=0)
        balanced_parts.append(pd.DataFrame(dict(feature_id=score.columns[start:end],
          descriptive_p025=quantiles[0], median=quantiles[1], descriptive_p975=quantiles[2],
          minimum=effects.min(axis=0), maximum=effects.max(axis=0),
          direction_agreement=np.mean(np.sign(effects) == np.sign(full_estimate[start:end]), axis=0),
          replicates=args.replicates)))
    pd.concat(balanced_parts).to_csv(out / "data/balanced_run_adjusted_effects.tsv.gz", sep="\t", index=False)
    print("Balanced 9-vs-9 run-adjusted effect stability completed", flush=True)
    estimates = []; qvalues = []; influence_rows = []
    original_hits = patient.q_value.to_numpy() < .05
    for i,pid in enumerate(patients):
        retained = np.arange(len(pm)) != i
        fit = hc3(x[retained], y[retained])
        estimates.append(fit.estimate.to_numpy()); qvalues.append(fit.q_value.to_numpy())
        valid = (patient.status == "ok").to_numpy() & (fit.status == "ok").to_numpy()
        hits = fit.q_value.to_numpy() < .05
        influence_rows.append(dict(omitted_patient=pid, condition=pm.condition.iloc[i], valid_features=valid.sum(),
          bh_hits=hits.sum(), original_hits_retained=(hits & original_hits).sum(),
          original_hit_fraction_retained=(hits & original_hits).sum()/max(1,original_hits.sum()),
          estimate_spearman=stats.spearmanr(full_estimate[valid],fit.estimate.to_numpy()[valid]).statistic,
          direction_agreement=np.mean(np.sign(full_estimate[valid]) == np.sign(fit.estimate.to_numpy()[valid]))))
    estimates=np.stack(estimates); qvalues=np.stack(qvalues)
    has_valid=np.any(np.isfinite(qvalues),axis=0)
    maximum_q=np.full(y.shape[1], np.nan)
    maximum_q[has_valid]=np.nanmax(qvalues[:,has_valid],axis=0)
    pd.DataFrame(dict(feature_id=score.columns, estimate_min=estimates.min(axis=0),
      estimate_max=estimates.max(axis=0), direction_agreement=np.mean(np.sign(estimates)==np.sign(full_estimate),axis=0),
      valid_omissions=np.sum(np.isfinite(qvalues),axis=0), bh_hit_omissions=np.sum(qvalues < .05,axis=0),
      maximum_q=maximum_q)).to_csv(out / "data/leave_one_patient_feature_summary.tsv.gz",sep="\t",index=False)
    pd.DataFrame(influence_rows).to_csv(out / "data/leave_one_patient_summary.tsv",sep="\t",index=False)
    mapping=pd.read_csv(args.figure_source / "feature_mapping.csv")
    mapping.to_csv(out / "data/feature_mapping.csv",index=False)
    selected=np.array([score.columns.get_loc(f) for f in mapping.feature_id])
    pd.DataFrame(estimates[:, selected],index=patients,columns=mapping.display_id).to_csv(
        out / "data/selected_leave_one_patient_estimates.tsv",sep="\t",index_label="omitted_patient")
    pd.DataFrame(qvalues[:, selected],index=patients,columns=mapping.display_id).to_csv(
        out / "data/selected_leave_one_patient_qvalues.tsv",sep="\t",index_label="omitted_patient")
    # Compare the vectorized sandwich calculation with an independent library.
    import statsmodels.api as sm
    verification=[]
    for j in selected:
        fit=sm.OLS(y[:,j],x).fit(cov_type="HC3",use_t=True)
        verification.append(dict(feature_id=score.columns[j],
          beta_difference=abs(fit.params[1]-patient.estimate.iloc[j]),
          se_difference=abs(fit.bse[1]-patient.std_error.iloc[j]),
          p_difference=abs(fit.pvalues[1]-patient.p_value.iloc[j])))
    verification=pd.DataFrame(verification)
    assert (verification[["beta_difference","se_difference","p_difference"]].max() < 1e-10).all()
    verification.to_csv(out / "reports/hc3_independent_validation.tsv",sep="\t",index=False)
    support=[]; distributions=[]
    for _,r in mapping.iterrows():
        j=score.columns.get_loc(r.feature_id)
        for condition in ["Control","Disease"]:
            sl=meta.condition.to_numpy() == condition
            pat=pm.condition.to_numpy() == condition
            support.append(dict(feature_id=r.feature_id, display_id=r.display_id, condition=condition,
              slices=int(sl.sum()), nonzero_slices=int((score.iloc[sl,j]>0).sum()),
              nonzero_slice_fraction=float((score.iloc[sl,j]>0).mean()),
              patients=int(pat.sum()), nonzero_patients=int((y[pat,j]>0).sum()),
              patient_mean_log1p=float(y[pat,j].mean()), patient_median_log1p=float(np.median(y[pat,j])),
              slice_median_raw_score=float(np.median(score.iloc[sl,j])), slice_mean_raw_score=float(score.iloc[sl,j].mean())))
        for i,pid in enumerate(patients):
            distributions.append(dict(feature_id=r.feature_id,display_id=r.display_id, patient_id=pid,
                                       condition=pm.condition.iloc[i], mean_log1p=y[i,j]))
    pd.DataFrame(support).to_csv(out / "data/selected_support_and_distributions.tsv",sep="\t",index=False)
    pd.DataFrame(distributions).to_csv(out / "data/selected_patient_distributions.tsv",sep="\t",index=False)
    versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,
                  scipy=scipy.__version__,statsmodels=sm.__version__ if hasattr(sm,"__version__") else __import__("statsmodels").__version__,seed=SEED,replicates=args.replicates)
    (out / "logs/python_environment.json").write_text(json.dumps(versions,indent=2)+"\n")
    (out / "logs/PATIENT_SENSITIVITY_COMPLETE").write_text("complete\n")
    print("Patient HC3 hits:",int(original_hits.sum()),"Welch hits:",int((bh(p)<.05).sum()),flush=True)
    print("Patient sensitivity COMPLETE",flush=True)


if __name__ == "__main__":
    main()
