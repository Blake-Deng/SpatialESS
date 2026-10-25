"""Independent validation of the saved classification using patient-level PMFs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def count_pmf(patient_probabilities):
    # Independent Poisson-binomial dynamic program; no enumeration of choices.
    n, f = patient_probabilities.shape
    pmf = np.zeros((n + 1, f))
    pmf[0] = 1
    for p in patient_probabilities:
        old = pmf.copy()
        pmf[0] = old[0] * (1 - p)
        pmf[1:] = old[1:] * (1 - p) + old[:-1] * p
    return pmf


def main():
    d = pd.read_csv(ROOT / "results/all_50040_feature_classification.tsv.gz", sep="\t", low_memory=False)
    original = pd.read_csv(ROOT / "inputs/SD4_Lung_Support.tsv.gz", sep="\t", low_memory=False)
    score = pd.read_csv(ROOT / "inputs/slice_scores.tsv.gz", sep="\t", index_col=0).loc[:,d.feature_id]
    meta = pd.read_csv(ROOT / "inputs/sample_metadata.tsv", sep="\t").set_index("sample_id").loc[score.index]
    patients = meta.patient_id.drop_duplicates()
    probabilities, condition = [], []
    for p in patients:
        rows = meta.patient_id.eq(p).to_numpy()
        probabilities.append(score.to_numpy()[rows].astype(bool).mean(axis=0))
        condition.append(meta.loc[rows,"condition"].iloc[0])
    probabilities = np.stack(probabilities)
    disease = np.array(condition) == "Disease"
    pf = count_pmf(probabilities[disease])
    c = count_pmf(probabilities[~disease])
    errors = {}
    for pct in [40,50,60]:
        pmin = (pct * 26 + 99) // 100
        cmin = (pct * 9 + 99) // 100
        expected = {
            "pf_restricted": pf[pmin:].sum(axis=0) * c[0],
            "control_restricted": c[cmin:].sum(axis=0) * pf[0],
            "shared": pf[pmin:].sum(axis=0) * c[cmin:].sum(axis=0),
        }
        expected["other"] = 1 - sum(expected.values())
        errors[pct] = {}
        for label, probability in expected.items():
            error = float(abs(probability - d[f"one_slice_exact_p_{label}_{pct}pct"]).max())
            assert error < 2e-14, (pct,label,error)
            errors[pct][label] = error
    copied_columns = {}
    for column in original.columns:
        if pd.api.types.is_numeric_dtype(original[column]) and not pd.api.types.is_bool_dtype(original[column]):
            a, b = original[column].to_numpy(), d[column].to_numpy()
            assert np.array_equal(pd.isna(a), pd.isna(b))
            assert np.allclose(a,b,rtol=2e-14,atol=0,equal_nan=True), column
            copied_columns[column] = float(np.nanmax(np.abs(a-b)))
        else:
            pd.testing.assert_series_equal(original[column],d[column],check_names=False)
    assert original.q_value.lt(.05).equals(d.q_value.lt(.05))
    assert original.estimate.gt(0).equals(d.estimate.gt(0))
    stability = pd.read_csv(ROOT / "results/one_slice_stability_summary.tsv",sep="\t")
    counts = pd.read_csv(ROOT / "results/threshold_class_counts.tsv",sep="\t")
    for pct in [40,50,60]:
        assert counts.query("threshold_percent == @pct").features.sum() == 50040
        for row in stability.query("threshold_percent == @pct").itertuples():
            selected = d[f"support_class_{pct}pct"].eq(row.baseline_class)
            assert selected.sum() == row.features
            assert (selected & d[f"one_slice_exact_retention_{pct}pct"].eq(1)).sum() == row.always_same_label
    output = {
        "passed": True,
        "independent_method": "Poisson-binomial PMF from each patient's fraction of positive saved sections",
        "maximum_absolute_probability_errors_vs_768_choice_enumeration": errors,
        "original_numeric_column_roundtrip_max_abs_errors": copied_columns,
        "original_statistical_calls_and_effect_signs_unchanged": True,
        "all_50040_features_checked_at_all_3_thresholds": True,
        "note": "Tiny floating-point TSV serialization differences do not represent new estimation; original input is included byte-for-byte."
    }
    (ROOT / "reports/independent_validation.json").write_text(json.dumps(output,indent=2))
    print(json.dumps({k:v for k,v in output.items() if k != "original_numeric_column_roundtrip_max_abs_errors"},indent=2))


if __name__ == "__main__":
    main()
