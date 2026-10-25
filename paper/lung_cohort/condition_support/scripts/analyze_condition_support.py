"""Classify existing saved scores only. No inference, model fit, or BH computation."""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import math
import platform
import time

import numpy as np
import pandas as pd


CATEGORIES = np.array([
    "PF-restricted", "Control-restricted", "Cross-condition shared", "Other"
])


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_tsv(frame, path):
    if path.suffix == ".gz":
        frame.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    else:
        frame.to_csv(path, sep="\t", index=False)


def classify(pf, control, percentage, npf, nc):
    # Integer arithmetic implements >= threshold without a rounding ambiguity.
    pf_min = (percentage * npf + 99) // 100
    c_min = (percentage * nc + 99) // 100
    labels = np.full(pf.shape, 3, dtype=np.uint8)
    labels[(pf >= pf_min) & (control == 0)] = 0
    labels[(control >= c_min) & (pf == 0)] = 1
    labels[(pf >= pf_min) & (control >= c_min)] = 2
    return labels


def choice_support(presence, choices, patient_rows, disease):
    shape = (len(choices), presence.shape[1])
    pf = np.zeros(shape, dtype=np.uint8)
    control = np.zeros(shape, dtype=np.uint8)
    for i, rows in enumerate(patient_rows):
        target = pf if disease[i] else control
        if len(rows) == 1:
            target += presence[rows[0]]
        else:
            target += presence[choices[:, i]]
    return pf, control


def selection_ledger(choices, patients, metadata, kind):
    repeats, npatients = choices.shape
    return pd.DataFrame({
        "selection_type": kind,
        "selection_id": np.repeat(np.arange(1, repeats + 1), npatients),
        "patient_id": np.tile(patients, repeats),
        "sample_id": metadata.index.to_numpy()[choices.reshape(-1)],
        "condition": metadata.condition.to_numpy()[choices.reshape(-1)],
    })


def main(root):
    started = time.perf_counter()
    config = json.loads((root / "config/design.json").read_text())
    results = root / "results"
    reports = root / "reports"
    source = root / "inputs"
    for folder in (results, reports):
        folder.mkdir(exist_ok=True)
    checks = []

    def check(name, condition, detail=""):
        checks.append({"check": name, "passed": bool(condition), "detail": detail})
        if not condition:
            write_tsv(pd.DataFrame(checks), reports / "validation.tsv")
            raise AssertionError(name + ": " + detail)

    print("Loading saved inputs", flush=True)
    sd4 = pd.read_csv(source / "SD4_Lung_Support.tsv.gz", sep="\t", low_memory=False)
    scores = pd.read_csv(source / "slice_scores.tsv.gz", sep="\t", index_col=0)
    metadata = pd.read_csv(source / "sample_metadata.tsv", sep="\t").set_index("sample_id")
    check("unique full feature keys", not sd4.feature_id.duplicated().any() and scores.columns.is_unique)
    check("fixed feature universe", len(sd4) == 50040 and set(sd4.feature_id) == set(scores.columns))
    check("unique slice identifiers", scores.index.is_unique and metadata.index.is_unique)
    check("identical slice identifiers", set(scores.index) == set(metadata.index))
    scores = scores.loc[:, sd4.feature_id]
    metadata = metadata.loc[scores.index]
    values = scores.to_numpy()
    check("45 x 50040 scores", values.shape == (45, 50040))
    check("complete finite nonnegative saved scores", np.isfinite(values).all() and (values >= 0).all())
    presence = values > 0
    check("saved nonzero record count", int(presence.sum()) == 692338)
    patients = metadata.patient_id.drop_duplicates().to_numpy()
    patient_rows = [np.flatnonzero(metadata.patient_id.to_numpy() == p) for p in patients]
    check("constant condition within patient", all(metadata.iloc[r].condition.nunique() == 1 for r in patient_rows))
    patient_condition = np.array([metadata.iloc[r[0]].condition for r in patient_rows])
    check("only two original condition labels", set(patient_condition) == {"Control", "Disease"})
    disease = patient_condition == "Disease"
    npf, nc = int(disease.sum()), int((~disease).sum())
    check("35 patients: 26 PF and 9 Control", len(patients) == 35 and npf == 26 and nc == 9)
    patient_any = np.stack([presence[r].any(axis=0) for r in patient_rows])
    patient_all = np.stack([presence[r].all(axis=0) for r in patient_rows])
    pf_base, c_base = patient_any[disease].sum(axis=0), patient_any[~disease].sum(axis=0)
    for name, actual in [
        ("pf_nonzero_patients", pf_base), ("control_nonzero_patients", c_base),
        ("pf_nonzero_slices", presence[metadata.condition.eq("Disease")].sum(axis=0)),
        ("control_nonzero_slices", presence[metadata.condition.eq("Control")].sum(axis=0)),
    ]:
        check("SD4 support reproduced: " + name, np.array_equal(actual, sd4[name].to_numpy()))
    check("SD4 denominators", sd4.pf_total_patients.eq(npf).all() and sd4.control_total_patients.eq(nc).all()
          and sd4.pf_total_slices.eq(35).all() and sd4.control_total_slices.eq(10).all())

    output = sd4.copy(deep=True)
    output["pf_patient_support_fraction"] = pf_base / npf
    output["control_patient_support_fraction"] = c_base / nc
    output["pf_patients_nonzero_in_all_available_slices"] = patient_all[disease].sum(axis=0)
    output["control_patients_nonzero_in_all_available_slices"] = patient_all[~disease].sum(axis=0)
    valid = sd4.status.eq("ok").to_numpy()
    hit = valid & sd4.q_value.lt(.05).to_numpy()
    heterotypic = sd4.sender_group_name.ne(sd4.receiver_group_name).to_numpy()
    output["primary_association_label"] = np.select([
        ~valid, hit & sd4.estimate.gt(0).to_numpy(), hit & sd4.estimate.lt(0).to_numpy()
    ], ["invalid_or_ineligible_fit", "PF-associated increase", "PF-associated decrease"],
        default="no_detected_differential_association")
    check("original primary counts", int(valid.sum()) == 39570 and int(hit.sum()) == 5147)
    baseline = {}
    count_rows, transition_rows = [], []
    for pct in config["threshold_percentages"]:
        label = classify(pf_base, c_base, pct, npf, nc)
        baseline[pct] = label
        output[f"support_class_{pct}pct"] = CATEGORIES[label]
        for k, category in enumerate(CATEGORIES):
            selected = label == k
            count_rows.append(dict(threshold_percent=pct, category=category,
                pf_required_patients=math.ceil(pct*npf/100), control_required_patients=math.ceil(pct*nc/100),
                features=int(selected.sum()), fraction_of_50040=float(selected.mean()),
                valid_primary_fits=int((selected & valid).sum()), primary_q_lt_005=int((selected & hit).sum()),
                primary_q_lt_005_beta_positive=int((selected & hit & sd4.estimate.gt(0).to_numpy()).sum()),
                primary_q_lt_005_beta_negative=int((selected & hit & sd4.estimate.lt(0).to_numpy()).sum()),
                heterotypic_features=int((selected & heterotypic).sum()),
                heterotypic_primary_q_lt_005=int((selected & heterotypic & hit).sum())))
    main_pct = config["main_threshold_percentage"]
    for pct in config["threshold_percentages"]:
        output[f"same_as_50pct_at_{pct}pct"] = baseline[pct] == baseline[main_pct]
        for i, origin in enumerate(CATEGORIES):
            for j, destination in enumerate(CATEGORIES):
                selected = (baseline[main_pct] == i) & (baseline[pct] == j)
                transition_rows.append(dict(from_threshold=main_pct, to_threshold=pct,
                    from_class=origin, to_class=destination, features=int(selected.sum()),
                    primary_q_lt_005=int((selected & hit).sum())))
    output["same_label_across_40_50_60pct"] = (baseline[40] == baseline[50]) & (baseline[60] == baseline[50])
    count_table = pd.DataFrame(count_rows)
    transitions = pd.DataFrame(transition_rows)
    write_tsv(count_table, results / "threshold_class_counts.tsv")
    write_tsv(transitions, results / "threshold_transitions.tsv")

    total_choices = math.prod(len(r) for r in patient_rows)
    check("finite exact one-slice choice space", total_choices == 768)
    exact = np.array(list(itertools.product(*patient_rows)), dtype=np.int16)
    rng = np.random.Generator(np.random.PCG64(config["random_resampling"]["seed"]))
    repetitions = config["random_resampling"]["repetitions"]
    monte_carlo = np.column_stack([rng.choice(rows, size=repetitions) for rows in patient_rows])
    selection_order = {tuple(row): i for i, row in enumerate(exact)}
    mc_exact_index = np.array([selection_order[tuple(row)] for row in monte_carlo])
    check("each resample has exactly one section per patient", all(
        np.isin(monte_carlo[:, i], rows).all() for i, rows in enumerate(patient_rows)))
    write_tsv(selection_ledger(exact, patients, metadata, "exact"), results / "exact_slice_selections.tsv.gz")
    write_tsv(selection_ledger(monte_carlo, patients, metadata, "monte_carlo"), results / "random_slice_selections.tsv.gz")
    patient_table = pd.DataFrame(dict(patient_id=patients, condition=patient_condition,
                                     available_sections=[len(r) for r in patient_rows]))
    patient_table["sample_ids"] = [";".join(metadata.index[r]) for r in patient_rows]
    write_tsv(patient_table, results / "patient_sampling_units.tsv")
    print(f"Enumerating {total_choices} exact choices; {repetitions} random repetitions", flush=True)
    pf_choices, c_choices = choice_support(presence, exact, patient_rows, disease)
    for condition, support in [("pf", pf_choices), ("control", c_choices)]:
        output[f"one_slice_{condition}_nonzero_patients_min"] = support.min(axis=0)
        output[f"one_slice_{condition}_nonzero_patients_max"] = support.max(axis=0)
        output[f"one_slice_{condition}_nonzero_patients_exact_mean"] = support.mean(axis=0)
    summary_rows, stability_rows, per_choice_rows, sampling_transition_rows = [], [], [], []
    exact_main = None
    for pct in config["threshold_percentages"]:
        labels = classify(pf_choices, c_choices, pct, npf, nc)
        frequencies = np.stack([(labels == i).sum(axis=0) for i in range(4)])
        mc_frequencies = np.stack([(labels[mc_exact_index] == i).sum(axis=0) for i in range(4)])
        check(f"exact probabilities sum to one at {pct}%", np.all(frequencies.sum(axis=0) == total_choices))
        check(f"random probabilities sum to one at {pct}%", np.all(mc_frequencies.sum(axis=0) == repetitions))
        retention_counts = frequencies[baseline[pct], np.arange(len(sd4))]
        mode = frequencies.argmax(axis=0)
        ties = (frequencies == frequencies.max(axis=0)).sum(axis=0) > 1
        output[f"one_slice_exact_retention_{pct}pct"] = retention_counts / total_choices
        output[f"one_slice_exact_modal_class_{pct}pct"] = CATEGORIES[mode]
        output[f"one_slice_exact_modal_tie_{pct}pct"] = ties
        for i, category in enumerate(CATEGORIES):
            token = ["pf_restricted", "control_restricted", "shared", "other"][i]
            output[f"one_slice_exact_p_{token}_{pct}pct"] = frequencies[i] / total_choices
            if pct == main_pct:
                output[f"one_slice_mc_p_{token}_{pct}pct"] = mc_frequencies[i] / repetitions
            exact_counts = (labels == i).sum(axis=1)
            mc_counts = exact_counts[mc_exact_index]
            summary_rows.append(dict(threshold_percent=pct, category=category,
                all_slices_any_patient_count=int((baseline[pct] == i).sum()),
                exact_combinations=total_choices, exact_count_mean=float(exact_counts.mean()),
                exact_count_min=int(exact_counts.min()), exact_count_max=int(exact_counts.max()),
                exact_count_descriptive_p025=float(np.quantile(exact_counts,.025)),
                exact_count_descriptive_p975=float(np.quantile(exact_counts,.975)),
                mc_repetitions=repetitions, mc_count_mean=float(mc_counts.mean()),
                mc_count_min=int(mc_counts.min()), mc_count_max=int(mc_counts.max())))
            for kind, counts in [("exact", exact_counts), ("monte_carlo", mc_counts)]:
                per_choice_rows.extend(dict(selection_type=kind, selection_id=j+1,
                    threshold_percent=pct, category=category, features=int(n)) for j,n in enumerate(counts))
            selected = baseline[pct] == i
            stability_rows.append(dict(threshold_percent=pct, baseline_class=category,
                features=int(selected.sum()), always_same_label=int((selected & (retention_counts == total_choices)).sum()),
                retention_at_least_095=int((selected & (retention_counts / total_choices >= .95)).sum()),
                ever_changed_label=int((selected & (retention_counts < total_choices)).sum()),
                mean_exact_retention=float((retention_counts[selected]/total_choices).mean()) if selected.any() else None,
                median_exact_retention=float(np.median(retention_counts[selected]/total_choices)) if selected.any() else None,
                max_mc_probability_error=float(np.abs(mc_frequencies[:,selected]/repetitions - frequencies[:,selected]/total_choices).max()) if selected.any() else None))
            for j, target in enumerate(CATEGORIES):
                sampling_transition_rows.append(dict(threshold_percent=pct, from_class=category, to_class=target,
                    exact_expected_features_per_choice=float(frequencies[j,selected].sum()/total_choices),
                    features_ever_in_target=int((selected & (frequencies[j]>0)).sum())))
        if pct == main_pct:
            exact_main = labels
            output["one_slice_mc_retention_50pct"] = mc_frequencies[baseline[pct],np.arange(len(sd4))]/repetitions
            output["one_slice_exact_retained_in_all_768_choices_50pct"] = retention_counts == total_choices
    # Independent direct check of a subset of the MC selections against exact lookup.
    direct_pf, direct_control = choice_support(presence, monte_carlo[:8], patient_rows, disease)
    check("random labels independently match exact selection lookup", np.array_equal(
        classify(direct_pf,direct_control,50,npf,nc),exact_main[mc_exact_index[:8]]))
    for column in sd4.columns:
        check("original column preserved: " + column, output[column].equals(sd4[column]))
    check("four classes exhaust main universe", sum(row["features"] for row in count_rows if row["threshold_percent"] == 50) == 50040)
    write_tsv(output, results / "all_50040_feature_classification.tsv.gz")
    write_tsv(output.loc[heterotypic], results / "heterotypic_feature_classification.tsv.gz")
    write_tsv(pd.DataFrame(summary_rows), results / "one_slice_class_count_summary.tsv")
    write_tsv(pd.DataFrame(stability_rows), results / "one_slice_stability_summary.tsv")
    write_tsv(pd.DataFrame(per_choice_rows), results / "one_slice_counts_by_selection.tsv.gz")
    write_tsv(pd.DataFrame(sampling_transition_rows), results / "one_slice_expected_transitions.tsv")
    # Preserve the complete universe. Selected feature rows are provided only as an index.
    representatives = output[output.display_id.notna() | (
        output.interaction_name.eq("SPP1_CD44") & output.sender_group_name.eq("Activated Fibrotic FBs")
        & output.receiver_group_name.eq("CD8+ T-cells"))]
    write_tsv(representatives, results / "existing_examples_classification.tsv")
    excel_tables = {"Threshold_counts":count_table,"Threshold_transitions":transitions,
        "One_slice_counts":pd.DataFrame(summary_rows),"One_slice_stability":pd.DataFrame(stability_rows),
        "Sampling_transitions":pd.DataFrame(sampling_transition_rows),"Patient_units":patient_table,
        "Existing_examples":representatives}
    with pd.ExcelWriter(results / "Condition_support_summary.xlsx", engine="openpyxl") as writer:
        for sheet, frame in excel_tables.items():
            frame.to_excel(writer,sheet_name=sheet,index=False)
            ws=writer.sheets[sheet]
            ws.freeze_panes="A2"
            ws.auto_filter.ref=ws.dimensions
            for col in ws.columns:
                letter=col[0].column_letter
                ws.column_dimensions[letter].width=min(48,max(12,len(str(col[0].value))+2))
    write_tsv(pd.DataFrame(checks),reports / "validation.tsv")
    provenance = dict(analysis=config["analysis"], timestamp=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,
        platform=platform.platform(), elapsed_seconds=time.perf_counter()-started,
        n_features=len(sd4),n_sections=len(metadata),n_patients=len(patients),
        npf=npf,ncontrol=nc,exact_choices=total_choices,random_repetitions=repetitions,
        random_seed=config["random_resampling"]["seed"],
        inputs={p.name:dict(sha256=sha256(p),bytes=p.stat().st_size,
                             original_source=config["original_sources"][p.name]) for p in source.iterdir() if p.is_file()},
        script_sha256=sha256(Path(__file__)),
        thresholds_agree_all_three=int(output.same_label_across_40_50_60pct.sum()),
        main_exact_always_same=int(output.one_slice_exact_retained_in_all_768_choices_50pct.sum()),
        checks_passed=len(checks))
    (reports / "run_provenance.json").write_text(json.dumps(provenance,indent=2))
    print(count_table.to_string(index=False),flush=True)
    print(pd.DataFrame(stability_rows).query("threshold_percent == 50").to_string(index=False),flush=True)
    print(f"Completed: {len(checks)} checks; {provenance['elapsed_seconds']:.2f} seconds",flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    main(parser.parse_args().root.resolve())
