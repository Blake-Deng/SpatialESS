"""Write supplementary descriptions from verified saved outputs."""
from pathlib import Path
import json
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def markdown_table(rows, columns):
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    lines += ["| " + " | ".join(str(x) for x in row) + " |" for row in rows]
    return "\n".join(lines)


def main():
    r = ROOT / "results"
    d = pd.read_csv(r / "all_50040_feature_classification.tsv.gz", sep="\t", low_memory=False)
    counts = pd.read_csv(r / "threshold_class_counts.tsv",sep="\t")
    stable = pd.read_csv(r / "one_slice_stability_summary.tsv",sep="\t")
    selection = pd.read_csv(r / "one_slice_class_count_summary.tsv",sep="\t")
    provenance = json.loads((ROOT / "reports/run_provenance.json").read_text())
    independent = json.loads((ROOT / "reports/independent_validation.json").read_text())
    assert independent["passed"]
    classes = ["PF-restricted","Control-restricted","Cross-condition shared","Other"]
    pivot = counts.pivot(index="category",columns="threshold_percent",values="features")
    threshold_table = markdown_table([[c,int(pivot.loc[c,40]),int(pivot.loc[c,50]),int(pivot.loc[c,60])]
                                      for c in classes],["观察支持类别","40%","50% 主分类","60%"])
    stability_table_rows = []
    for category in classes:
        row = stable[(stable.threshold_percent==50)&stable.baseline_class.eq(category)].iloc[0]
        ct = counts[(counts.threshold_percent==50)&counts.category.eq(category)].iloc[0]
        sm = selection[(selection.threshold_percent==50)&selection.category.eq(category)].iloc[0]
        ss = d[d.support_class_50pct.eq(category)]
        strong = ss.same_label_across_40_50_60pct.copy()
        for pct in [40,50,60]:
            strong &= ss[f"one_slice_exact_retention_{pct}pct"].eq(1)
        stability_table_rows.append([category,int(row.features),int(ct.primary_q_lt_005),
            int(row.always_same_label),int(strong.sum()),f"{int(sm.exact_count_min)}–{int(sm.exact_count_max)}"])
    stable_table = markdown_table(stability_table_rows,["50% 类别","主分类数量","原 primary q<0.05",
        "50% 下全部取片组合保留原标签","三阈值及全部取片组合均同类","各取片组合中该类总数范围"])
    changed40 = int(d.support_class_40pct.ne(d.support_class_50pct).sum())
    changed60 = int(d.support_class_60pct.ne(d.support_class_50pct).sum())
    spp = d[d.interaction_name.eq("SPP1_CD44") & d.sender_group_name.eq("Activated Fibrotic FBs")
            & d.receiver_group_name.eq("CD8+ T-cells")].iloc[0]
    chinese = f"""# 条件观察支持分类：结果与解释

本次是新增的探索性描述分析，使用 343 基因/37 LR 肺队列的既有 45×50,040 得分、35 位患者映射和最终 SD4。没有重新推断通信、重拟合 primary 模型或重做 BH。全部原特征均保留，包括无效/不合格拟合；拟合状态与分类并列。

## 定义

切片支持为保存的原始 score>0；主分类中的患者支持为该患者至少一张切片非零。分母固定为 PF 26 人、Control 9 人。50% 阈值分别需要至少 13 人、5 人；40% 为 11 人、4 人，60% 为 16 人、6 人。

- PF-restricted：PF 达支持阈值，Control 支持数为 0。
- Control-restricted：Control 达支持阈值，PF 支持数为 0。
- Cross-condition shared：两组均达到各自支持阈值；这是本分析明确采用的跨条件 consensus 定义。
- Other：其余，包含低支持、不均匀支持以及支持比例未满足阈值的特征，不能解释成“没有通信”。

“restricted”始终指本数据内观察支持的限制，不代表另一条件生物学缺失；“shared”不代表得分相等。50%、40%、60% 为现在明确的新规则，不能称历史预设。

## 主分类及阈值检查

{threshold_table}

相对 50%：40% 改变 {changed40:,} 条标签，60% 改变 {changed60:,} 条；{int(d.same_label_across_40_50_60pct.sum()):,}/50,040 条在三个阈值下标签一致。这项总体比例包含大量 Other，不能代替各目标类别的稳定性。

## 每人一张切片的敏感性

完成 1,000 次独立均匀抽取（每人每次一张，seed=20261024，NumPy PCG64）。每次保留全部 35 位患者，不重采样患者、不重新检验。另穷举全部 768 种组合（1 人有 3 张，8 人各 2 张，26 人各 1 张），等权得到精确的条件标签保留比例。它是对现有切片选择的稳定性，不是生物学真实性概率或患者总体置信区间。

{stable_table}

“保留原标签”的列按原主分类特征逐条检查；末列是在各取片组合中重新分类后的类别总数，可能包含原来属于其他类的特征，因此两列不可直接相减。原主分类 PF-restricted 的 2,061 条中，仅 1,081 条在 50% 下对所有取片选择保持该标签，提示不能把全部候选描述为稳健。

## 原有差异关联保持独立

仍为 39,570 个有效 primary fits、5,147 个原 primary BH 发现。50% 分类中 PF-restricted 的 1,311 条、Control-restricted 的 1 条、shared 的 930 条和 Other 的 2,905 条属于原发现。shared 中存在显著得分差异并不矛盾；q≥0.05 不证明等效。没有对分类内重做 BH。

## 现有例子的核实

Activated Fibrotic FBs → CD8+ T-cells 的 SPP1–CD44：PF 支持 29/35 切片、22/26 患者；Control 为 0/10、0/9。原 beta={spp.estimate:.8g}，q={spp.q_value:.8g}。三个阈值均为 PF-restricted；全部 768 种取片组合均保持，取片后 PF 支持为 20–22 人。可以作为具有切片选择稳健性的观察支持候选，但不证明 Control 生物学缺失或因果信号。

50% 下唯一 Control-restricted 为 AT2 → Monocytes/MDMs 的 CD8A–CEACAM5，Control 支持 5/9、PF 0/26。它只在 384/768 种组合保留标签，60% 时为 Other，不能突出为稳健的条件限制候选。

I15（FN1–integrin）及 I19 的 50% 标签是 Other，不能为了与正文叙事一致而重新标为 restricted。I06（COL1A1–integrin）在 40%/50% 为 shared、60% 为 Other，且 50% 的取片保留比例为 0.5。完整身份及各项数值见 existing_examples_classification.tsv。

## 检查、限制与交付

主脚本通过 {provenance['checks_passed']} 项检查，逐项重现 SD4 两组切片与患者支持数，并验证原列未被更改。独立验证通过患者支持数的 Poisson-binomial 动态规划核对全部特征、全部阈值、四类概率，与枚举最大差低于 2×10⁻¹⁴。文本保存/重读造成的极小浮点差记录在 independent_validation.json；原 SD4 输入按字节原样附带。

本分析没有区分零分数是否来自细胞群缺失、邻接缺失、表达覆盖或过滤；没有重新审核细胞注释，也不解决剩余的 run/region 混杂。不能把类别标签当作独立生物验证或新显著性检验。

入口：results/Condition_support_summary.xlsx；完整 50,040 条主表：results/all_50040_feature_classification.tsv.gz；异型细胞群方便子表：results/heterotypic_feature_classification.tsv.gz。保留脚本、参数、输入、SHA-256、随机及穷举取片清单和运行记录。主分析实测约 {provenance['elapsed_seconds']:.2f} 秒，不包括独立验证、文档和打包。
"""
    (ROOT / "reports/RESULTS_zh.md").write_text(chinese)
    english = f"""# Supplementary methods and results: observed condition support

## Methods

We performed a post hoc descriptive classification of the 50,040 sender-group–receiver-group–ligand–receptor features in the frozen 343-gene/37-record GSE250346 cohort analysis. We reused the saved nonnegative section scores and metadata for 45 sections from 35 patients (26 PF and 9 control), without rerunning communication inference, fitting models or recalculating Benjamini–Hochberg corrections. A section was counted as supporting a feature when its saved score was greater than zero. In the primary descriptive classification, a patient was counted as supported if at least one available section supported the feature. All patients in each condition remained in the denominator.

At the 50% patient-support threshold, a feature was classified as PF-restricted observed communication when it was supported in at least 13/26 PF patients and no control patients; the reciprocal criterion defined control-restricted observed communication (at least 5/9 controls and no PF patients). Cross-condition shared communication required support in at least 13/26 PF and 5/9 control patients; all other features were classified as Other. This defines consensus as shared support across conditions, rather than equality of scores or recurrence within a single condition. The four labels were mutually exclusive and covered the complete feature universe. Threshold checks used 40% (11 PF and 4 controls) and 60% (16 PF and 6 controls).

To assess dependence on the number of available sections, we independently selected one existing section uniformly per patient in each of 1,000 repetitions (NumPy Generator/PCG64, seed 20261024). All 35 patients remained in every selection. We additionally enumerated all 768 possible selections, giving equal weight to every combination; one patient had three sections, eight had two, and 26 had one. Per-feature category probabilities and the fraction of selections retaining the primary support label were calculated at each threshold. An independent Poisson-binomial calculation based on each patient's fraction of positive sections verified the exact enumeration. Selection-distribution percentiles describe variation among the existing section choices and are not inferential confidence intervals.

Original primary estimates, standard errors, test statistics, raw P values, BH q values, confidence intervals, fit status and singular-fit flags were retained. Differential-association labels use the original status-ok primary q < 0.05 and effect direction; no category-specific BH correction was applied. Classification thresholds were introduced in this supplementary analysis and were not represented as historical prespecified rules.

## Results

At the 50% threshold, 2,061 features were PF-restricted, one was control-restricted, 7,567 were cross-condition shared and 40,411 were Other. Of these, 1,311, one, 930 and 2,905, respectively, belonged to the original 5,147 primary BH discoveries. At 40%, the respective counts were 2,994, nine, 9,353 and 37,684; at 60%, they were 1,170, zero, 5,906 and 42,964. Compared with 50%, 2,727 labels changed at 40% and 2,553 at 60%.

Across every possible one-section-per-patient selection, 1,081 of the 2,061 primary PF-restricted features, zero of the one control-restricted feature, 6,365 of the 7,567 shared features and 40,176 of the 40,411 Other features retained their 50% labels. Reclassified per-selection counts ranged from 1,502 to 1,770 PF-restricted, zero to three control-restricted, 6,719 to 7,226 shared and 41,238 to 41,645 Other features. These count ranges include features transitioning from other baseline categories.

SPP1–CD44 from activated fibrotic fibroblasts to CD8+ T cells showed nonzero scores in 29/35 PF sections from 22/26 PF patients and in 0/10 control sections from 0/9 control patients. It remained PF-restricted at all three support thresholds and in all 768 section selections; 20–22 PF patients supported it after selection. Its original condition coefficient ({spp.estimate:.8g}) and primary BH q ({spp.q_value:.8g}) were retained. This supports a reproducible observed-support pattern within the available sections, without establishing biological absence in controls or causal signaling. The sole 50% control-restricted feature retained its label in only half of the exact selections and was Other at 60%, indicating sensitivity to the descriptive rules.

## Interpretation

Observed restriction and differential association measure distinct properties. A shared feature can have a differential score, and a nonsignificant association is not evidence of equivalence. Zero scores can reflect cell-composition, spatial-adjacency, assay-coverage or filtering effects. The section-selection analysis does not resolve these mechanisms, validate annotations, remove run/region confounding or provide independent biological validation.
"""
    (ROOT / "reports/METHODS_and_RESULTS.md").write_text(english)
    notes = [
        ("analysis", "New post hoc descriptive analysis; no inference, refitting or BH recalculation"),
        ("unit", "Complete sender -> receiver x LR feature_id, not LR alone"),
        ("patient_support", "At least one saved section score > 0; denominator all 26 PF / 9 Control patients"),
        ("50% threshold", "PF >=13/26; Control >=5/9"),
        ("40% threshold", "PF >=11/26; Control >=4/9"),
        ("60% threshold", "PF >=16/26; Control >=6/9"),
        ("PF-restricted", "PF reaches threshold, zero supported controls; observed support only"),
        ("Control-restricted", "Controls reach threshold, zero supported PF patients; observed support only"),
        ("Cross-condition shared", "Both conditions reach threshold; does not imply equal scores"),
        ("Other", "Remainder; not an absence-of-communication category"),
        ("primary_q_lt_005", "Original status-ok primary global BH q<0.05; no new BH"),
        ("sampling", "1000 one-section-per-patient random repetitions, seed 20261024, plus all 768 exact choices"),
        ("exact_retention", "Fraction of choices retaining the full-section baseline label at the same threshold"),
        ("exact_probability", "Probability of a descriptive label over equally weighted existing section choices, not a biological probability"),
        ("exact_count_range", "Reclassified count in each choice; may contain transitions from another baseline class"),
        ("descriptive percentiles", "Variation over saved section selections, not inferential confidence intervals"),
        ("no_detected_differential_association", "Original valid fit without q<0.05; not evidence of equivalence"),
        ("limitations", "No decomposition of zeros or new annotation/biological validation; no resolution of run/region confounding"),
    ]
    with pd.ExcelWriter(r / "Condition_support_summary.xlsx",mode="a",if_sheet_exists="replace",engine="openpyxl") as writer:
        pd.DataFrame(notes,columns=["item","definition"]).to_excel(writer,sheet_name="Read_me",index=False)
        ws=writer.sheets["Read_me"];ws.freeze_panes="A2"
        ws.column_dimensions["A"].width=38;ws.column_dimensions["B"].width=110
        # Put the definitions before the numeric sheets for standalone readers.
        writer.book.move_sheet(ws,offset=-writer.book.index(ws))
    (ROOT / "reports/report_metrics.json").write_text(json.dumps(dict(
        threshold_40_changed_vs_50=changed40,threshold_60_changed_vs_50=changed60,
        stability_summary_columns=["category","main_count","original_primary_hits","50pct_always_same",
                                   "all_thresholds_all_choices_same","one_slice_class_count_range"],
        stability_summary=stability_table_rows),indent=2))
    print("Wrote Chinese and English reports and workbook definitions.")


if __name__ == "__main__":
    main()
