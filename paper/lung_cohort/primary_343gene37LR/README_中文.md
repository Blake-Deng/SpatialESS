# 肺纤维化资料：正文 343 基因 / 37 LR

从 `data/complete_lmm_all_50040_features.tsv.gz` 开始筛选：完整 50,040 个特征，保留全部拟合状态，未只保留显著项。原始模型原封不动保存在 `original_primary_model/`。

- `estimate` 为 PF−Control 的 β，`std_error` 为 SE；`ci95_lower/upper` 直接取同一版本已保存的绘图模型表 lower/upper，原绘图使用 β ± 1.96×SE 的 95% Wald 区间。没有重新拟合模型或改动已有区间定义。
- `pathway_name/annotation` 来自原 37-LR 表按 `interaction_name` 一对一连接。`feature_id` 为完整键；只有原图已选的 28 条有稳定 `display_id`，其余不新编 I 编号。
- `q_value` 沿用完整检验集合的原 BH 值。`status=ok` 的 singular 拟合仍在主 BH 集合内；其他状态的 q 为 NA。不要在筛出的子表上重新 BH 后称为原结果。
- `standardized_effect_beta_over_SD` 是 β / 含零 log1p(score) 跨切片标准差，不是 fold change。

## 切片得分和补零

`data/slice_scores_nonzero_long.tsv.gz` 保存非零原始 score，键为 sample_id + feature_id。完整特征集合为 `complete_feature_universe.tsv.gz`，切片集合为 `slice_metadata.tsv`。只对这 45×50,040 的笛卡尔积中缺失的键补 0；存在的 score 不变。超出此固定集合的项目不能补零、不能由此推断支持与否。得分尚未 log1p，模型响应为 log1p(score)。本表没有按置换 p 或模型 q 过滤。原输入中未出现的组对/记录按原 prepare_multisample_communication(fill=0) 规则进入模型，零不等于已证明无生物学通信。

原始 45 张切片的已保存 records 已逐键核对，score 差值为 0；非零支持与原模型相同。`historical_slice_score_audit.tsv`、`score_validation_summary.tsv` 保存结果。

## 样本、注释与 pDC

45 切片、35 患者；Control=10 切片/9 患者，Disease（显示为 PF）=35 切片/26 患者。原始 condition 名称保留。提供 sample_affect、run、tma 和每切片各类型细胞数；sample_affect 是原影响程度标签，不能改称精细 tissue-region 坐标注释。未找到与正文输入同版本的切片面积表。

原 RDS 的 `final_CT` 被直接用于分组：对照保存的 45 个 lung_paper fixture 逐细胞验证，无重新聚类/重注释。`provenance/teacher_pipeline_snapshot/metadata_mapping_run.log` 也记录 `cell_type -> final_CT`。这能确认标签传入来源，不能代替原作者构建 final_CT 的 marker/注释算法证据。`pDC_counts_by_slice.tsv` 给出每张切片的 pDC 数量。

检索到的是原有队列图中的 pDC 通信示例，没有找到独立的 pDC marker/表达核验图，也没有找到已经分解 same-cell self 与 different-cell homotypic 贡献的保存结果。pDCs→pDCs 是细胞群层面同型特征，不能仅据此称为跨细胞 crosstalk 或已证实的自分泌。此次未补做 marker 分析或 self 分解。

## 版本边界

主模型仍是 2026-09-07 完成的主分析；0.1.4 是后续安全修复后的推断引擎。`existing_sensitivity_20261012/` 是已经存在的保存得分敏感性分析，包含 run/TMA 调整等，不覆盖主模型，也没有在此次重跑。主模型未调整 run/TMA/sample_affect，不能因元数据中有这些列而称已调整。

99 基因的官方配对材料独立放在 `../01b_lung_99gene_reference_validation/`，不能把其一致性结果当作 343 基因主队列逐一官方配对。保留原 run.log 与 finalize.log：前者末尾的历史校验错误由后者收尾，不能只凭 run.log 认定最终表未完成。
