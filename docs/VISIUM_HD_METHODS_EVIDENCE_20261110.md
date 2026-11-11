# Visium HD 不同 bin size 的 RCTD：方法依据与发表边界

核查日期：2026-11-10。依据为 10x 官方文档、固定版本官方源码、RCTD 原始论文/作者说明及同行评议的 Visium HD 研究。本文区分软件支持、数学等价条件与真实数据效果；合成测试不代表后两项均已完成。

## 结论与推荐流程

**将目标 bin 内的原始 UMI 按基因求和，再在该分辨率重新拟合 RCTD，原理有依据。** 32 µm 是 Space Ranger 正式支持的自定义方形 bin；在原有 `spaceranger count` 命令中增加 `--custom-bin-size=32` 可请求官方输出，其他参考、探针过滤和组织检测设置应保持一致。[官方 CLI](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/running-pipelines/command-line-arguments)

用于文章时，建议优先使用官方原生 `square_032um`。若使用本包的计数聚合作为替代，应先对同一样本的官方输出逐 feature ID、逐 barcode 检查整数计数、完整组织掩膜及 feature 集合。没有原生对照时，不能将数学推导或合成测试写成“已实证证明与官方输出完全一致”。

## 一手依据

| 来源 | 可支持的结论 | 不支持的扩展结论 |
| --- | --- | --- |
| [10x 自定义分箱指南](https://www.10xgenomics.com/analysis-guides/segmentation-visium-hd) | 方形 bin 的表达是来源 2 µm squares 的逐基因 UMI 累加。 | 不能证明某一组织下 32 µm 比 8/16 µm 更准确。 |
| [10x 官方 binning 源码，固定提交 974e49c](https://github.com/10XGenomics/spaceranger/blob/974e49c1b641fa9c1d48d24989a7d5aa72e0b4ca/lib/rust/spatial/src/stages/bin_count_matrix.rs) | 用整数行列除法映射网格，累计计数；组织掩膜按任一来源 square 在组织内传播。raw 矩阵与完整位置表的 barcode 集合可以不同。 | 只保存计数矩阵中出现的位置，不能保证保留正确的粗 bin 组织掩膜。 |
| [官方 filtered matrix 写入源码](https://github.com/10XGenomics/spaceranger/blob/974e49c1b641fa9c1d48d24989a7d5aa72e0b4ca/mro/rna/stages/spatial/write_binned_matrices/__init__.py)；[矩阵定义](https://www.10xgenomics.com/support/software/space-ranger/latest/tutorials/outputs/space-ranger-feature-barcode-matrices) | filtered 除了组织 barcode 过滤，还可能限制探针对应的 target feature。 | raw 全部基因求和不自动等同于官方 filtered bin32；按 symbol 筛选也不能代替按 feature ID 筛选。 |
| [Cable et al., Nature Biotechnology 40, 517–526 (2022)](https://doi.org/10.1038/s41587-021-00830-w)；[原版 SpatialRNA](https://github.com/dmcable/spacexr/blob/ccc43442cd5b4e3b05cddeeff32afdfce0d9d483/R/SpatialRNA.R) | RCTD 对原始计数建模，估计像素内混合贡献；原文包含多类型混合模拟。 | 原论文早于 Visium HD，不是 HD bin32 准确率的验证。 |
| [de Oliveira et al., Nature Genetics 57, 1512–1523 (2025)](https://doi.org/10.1038/s41588-025-02193-3) | 已发表的 HD 8 µm RCTD 应用，并结合独立技术；小细胞混入同一 bin 会影响注释。 | bin 主导标签不是单细胞身份，8 µm 的结果不能直接代替 32 µm 验证。 |
| [Guo et al., Nature Communications 15, 10413 (2024)](https://doi.org/10.1038/s41467-024-54685-1) | Methods 生成了 2/8/16/32 µm 数据，最终选择 16 µm 并使用 RCTD。 | 不是 RCTD 在 32 µm 最优的证据。 |
| [Bioconductor 作者的 HD 工作流](https://bioconductor.org/books/release/OSTA/pages/seq-workflow-visium-hd-bin.html) | 16 µm 计数重新拟合 RCTD；模式和 QC 阈值是显式的数据分析选择。 | 不存在通用的“8 必须 doublet、32 必须 full”或所有分辨率统一 100 UMI 的已验证规则。 |
| [RCTD 原作者的权重单位说明](https://github.com/dmcable/spacexr/issues/108) | 权重表示各类型的 UMI/count 相对贡献。 | 未校正细胞 RNA 含量时，不应直接写成细胞数占比或绝对细胞数量。 |

## 计数等价的条件与推导

同一个 run、同一原始网格原点、同一 feature ID 集合和完整 raw 计数来源下：

\[
Y^{32}_{g,b}=\sum_{u\subset b}Y^{2}_{g,u}
=\sum_{v\subset b}Y^{8}_{g,v}
=\sum_{w\subset b}Y^{16}_{g,w}.
\]

即 16→32 汇总 2×2 来源 bins，8→32 汇总 4×4 来源 bins。这是嵌套整数分箱和求和结合律的推导，**不是已经完成真实样本的实证比较**。基因筛选必须先按 feature ID 进行，再映射 symbol、合并重复 symbol；否则被排除的 feature 可以通过同名基因重新混入。

空间计数矩阵可能不含无 reads 的来源 barcode，但完整位置表仍有它的 `in_tissue`。因此组织掩膜必须由完整位置表单独传播。只有 filtered 来源时，来源组织外 barcode 已被去掉，无法恢复粗 bin 内应纳入的背景计数。修订版默认拒绝这种自动补全，只有显式允许时才生成标注为不完整的自定义聚合。

**raw 矩阵按相同 feature 集合聚合前后应保持 UMI 总量；不同分辨率的 filtered 总量则不必相同。** 粗 bin 只要与组织有交集就可能保留，而该 tile 中组织外来源的计数也随之进入，不能用 filtered 总量必须相同作为通用验收条件。[官方组织规则](https://www.10xgenomics.com/support/software/space-ranger/latest/algorithms-overview/gene-expression)

## RCTD 与下游通讯分析的解释

- 每种分辨率的目标原始计数分别输入 RCTD；不平均细 bin 的拟合权重来冒充粗 bin 拟合。RCTD 是带平台效应的非线性计数模型，此做法是据其模型作出的分析选择。
- `full` 允许多种参考类型贡献；`doublet` 最多拟合两类；选择应符合组织与 bin 的混合程度，并做敏感性检查。增大 bin 不自动提供正确的细胞数模型。
- 本包的默认 UMI/基因下限是运行默认值，不能直接称为论文验证过的最优 QC。可通过 `fit_args_by_bin` 逐分辨率指定，保留解析后的参数和上游有效配置。
- doublet 返回的主导标签/单类型 one-hot 是由分类得到的注释权重；原始 `results_df` 和 `weights_doublet` 单独保留，不能把派生注释当作未经改动的全部 RCTD 权重。
- 可选 SpatialESS bridge 以观测到的混合 bin 表达和主导 RNA 贡献类型进行 **bin 级** 分析。它没有拆出类型特异表达，也不能由 32 µm 中心坐标恢复 bin 内细胞接触。不要将这部分表述为经过验证的细胞级或去混合通讯推断。

## 服务器上仍需完成的真实数据验证

1. 对同一 run 用本包从 raw8 生成 16，与已有官方原生 16 比较；若可获得官方 custom32，再比较 8→32 或 16→32。按 feature ID 验证计数、barcode 集合及完整组织掩膜。
2. 同一单细胞参考、相同注释体系、明确 QC/模式分别拟合 8/16/32；报告保留 bins、UMI/基因深度分布、未拟合/拒绝/不确定比例、参考实际保留细胞数和运行资源。
3. 结合 H&E、已知 marker 或独立成像/测序技术检查生物学合理性。使用参与拟合的 marker 只能提供内部一致性证据，不能替代独立真值。
4. 检查主要细胞类型和下游结论对 bin size、QC、RCTD 模式及参考选择是否稳定。不要把同一样本的不同 bin size 当作独立生物学重复。

已提供原生输出核对命令，**它要求真实存在的原生目标目录，不会用自动回退与自己比较**：

```bash
# 先核对通常已经存在的官方 16 µm 输出。
Rscript scripts/validate_visium_hd_binning.R /data/sample/outs 8 16 check_8_to_16.tsv filtered

# 有官方 --custom-bin-size=32 输出后再核对 32 µm。
Rscript scripts/validate_visium_hd_binning.R /data/sample/outs 16 32 check_16_to_32.tsv filtered
```

没有真实 Space Ranger 对照数据时，以上检查保持“待验证”，不得写成已通过。

## Methods 英文草稿：仅在实际完成相应步骤后使用

> Native 8-, 16- and 32-µm square-bin outputs were generated using Space Ranger [version], with the 32-µm output requested using `--custom-bin-size=32`. Untransformed Gene Expression UMI counts were analysed independently at each resolution using RCTD (spacexr [version/commit]) and a common annotated single-cell reference. The feature inclusion policy, tissue masks, quality-control thresholds and RCTD mode were recorded for each resolution. Estimated weights were interpreted as relative RNA/UMI contributions rather than direct cell-number fractions.

如果实际采用本包聚合，应将第一句替换为实际来源分辨率、网格映射和 feature-ID/组织掩膜处理，并在**完成真实对照之后**报告最大计数差异及掩膜一致性；未运行的实验或未获得的结果不得写入方法或结果。
