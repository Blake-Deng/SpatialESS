# Visium HD：bin8、bin16、bin32 的 RCTD 工作流程

修订版 0.1.4.9002 增加了严格的原始计数与 feature-ID 策略、完整位置表组织掩膜，
并保留 Space Ranger 方形 bin 数据读取、稀疏计数聚合、逐分辨率
RCTD 反卷积，以及反卷积结果到 SpatialESS 的接口。`bin_size = 8` 表示
8 × 8 µm，而不是 8 个原始条码。各分辨率分别拟合；不会把 bin8 的权重平均
成 bin16/bin32 的结果。

## 1. 安装

在 SPARKLE 分支源码目录运行：

```bash
Rscript -e 'install.packages(c("Rcpp", "Matrix", "testthat", "remotes"))'
# SPARKLE 分支：先安装内置引擎，再安装 H5AD bridge。
R CMD INSTALL vendor/SpatialESS
R CMD INSTALL .

# 只有读取 H5 和 Parquet 时才需要对应包。
Rscript -e 'install.packages(c("hdf5r", "arrow"))'

# 本接口使用原版 dmcable/spacexr API。
Rscript -e 'options(timeout = 1800); remotes::install_github("dmcable/spacexr@ccc43442cd5b4e3b05cddeeff32afdfce0d9d483", upgrade = "never")'

Rscript scripts/run_visium_hd_rctd.R --help
```

本地验证使用 R 4.5.2，真实服务器验证使用 R 4.5.3；均使用原版 spacexr 2.2.1（上方固定提交）。完整切片计数与真实 ROI 数值对照见 [服务器验证摘要](../vendor/SpatialESS/docs/visium_hd_validation/20261110/README.md)，其中明确列出 native32 和全切片 RCTD 尚未完成的检查。

原版 spacexr 必须导出 `SpatialRNA`、`Reference`、`create.RCTD` 和
`run.RCTD` 等接口。Bioconductor 中采用 `createRctd` / `runRctd` 的同名包
接口不同，不能直接替换。读取 MEX 计数加 CSV 坐标、或者直接传入矩阵和
坐标表时，不需要 `hdf5r` / `arrow`。原版 spacexr 的其他依赖按其安装说明准备。

## 2. 准备单细胞参考

```r
library(Matrix)

# ref_counts：gene × reference-cell 的原始 UMI 计数稀疏矩阵。
# 注释必须有细胞名，且与矩阵的列名一一对应；函数会按名称对齐。
# 例如从已有参考对象提取计数和注释后：
ref_types <- setNames(as.character(reference_metadata$cell_type),
                      rownames(reference_metadata))
saveRDS(list(counts = ref_counts, cell_types = ref_types), "reference.rds")
```

空间数据和参考需要使用同一种基因标识体系（例如都使用 gene symbol）。
不要输入 log-normalized、CP10K 或已取整的校正表达值。参考类型的最小细胞
数量等约束仍由 spacexr 检查，可通过 `create_args` 显式调整其支持的参数。

## 3. 准备配置并运行

Space Ranger 的典型目录是：

```text
sample/outs/
  binned_outputs/
    square_008um/
      filtered_feature_bc_matrix.h5  # 或 filtered_feature_bc_matrix/ MEX
      raw_feature_bc_matrix.h5       # 原生目标读取可选；自动聚合时必需
      spatial/tissue_positions.parquet
    square_016um/
      ...
```

`data_dir` 可指向运行目录、`outs`、`binned_outputs`，或单个
`square_NNNum` 目录（实际命名如 `square_008um`）。单分辨率目录只可读取
该分辨率或向更粗的整数倍分辨率聚合；不能从 bin16 推算 bin8。

修改 [配置示例](../examples/visium_hd_rctd_config.R) 中的服务器路径，然后：

```bash
Rscript examples/visium_hd_rctd_config.R
Rscript scripts/run_visium_hd_rctd.R config.rds
```

配置保存为 RDS 数据；命令不会执行配置文件内的代码。相对路径以
`config.rds` 所在目录为起点。默认拒绝向非空输出目录写入，可更换目录；
明确设置 `overwrite = TRUE` 才允许覆盖同名结果文件。

方法依据、代码审计与论文表述边界见 [方法核查记录](VISIUM_HD_METHODS_EVIDENCE_20261110.md)。

关键参数：

| 配置项 | 默认 / 示例 | 含义 |
| --- | --- | --- |
| `bin_sizes` | `c(8,16,32)` | 逐个独立读取并拟合的 bin 边长，单位 µm |
| `doublet_mode` | `"full"` | 允许多种类型混合；也可显式选 `"doublet"` 或 `"multi"` |
| `min_umi` | `100` | 输入 bin 的 UMI 下限；低深度 bin8 可试较低值并检查结果 |
| `min_genes` | `10L` | 输入 bin 的非零基因数下限 |
| `max_cores` | `8L`（示例） | RCTD 核数，按服务器分配资源设置 |
| `batch_size` | `5000L` | 固定全局拟合后，每批像素拟合的 bin 数 |
| `seed` | `1L` | RCTD 随机种子 |
| `save_inputs` | `TRUE` | 保存实际读取/聚合后的稀疏计数与坐标，便于复查 |
| `keep_fit` | `FALSE` | 默认不保留较大的内部拟合对象 |
| `read_args` | `list(...)` | 传给所有分辨率的读取参数 |
| `read_args_by_bin` | `list(bin32 = list(...))` | 对单个目标分辨率覆盖读取参数 |
| `fit_args_by_bin` | `list(bin16 = list(min_umi = 400))` | 显式覆盖该分辨率的拟合/QC参数；示例不是通用最优阈值 |

`create_args` 和 `run_args` 只用于传递底层 spacexr 支持的参数；核心包装器
已经管理的参数不能重复指定。`run_args` 支持 `MAX_MULTI_TYPES`、
`CONFIDENCE_THRESHOLD`、`DOUBLET_THRESHOLD`、`MIN_CHANGE_REG` 四项拟合
配置；其余可调整的建模阈值通过 `create_args` 传入。`fit_args_by_bin` 可按分辨率覆盖这些参数，实际配置同时写入每个 bin 目录的 `resolved_fit_args.rds`。输入过滤和 spacexr 自身过滤都会留在逐 bin
注释中，不会把未拟合的 bin 伪装成成功预测。

## 4. bin32、组织掩膜与基因集合

用于文章时优先使用同一 Space Ranger 流程通过 `--custom-bin-size=32` 生成的
原生 `square_032um`。如果缺少原生目标，则使用可整除目标边长的较小 bin 的
**raw** 计数求和，例如 16→32 按 2×2 网格分组。基因表达不取平均、不先做
CP10K/log 变换，也不以细 bin 的权重平均代替目标分辨率重新拟合。

读取默认 `feature_policy = "filtered"`：原生目标读取 filtered；重分箱从 raw
取计数，但先用来源 filtered 矩阵的 feature ID 元数据保持相同探针/基因排除
规则，再映射 gene symbol。只有 metadata 的 features.tsv/H5 也可用，无需加载
第二份完整表达矩阵。默认同时要求 raw 来源和可核查的 filtered feature 元数据。
同名 symbol 对应的不同 ID 不能在筛选前合并，否则会重新引入被排除的 feature。

若分析有意使用全部 raw Gene Expression features，显式设置
`read_args = list(feature_policy = "raw")`，在所有分辨率一致使用 raw；这与官方
filtered 结果的 feature 定义不同，应在 Methods 说明。缺少 raw 时，默认拒绝
自动粗分箱。`allow_filtered_aggregation = TRUE` 只适用于明确接受来源不完整的
自定义探索，函数会警告并记录标记，不应据此宣称与官方 coarse raw/filtered 输出等价。

组织规则为任一来源位置在组织内，则目标 bin 在组织内，符合官方规则。
完整 positions 表独立存于 `hd$all_positions`，包括未出现在计数矩阵的零 reads
barcode；组织掩膜由完整位置表聚合，不能仅看有计数的列。缺少组织标记时不能
静默执行组织过滤。`source_position_count` 等字段辅助核查网格覆盖，不自动代表
已证明整张 slide 位置表完整或已与原生输出验证一致。

坐标是阵列网格中心的 µm 坐标：
`x = (array_col + 0.5) * bin_size`，`y = (array_row + 0.5) * bin_size`。
图像像素坐标不直接作为物理距离；聚合后的中心也不取幸存来源位置的平均。

raw 计数在相同 feature 集合下聚合前后守恒；不同分辨率 filtered 总 UMI 在组织
边界不必相同。应使用脚本逐 feature ID 核对真实官方原生目标，避免用总量相等
代替逐项正确性：

```bash
Rscript scripts/validate_visium_hd_binning.R /data/sample/outs 8 16 check16.tsv filtered
# 获得官方 custom32 后：
Rscript scripts/validate_visium_hd_binning.R /data/sample/outs 16 32 check32.tsv filtered
```

脚本只接受真实存在的原生目标目录；缺失时失败，不会用自动回退结果与自己比较。
本地合成检查与源码推导不代替该真实数据核查。

## 5. 输出与服务器检查

```text
output_dir/
  resolved_config.rds
  sessionInfo.txt
  run_summary.tsv
  bin8/
    visium_hd_input.rds
    rctd_result.rds
    proportions.rds
    annotations.tsv
  bin16/...
  bin32/...
```

`proportions.rds` 是稀疏的 **bin × reference-type** 相对 RNA/UMI 贡献矩阵（doublet 的 singlet 为派生 one-hot 注释权重），行顺序与输入一致。其单位不是直接的细胞数量或细胞数占比。
过滤、未拟合或拒绝的 bin 保留在审计中，相应权重行可能全零。不要把全零
行理解为“该组织没有细胞类型”。`rctd_result.rds$classification` 另存 doublet 的原始分类表及 pair weights（包括置信/收敛字段），即使 `keep_fit = FALSE` 也保留。`parameters$upstream_config` 与 `diagnostics$reference` 记录上游有效参数及参考构造器实际保留细胞数。`annotations.tsv` 保留 RCTD 分类、主导
类型和最大权重。运行失败会返回非零退出状态，在 `run_summary.tsv` 和
对应目录的 `ERROR.txt` 记录原因，后续分辨率不会继续静默运行。

每个分辨率至少检查：

```r
x <- readRDS("/server/results/sample_visium_hd_rctd_test/bin32/visium_hd_input.rds")
r <- readRDS("/server/results/sample_visium_hd_rctd_test/bin32/rctd_result.rds")
stopifnot(identical(colnames(x$counts), rownames(r$weights)))
rs <- Matrix::rowSums(r$weights)
stopifnot(all(rs == 0 | abs(rs - 1) < 1e-8))
table(r$annotations$status, useNA = "ifany")
summary(r$annotations$max_weight)
x$provenance
r$diagnostics
```

建议先用一个小样本/ROI 跑通，再跑完整切片，并检查低 UMI、未拟合、拒绝
比例和参考类型是否合理。完整切片的运行时间、内存和真实组织准确率需要
在服务器实测；本地合成数据检查不能替代这一验证。

## 6. R API 与可选 SpatialESS 分析

```r
library(SpatialESS)
ref <- readRDS("reference.rds")
hd <- read_visium_hd("/server/data/sample/outs", bin_size = 32)
rctd <- run_visium_hd_rctd(
  hd, ref$counts, ref$cell_types,
  doublet_mode = "full", max_cores = 8L,
  min_umi = 100, min_genes = 10L, batch_size = 5000L, seed = 1L
)

# 已在内存中的原始计数和 Space Ranger 坐标表也可直接使用：
# hd8 <- prepare_visium_hd(counts, positions, bin_size = 8)
# hd32 <- aggregate_visium_hd(hd8, bin_size = 32)

input <- prepare_visium_hd_spatialess(hd, rctd, min_weight = 0.5)
db <- readRDS("lr_database.rds") # list(lr = ..., complex = ..., cofactor = ...)
result <- spatialess_from_visium_hd(
  hd, rctd, lr = db$lr, complex = db$complex, cofactor = db$cofactor,
  min_weight = 0.5, lr_missing = "filter", nboot = 100L, seed.use = 1L
)
```

此接口将通过筛选的 bin 赋予主导细胞类型标签，并对 **bin 级表达与空间关系**
运行现有 SpatialESS。RNA/UMI 混合贡献与原始分类审计单独保留；不把单个 bin 的表达拆成
多个细胞，不实现比例加权的多类型通讯模型。尤其在 bin16/bin32，主导标签
不能解释为已解析的单细胞身份；bin 中心也无法恢复 bin 内或跨边界细胞的真实接触距离。`min_weight` 是下游保留阈值，不是 RCTD
的概率校准。`doublet` 模式会保留其分类语义，拒绝或不确定结果不能当作
可靠主导标签。

包装器在每个分辨率上只估计一次全局平台效应和 sigma，再分批进行像素
拟合，避免对每批重新拟合全局参数。该策略限制的是像素拟合阶段的密集
工作集；参考构建、平台拟合以及最终结果仍需相应内存。不同分辨率改变了
观测单位及计数深度，其权重和通讯结果不能直接当作可互换的重复测量。

## 数据格式与算法依据

- [10x Space Ranger 输出结构](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview)：标准 HD 目录与默认 2、8、16 µm 分辨率。
- [10x 自定义 bin 指南](https://www.10xgenomics.com/analysis-guides/segmentation-visium-hd)：按 bin 逐基因累加 UMI 的数据定义。
- [原版 spacexr](https://github.com/dmcable/spacexr)：本包装器使用的 RCTD 实现与依赖。
- [原版 RCTD 像素拟合](https://github.com/dmcable/spacexr/blob/master/R/runRCTD.R)与[结果整理](https://github.com/dmcable/spacexr/blob/master/R/postProcessing.R)：full/doublet/multi 模式、权重与分类语义。
