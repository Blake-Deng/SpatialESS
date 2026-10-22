# Atera 预处理与运行证据

来源为 10x Atera WTA **preview** FFPE breast/cervical 数据，公开页和下载审计保存在 data_provenance/。这是开发版本预览数据，不应写成最终商业流程或队列生物学验证。

乳腺 170,057→169,940，宫颈 717,576→715,512。只保留 18,028 个 Gene Expression 行；使用这些行计算每细胞 counts/features。保留 counts≥20 且 features≥20、且存在官方 group 标签的细胞。缺注释原始数分别 12、293；与 QC 不合格可能重叠，不能相加推断排除数。原始 excluded_cells.tsv 与 matrix_complete.json 已提供，精确排除组合见 qc_exclusion_combinations.tsv。

表达为 log1p(10,000×count/该细胞全部 measured-GEX counts)。group 直接匹配 10x 的 *_cell_groups.csv 官方 group 字段，不是本项目补做的细胞类型注释。x_centroid/y_centroid 原为 µm，ratio=1，不重复乘 pixel_size。

central 取到所有 eligible 细胞坐标算术平均中心距离最近的 N 个，距离并列按原行顺序，无随机数。stratified 使用全部 eligible group 比例、向下取整后按余数从大到小补齐（并列按排序 group 顺序），组内无放回 sample.int 并最终打乱；三次 seed=20261013/14/15。block 使用 x 的 25/50/75 分位与 y 中位附近的实际 centroid 为锚点，再选最近 N 个。精确算法见 prepare_fixtures.R。200k/400k 延伸保持 central 嵌套；代码与已有校验报告保留。barcode_lists/ 是只保留 cell_id/subset_order 的压缩导出，原完整 barcode 文件的 SHA256 保存在 manifest.tsv。

两引擎输入同一固定 fixture：同表达、坐标、group、LR 顺序、参数、100 次置换及 seed.use=1。contract.rds 保存基因、barcode、group counts 与参数，summary 中有 input_md5；配对检查见 paired_input_audit.tsv 和原 agreement 表。预处理不做 overexpression 筛选；3,234 条数据库记录中有 3,081 条 mandatory-component-complete LR，属于 fixed-LR 比较，不能称完整默认 LR-selection pipeline。

interaction.range=30 µm，tol=5，CSR 候选半径=35 µm；contact.range=10 µm，contact.dependent=TRUE、forced=FALSE，use.AGAN=TRUE，min.percent=.1，min.cells.sr=5，Kh=.5，n=1。scale.distance 按组织固定为 2/minimum positive cell distance，乳腺约 2.04351143947、宫颈约 2.29750351419；两引擎共用。精确距离权重与 self 约定见归档 source/SpatialESS/R、src 和 worker.R，不将 CSR 无 self edge 等同于去除官方 self-score 约定。

CPU 绑定及线程在各 config/runner_environment 里。每 worker 单 CPU、OMP/BLAS/MKL/BLIS=1；96 GiB 是虚拟地址空间限制而非 peak-RSS 上限。最初参考超时 1h，后续 scale/200k/400k/700k 参考超时 12h；ESS 超时 24h。各运行 sessionInfo 和 source_snapshot 哈希随包附上。硬件现状另存 hardware_snapshot_at_collection.json；若缺逐运行硬件型号快照，不将这次采集时间冒充运行时间。

whole-process 来自 GNU time（包含载入/写盘），inference_seconds 是 R 内推断段，峰值 RSS 为进程实测；CSR bytes 与理论 dense bytes 是别的指标。resource_summary_all_source_rows.tsv 刻意保留各阶段来源，复用 ESS 行可跨阶段重复，不当新重复实验；all_attempts_ledger.tsv 另保留历史失败/重试。失败运行的耗时/RSS 是尝试消耗，不当成功性能或加速比。

已找到 200k/400k 数值日志，非只有状态。200k ESS：推断 1308.950 s、全进程 1319.80 s、峰值 3.529719808 GB；400k：2423.322 s、2440.91 s、6.170734592 GB。两个参考尝试均 allocation failure，实测尝试全进程分别 36.49/84.85 s，RSS 33.293033472/65.585041408 GB。700k 类别实际 N=715,512，参考尝试同样失败，31.89 s、12.961091584 GB，申请 96.1 GiB 对象失败。失败原因与原日志保持，不推断无限资源下也必定失败。此次没有重跑。
