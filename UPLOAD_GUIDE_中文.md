# GitHub 补充材料放置说明

适用仓库：<https://github.com/Blake-Deng/SpatialESS>

本包针对 2026-10-22 核对时的公开 main：
`f92dff98f9019300f8fe5086d8a011eed4f089d8`，软件版本 0.1.4。
本次整理没有更改 R/C++ 引擎、重新计算显著性或推送 GitHub。

## 上传到哪里

解压后，将本目录内部的文件按相对路径合并到仓库根目录，
不要把整个 `SpatialESS_GitHub_addendum_20261022` 再套一层传上去。

| 内容 | 仓库目标位置 | 操作 |
| --- | --- | --- |
| 最新 Figure 1，包括 Atera 和邻近图分析 | `paper/figure_package/updates/Figure1_Atera_20261016/` | 新增；原始流程图资产不改 |
| 最新 Figure 2，包括 I15/I06/I19 与分组非零数 | `paper/figure_package/updates/Figure2_20261012/` | 新增并作为当前主图；旧主图保留作归档 |
| 当前主图入口 | `paper/figure_package/README.md` | 替换为本包索引，指向上面两个最新目录 |
| Atera 预处理、抽样、参数、版本、数据库和运行记录 | `paper/atera_validation/` | 新增 |
| 肺队列完整主模型、逐切片分数、样本和 pDC 数量 | `paper/lung_cohort/primary_343gene37LR/` | 新增 |
| 独立的 99 基因官方比较 | `paper/lung_cohort/reference_99gene_validation/` | 新增，与 343 基因主模型分开放 |
| 消融方法、源代码、170 次运行和数值验证 | `paper/ablation/` | 新增；实验变体不覆盖根目录引擎 |
| RCTD 环境证据与原始环境缺失说明 | `paper/provenance/ovarian/` | 合并；已有 DESCRIPTION 内容一致 |
| 全部材料总索引 | `docs/MANUSCRIPT_MATERIALS_20261022.md` | 新增 |
| 入口与校验说明 | `README.md`、`docs/SUBMISSION_CONTENTS.md` | 使用本包提供的更新版 |
| RDS 上传例外和 R 构建排除 | `.gitignore`、`.Rbuildignore` | 使用本包提供的更新版 |
| 更新材料校验记录 | `paper/provenance/github_update_20261022/` | 新增 |
| 合并后的整仓库校验清单 | `SHA256SUMS.txt`、`SUBMISSION_SHA256.tsv` | 使用本包提供的更新版 |

补图 S4 已与最终核对包一致，不需再次替换。S4D 仍为 I03/I21/I19，
与主图 2D 的 I15/I06/I19 不同，图注必须按实际图说明。
现有 `paper/cohort_sensitivity/` 也保留。

## 核验后提交

建议在独立克隆中合并。本包文件已经按仓库层级摆好；GitHub 页面上传
时应先解压，或通过 Git 一次性提交。仅上传 ZIP 不能让读者直接浏览
材料，也不会更新现有图和 README。

合并后从仓库根目录运行：

```sh
python3 paper/provenance/github_update_20261022/verify_addendum.py
python3 paper/verify_manifest.py
sha256sum --check --quiet SHA256SUMS.txt
git status --short
```

确认变更后再提交、推送。本包没有替你执行远程写入。
如果此后 GitHub main 有新变更，先核对 README、忽略规则和校验清单，
不要用旧快照覆盖后续修改。

## 仍需作者确认的来源问题

- 卵巢归档的 DESCRIPTION 为 spacexr **2.2.1**，草稿写的 **2.2.16**
  与其不符；原 RCTD `sessionInfo()` 未找到。不能用后来的检查环境
  代替原始注释环境，也不能未经证据确定原运行精确版本。
- 肺队列 `final_CT` 的原始 marker/人工注释过程仍未完整归档；本包有
  原输入标签转移记录及 pDC 数量，但没有独立 pDC marker 验证图或
  同细胞/不同细胞同型贡献拆分。
- Atera 的具体 CPU/物理内存快照收集于 2026-10-20。若原文要写实验
  当时的硬件型号，仍需确认期间硬件未变，不能把该快照说成逐运行记录。
- condition-specific/consensus 尚无预先确定阈值；现有支持数应作为
  描述性统计，不把对照零值写成生物学缺失。
- 如投稿要求固定永久归档，应在完成这次 GitHub 更新后建立相应 release
  并保存 Zenodo 等持久标识；本包不虚构 DOI。

原始大矩阵、稿件 Word 和本轮主文未使用的 TSTA/SecAct 探索未重复放入。
输入下载来源、校验值、抽样 barcode、配置和代码保留在对应实验目录。
