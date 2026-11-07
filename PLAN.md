# T577 六图复现与 GRN 结果整理方案

目标目录：`/data/hulei/HongtaoYu_Lab/Fetal_brain/Revision/GRN/reproducibility`。

状态：已按方案完成复制、路径适配和轻量复现验收；推断/降采样及外部矩阵分析保留历史输出，未重新训练。交付与验证记录见 README.md 和 manifests/validation.json。

## 交付范围

整理当前确认的六张图：`Relative_edge_count`、`Relative_regulon_count`、`Shared_regulon_Pearson`、`TF_expression_module_recovery`、`TF_expression_target_recovery`、`Topk_Jaccard_top100_subset`。每张保留 PNG、原生 PDF、绘图代码、实际 cell 输出和全部直接绘图输入。

同时收录完整组与五次 counts50/lowdown 的最终 GRN、小型原生结果文件、降采样代码、每次 GRN 推断代码，以及生成绘图表所需的定量分析代码。以当前 bin100 和五次重复为准，不混入旧 cellbin 流程、旧四分组图、旧固定-signature AUCell 图或三张统计标注图。

原目录文件保留，采用复制整理。包内默认单文件不超过 **10 MiB**，目标总体约 **50–70 MiB**；最终依据实际清单确认，避免超过 100 MiB。文件大小按原始字节统计，不把大矩阵压缩、切片或软链接后伪装为已随包提供。

## 建议目录

```text
reproducibility/
├── PLAN.md
├── README.md
├── T577_SpaGRN_visualization.ipynb
├── notebooks/
│   ├── downsampling/           # counts50 生成代码及实际输出
│   ├── grn/                    # full 和 rep01–rep05，共六次推断
│   └── analysis/               # 网络、Top-k、原生 Pearson、表达分组
├── data/
│   ├── figure_inputs/          # 绘图小表与必要的分析中间表
│   └── grn/
│       ├── full/
│       ├── rep01/
│       ├── rep02/
│       ├── rep03/
│       ├── rep04/
│       └── rep05/
├── figures/                    # 当前六图，各一份 PNG 和 PDF
├── reference/                  # TF 列表、参考资源配置与版本说明
├── environment/                # 推断/绘图两个环境版本及字体信息
└── manifests/                  # 文件、样本、外部大文件和图表来源清单
```

## 1. 六图的输入与代码

| 图 | 必须保留的直接绘图输入 |
|---|---|
| Relative edge count | `network_runs.tsv`、full 的 `regulons.json` |
| Relative regulon count | `network_runs.tsv`、full 的 `regulons.json` |
| Shared-regulon Pearson | `native_shared_regulon_pearson.tsv` |
| TF expression module recovery | `tf_expression_quintile_group_runs.tsv` |
| TF expression target recovery | `tf_expression_quintile_group_runs.tsv` |
| Top-k Jaccard | `network_topk_top100_subset.tsv` |

保留 `network_regulons.tsv`、Top100 逐 regulon 表、原生 Pearson 汇总/质量表、`full_gene_expression_summary.tsv`、五分组 TF 名单及全部 TF–repeat 端点等必要小表，支持追溯各绘图数值的形成过程。

绘图 notebook 从现有根目录版本提取六图及字体/PDF设置，保留实际图片输出，移除不属于这六图的历史绘图部分。保留 Arial、PDF/PS fonttype 42、16 pt 字体、4.5×4.5 英寸、2 pt 左下边框和 3 pt 刻度线。表达两图保持 cividis，Top-k 保持恢复后的紫色渐变。

包内副本改为相对路径，统一从包根目录运行。原路径和原文件校验值写入清单；原 notebook 不改动。单纯修改路径后的历史输出会注明来源，只有实际重新执行的 notebook 才标为本次验证完成。

## 2. 原始 full 与 lowdown GRN 结果

每个 run 保存：

- `regulons.json`：原生最终模块及靶基因集合，单文件约 0.55–1.14 MiB。
- `motifs.csv`：原生 motif pruning 结果，单文件约 5.4–7.7 MiB，六份合计约 35.7 MiB。
- `final_edges.tsv.gz`：从原生 motif/regulon 结果导出的轻量最终边表，列为 `run, regulon, TF, target, weight`，明确权重来源和导出代码。

`final_edges.tsv.gz` 表示最终网络，不冒充原始空间候选邻接表；不额外生成边级 p 值或 FDR。保留最终网络现有自指边口径，并核对 full 的 260 个 regulon、68,492 条边及五次结果计数。

## 3. 降采样、推断及定量代码

降采样以现有 `01_counts50.ipynb` 为主：完整保留 117,867 个 bin 和 15,541 个基因列；各 bin 总 counts 抽样到 `floor(total_counts/2)`，不强制检出基因数减半；记录五次 seed 及 QC。现有代码对 rep01 使用历史结果一致性核验，这一点在说明中保留。

GRN 推断保留 full、rep01 和 rep02–rep05 的六份已执行 notebook。参数、物种参考、常数基因处理、AUCell 设置及原有输出均可追溯。只做路径适配，不更改推断算法、不重新训练。

定量部分保留网络/Top-k 比较、原生共有 regulon AUCell Pearson、完整组表达特征生成及五分组恢复分析。原生 Pearson 使用各自推断模块的评分，在相同 117,867 个 bin 上配对；模块缺失不删除，恢复率记 0、Pearson 保持不可计算。与六图无关的旧分析不混入当前主流程。

## 4. 大文件仅登记，不复制

目前已确认的外部大文件包括：

| 类型 | 当前单文件大小 | 处理方式 |
|---|---:|---|
| full bin100 输入 H5AD | 290.6 MiB | 登记外部路径 |
| 五份 counts50 输入 H5AD | 约 186.9 MiB | 登记外部路径 |
| 六份完整推断结果 H5AD | 约 619–918 MiB | 登记外部路径 |
| 原始空间候选邻接表 | 约 141–224 MiB | 登记外部路径 |
| 全 bin AUCell 大矩阵 | 数百 MiB | 登记外部路径 |
| cisTarget 排名数据库 | 每份约 549–550 MiB | 登记路径、版本和已有 SHA256 |
| motif 注释、基因组及其它大型参考 | 约 15–811 MiB | 登记外部路径 |

外部清单记录文件用途、所属 run、绝对路径、字节数、校验值/校验状态以及下载或构建来源。已有哈希注明其来源，未重新验证的哈希不标为已验证。

小型 TF 列表、参考资源 manifest、相关自定义构建代码可随包保留。记录 SpaGRN 版本、关键本地源码校验信息及参数；如存在本地修改，保留差异。绘图与推断环境分开记录，不复制整个 Conda 环境。保留当前使用的 Arial 字体文件或其明确来源和加载配置，避免依赖临时缓存。

## 5. 复现能力与边界

**随包可完成：** 不访问原 GRN 工作目录的小文件即可重画六图；利用随包最终网络与分析表复算网络规模、模块/靶边恢复率、Top-k 及表达五分组结果。

**需要外部大文件：** 从表达矩阵重新降采样和训练 GRN；从全 bin 原生 AUCell 矩阵重新计算 Pearson；从原始表达矩阵重新计算全基因表达特征。保留这些操作的完整代码和输入路径，但不声称轻量包已包含全量训练/评分数据。

这次整理不会新跑 GRN。完整结果 H5AD 虽然不复制，其路径和与随包小型最终网络的关系会写清楚。

## 6. 执行顺序与验收

1. 固定六图与 full/rep01–rep05 清单，生成源文件和外部文件 manifest。
2. 复制必要小文件、推断/降采样 notebook 和六图；导出最终边表；整理环境、参考及字体记录。
3. 在副本中适配路径，生成六图来源对应表和 README 的运行顺序。
4. 只重新执行绘图及必要轻量复算，保留真实输出；比对 full 计数、五次重复、每组 52 个 TF、898 个有效 Pearson 和固定 Top100 子集等关键值。
5. 校验包内文件大小、校验和、链接和外部依赖边界；确认没有误复制 H5AD、大 AUCell、原始大邻接表或排名数据库。

验收结果包括：六张 PNG/PDF；可执行的六图 notebook；全部直接绘图输入；六组最终 GRN 与权重来源；降采样和六次推断代码及历史输出；README；图→代码→输入表→原生结果的来源链；缺失于包内的大文件清单。校验和只用于文件完整性检查，不代替数值或科学口径核验。
