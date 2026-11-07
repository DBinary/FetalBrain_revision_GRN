# T577 六图复现与 GRN 结果

本目录包含当前六张图的代码、全部直接绘图输入、full 和五次 counts50/lowdown 的最终 GRN，以及降采样和六次推断代码。原工作目录保持不变。大表达矩阵、完整结果 H5AD 和大型参考数据库只在清单中登记，没有复制进来。

当前版本使用五组 TF 表达分析。表达两图保留 cividis，Top-k 图保留紫色渐变；不包含旧四分组、固定 full-signature 评分图或统计标注图。

## 文件入口

| 路径 | 内容 |
|---|---|
| `T577_SpaGRN_visualization.ipynb` | 六图的全部读表及绘图代码、实际图片输出 |
| `figures/` | 六张最新单图 PNG、六张原生 PDF |
| `data/figure_inputs/` | 已确认的绘图输入与分析基准，不由复算步骤覆盖 |
| `data/grn/full/`、`rep01/`–`rep05/` | 六组原生 regulon JSON、motif CSV、带权最终边表 |
| `notebooks/analysis/` | 最终边导出、网络/Top-k、Pearson、表达特征、五分组、验收 |
| `notebooks/downsampling/`、`notebooks/grn/` | 降采样及六次推断代码，保留历史输出 |
| `environment/` | 两套 Python/包版本、Arial、SpaGRN Python 源码快照 |
| `reference/` | TF 列表、参考配置及相关构建代码 |
| `manifests/` | 来源、参数、文件大小/校验、外部依赖和验证记录 |

## 只使用包内小文件的运行方式

所有 notebook 均从包根目录运行。以下命令显式设置工作目录，避免 Jupyter 默认进入子目录：

```bash
git clone https://github.com/DBinary/FetalBrain_revision_GRN.git
cd FetalBrain_revision_GRN
python - <<'PY'
from pathlib import Path
import nbformat
from nbclient import NotebookClient

root = Path.cwd()
kernel = "hulei_spatrack"  # 本机绘图/定量环境；迁移后替换为相应 kernel
for name in [
    "notebooks/analysis/02_network_topk.ipynb",
    "notebooks/analysis/05_TF_expression_quintiles.ipynb",
    "T577_SpaGRN_visualization.ipynb",
    "notebooks/analysis/06_validate_package.ipynb",
]:
    path = root / name
    notebook = nbformat.read(path, as_version=4)
    NotebookClient(notebook, kernel_name=kernel, timeout=600).execute(cwd=str(root))
    nbformat.write(notebook, path)
PY
```

只重画图片时，执行根目录绘图 notebook 即可。它只读取包内表和字体，不读取外部大矩阵。Arial 显式注册并从包内加载，不依赖临时字体缓存或联网下载。字体为 16 pt、图幅为 4.5×4.5 英寸，PNG 为 600 dpi，PDF 嵌入 TrueType 字体。

轻量复算写入 `data/analysis_recomputed/`，与 `data/figure_inputs/` 的确认快照对照。`01_export_final_edges.ipynb` 从 motif 结果重导权重，需要 PySCENIC/SpaGRN 环境；带权边表已经提供，日常轻量复现可跳过。两套环境版本见 `environment/core_package_versions.tsv` 和 pip freeze 文件；版本记录不代表已在另一台机器重新安装验收。

## 六图与输入口径

| 图 | 包内直接输入 | 点的含义 |
|---|---|---|
| Relative edge count | `network_runs.tsv`、full `regulons.json` | 一次重复的总边或共同边 / full 的 68,492 条边 |
| Relative regulon count | 同上 | 一次重复的总模块或共同模块 / full 的 260 个模块 |
| Shared-regulon Pearson | `native_shared_regulon_pearson.tsv` | 一个共有模块的原生 AUCell Pearson |
| TF expression module recovery | `tf_expression_quintile_group_runs.tsv` | 一次重复、该组保留模块数 / 52 |
| TF expression target recovery | 同上 | 一次重复、该组保留的原始边 / 该组全部 full 边 |
| Top-k Jaccard | `network_topk_top100_subset.tsv` | 一次重复固定 Top100 子集的 Jaccard 中位数 |

Pearson 比较两组各自推断 regulon 的原生评分，不是固定 full-signature 评分。每个相关系数对应全部 117,867 个配对 bin。五次共有模块数为 183、181、181、178、175，共 898 行；full 260×5 配对中另有 402 个缺失模块，其恢复率为 0、Pearson 不可计算。负相关值保留。

表达组只由 full 的平均 CP10K 定义，零表达 bin 计入均值。每组 52 个 TF，图按 Q5→Q1 排列。靶基因恢复按 TF–target 边计数，同一个 target 在不同 TF 下分别计数；不要与逐 TF 等权恢复比例混淆。

Top-k 按原生权重降序、target 名升序排序，不取绝对值。双方最终靶标数均不少于 100 的固定子集依次为 67、60、61、64、62 个；同一 repeat 的五个 k 使用相同名单，跨 repeat 名单不同。这是条件性比较，不能单独代表全部模块的稳定性。

## 原生 GRN 与历史推断代码

每组 `regulons.json`、`motifs.csv` 是原生结果的复制件。`final_edges.tsv.gz` 从 motif 的 `gene2weight` 导出，并严格保留 JSON 最终成员；列为 `run, regulon, TF, target, weight`。自指边及原生负权重保留。它是最终网络，不是原始空间候选邻接表，weight 也不是边级 P 值或 FDR。

full 为 260 个模块、68,492 条边。六份推断 notebook 保留历史实际输出，仅适配路径，本次没有重新训练。来源及原文件 SHA256 见 `training_notebook_sources.tsv`；I/O 改动、排序来源等细节见 `training_notes.md`。full 和低采样在常数基因处理上的原有差异保留。

## 需要外部大文件的步骤

- 重新生成五份 counts50，或重新训练 GRN。
- 执行 `03_native_shared_Pearson.ipynb`，从六份原生结果 H5AD 的 AUCell 重新计算相关系数。
- 执行 `04_full_gene_expression.ipynb`，从 full 表达矩阵重新计算全基因特征。

后两本分析 notebook 保留原始实际输出，本次未读取大矩阵重新执行。已有 Pearson、配对质量表和全基因表达摘要随包提供。

`workflow_manifest.json` 的 JSON、motif、最终边路径为包内相对路径；`input`、`result`、`adjacencies` 为外部路径。训练/降采样使用 `training_manifest.json`：默认训练读取历史外部矩阵，重新生成输入另列为 `generated_input`，可显式切换。所有新大输出默认写包外 `../reproducibility_work/`，不会覆盖历史文件或塞入本包。

降采样保持逐 bin 无放回多元超几何抽样到 `floor(total_counts/2)`，seed 为 20261016–20261020，不要求检出基因数减半。rep01 的历史一致性核验保留，并可将五份新输入全部写到包外。

大型参考包括两份约 549–550 MiB 的 ranking 和约 57 MiB 的 motif 注释表。包内只有 TF 列表、配置和构建代码；原 reference manifest 中的路径属于来源记录，不代表大文件已经提供。外部文件见 `external_large_files.tsv`，已有哈希与现场复算状态分别注明，未现场重算的哈希不标为已验证。

## 验收与解释边界

`manifests/validation.json`、`06_validate_package.ipynb` 记录实际数值和文件检查；`files.tsv` 记录包内文件大小、SHA256；`figure_sources.tsv` 将图对应到代码与输入。文件清单不包含它自己的校验值，避免自引用。

counts 来自此前重建的整数，未核实为原始 UMI。五次 repeat 是同一 T577 输入的技术降采样，不是五个独立个体。文件校验只能证明完整性；相应数值复算只能证明本次统计口径一致，均不能代替跨个体科学验证。

## 最终版本与 caption

以用户最新确认的六张独立图为准：Top-k 为紫色渐变，两张 TF 表达图为 cividis。英文 caption 为 [CAPTION.md](CAPTION.md)，LaTeX 为 [CAPTION.tex](CAPTION.tex)；Fig. X 为待替换图号。caption 按最终排版 a–f 顺序：regulon count、edge count、Pearson、Top-k、module recovery、target recovery。用户提供的旧色阶拼版仅作为文字排版参考，来源见 `manifests/caption_layout_reference.json`，没有混入包内最终结果。

本次已为外部清单的 35 个文件计算当前 SHA256；其中两份 ranking 还与参考 manifest 的已有哈希核对一致。其余 33 项是本次版本指纹，不声称追溯验证了未记录的历史哈希。大文件未复制，也没有重新推断 GRN 或重算全量评分。

`files.tsv` 是交付时文件快照；重新执行 notebook 会更新输出及文件字节哈希，数值复现通过不要求重跑后 notebook 的字节哈希仍与交付快照相同。参考构建脚本为原始代码归档，重建参考仍需其上游工具和外部资源。
