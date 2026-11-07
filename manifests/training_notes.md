# 推断与降采样 notebook 整理记录

已整理 full、rep01–rep05 六份 GRN 推断 notebook，以及一份五次 counts50 降采样 notebook。本次没有执行这些 notebook、没有重新训练、没有生成新的降采样矩阵，也没有安装软件。每本前置说明与 notebook metadata 都明确标记历史输出；原 execution_count、所有输出块（含进度、警告和图像）逐项保留。原 notebook 和原数据没有修改。

从包根目录 `reproducibility` 运行。所有代码留在 ipynb 内，以 `Path.cwd()` 读取 `manifests/training_manifest.json`。清单中 `input`、`source_input` 默认指向外部历史矩阵；`generated_input` 指向将来重新抽样生成的矩阵。若要使用重新生成的矩阵训练，将相应 run 的 `input` 改成其 `generated_input`。所有新写入均位于包外 `../reproducibility_work`：`data/` 为五份降采样输入，`runs/<run>/` 为推断结果及 `tmp_files`，`reports/downsampling.tsv` 为抽样 QC。输出路径中不再包含原数据目录。

TF 列表读取包内 `reference/T577_macaca_fascicularis.Ensembl91.one2one.tfs.txt`。两份大型 cisTarget 排名库和 motif 注释仍读取原外部路径。原数据库 glob 改为精确的两个排名文件，逐文件调用原 `load_database` API 后合并其返回列表；已只读确认本地 SpaGRN 实现接受单文件路径并返回列表。两库顺序固定为本次整理时原 glob 的实际展开顺序（500 bp 窗口库、10 kbp 窗口库），两库集合与源参考 manifest 一致；历史运行时的 glob 顺序未另行记录，不能据此声称已核验历史顺序。参考 manifest 原声明顺序相反，原 notebook 直接使用 glob，并未使用该声明顺序。

训练代码保持原参数与处理差异：full 使用全部基因，五个 counts50 run 的空间统计排除零方差基因；AUCell 仍使用完整基因背景。Moran、30 邻居、四项 FDR 交集、CPU 16 / pruning 8 workers、`cuda:0`、rank threshold 1500、pruning AUC 0.05、NES 3.0、motif similarity FDR 0.05，以及 AUCell seed 20261014 / threshold 0.05 均未更改。保留历史 notebook 元数据中的原状态，在新增 `reproducibility.execution_status` 中区分本次未执行状态。

降采样函数 `thin_counts` 与原文件逐字一致，采用逐 bin 多元超几何无放回抽样，目标 `floor(total_counts/2)`，五次 seed 20261016–20261020。保留全部 117,867 个 bin 和 15,541 个基因列，不对基因列数或每 bin 检出基因数施加减半约束。原代码对 rep01 只核验历史矩阵、写出其余四份；本副本将五份矩阵都写到外部 `generated_input`。rep01 先对历史 `source_input` 做同 seed 的逐项矩阵比较，发现差异会明确报错。该 I/O 改动不改变抽样函数、随机数调用或抽样次数。保留的原 QC 输出因此仍展示历史输入路径和历史 rep01 只核验行为，不能作为改写后代码已运行的证据。来源计数的 `round(expm1(expr))` 溯源及原始 UMI 未确认声明保留。

已完成的静态检查：七份 notebook JSON/nbformat 校验、全部 code cell AST 语法检查、历史输出与 execution_count 逐项相等、抽样函数逐字相等、五个 seed 顺序、输出路径统一落在包外 work 目录。六份推断 notebook 的 `spg`、`get_modules`、`prune_modules`、`cal_auc` 调用按 AST 比较：除输出文件路径参数外，位置参数及全部关键词参数与原文件相等；原训练文件哈希未变。外部六份输入矩阵、两份排名库、motif 注释和包内 TF 列表均已确认存在，仅检查路径存在，未读取验证大型资源内容。源与副本 SHA256、输出块数、历史 error 输出数和修改的原 cell 索引记录在 `training_notebook_sources.tsv`。这些检查验证结构和保留状态，不验证当前环境能完成训练、不验证大型外部矩阵或参考数据库的内容，也不证明重新执行可数值复现历史结果。清单保留的原矩阵 SHA256 仅为来源记录，本次未重新哈希大型矩阵。

历史输出统计：

| notebook | 历史输出块数 | 历史 error 输出数 | 本次状态 |
|---|---:|---:|---|
| `notebooks/grn/full.ipynb` | 5414 | 0 | 未执行，保留历史输出 |
| `notebooks/grn/rep01.ipynb` | 4484 | 0 | 未执行，保留历史输出 |
| `notebooks/grn/rep02.ipynb` | 5019 | 0 | 未执行，保留历史输出 |
| `notebooks/grn/rep03.ipynb` | 4483 | 0 | 未执行，保留历史输出 |
| `notebooks/grn/rep04.ipynb` | 4454 | 0 | 未执行，保留历史输出 |
| `notebooks/grn/rep05.ipynb` | 4458 | 0 | 未执行，保留历史输出 |
| `notebooks/downsampling/01_counts50.ipynb` | 1 | 0 | 未执行，保留历史输出 |

逐本路径修改：

- `notebooks/grn/full.ipynb`：原 cell 1、2。输入、输出、包内 TF 列表及外部精确数据库路径适配。
- `notebooks/grn/rep01.ipynb`：原 cell 1、2。输入、输出、包内 TF 列表及外部精确数据库路径适配。
- `notebooks/grn/rep02.ipynb`：原 cell 1、2、4、5。输入、输出、包内 TF 列表及外部精确数据库路径适配。
- `notebooks/grn/rep03.ipynb`：原 cell 1、2、4、5。输入、输出、包内 TF 列表及外部精确数据库路径适配。
- `notebooks/grn/rep04.ipynb`：原 cell 1、2、4、5。输入、输出、包内 TF 列表及外部精确数据库路径适配。
- `notebooks/grn/rep05.ipynb`：原 cell 1、2、4、5。输入、输出、包内 TF 列表及外部精确数据库路径适配。
- `notebooks/downsampling/01_counts50.ipynb`：原 cell 1、3、0。增加历史/现行 I/O 说明，外部五份矩阵写入及 rep01 历史矩阵核验。
