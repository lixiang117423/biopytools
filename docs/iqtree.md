# IQ-TREE 系统发育树分析 | IQ-TREE Phylogenetic Tree Analysis

一句话理解：**把一堆序列的比对文件交给 IQ-TREE，让它自动挑最优进化模型、算出一棵「最可信的家族谱系树」，并给每个分支一个可信度打分**。

输入一个多序列比对文件（FASTA/PHYLIP/NEXUS），输出一棵带支持值的系统发育树，可选做一致性因子分析和祖先状态重建。

## 功能概述 | Overview { #overview }

- 基于 IQ-TREE 从多序列比对构建最大似然（ML）系统发育树
- 未指定模型时自动做模型选择（ModelFinder，MFP），无需用户懂进化模型
- 内置两种 Bootstrap：UFBoot（默认，快）与标准 Bootstrap，给每个分支算支持值
- 可选功能：外群指定、约束树、分区分析（partition）、基因一致性因子（gCF）、祖先状态重建
- 断点续传：树、一致性因子、祖先重建各自按输出文件存在性跳过（换参数重跑用 `--redo`）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools iqtree -i alignment.fasta -o tree_results -p my_tree
```

最小输入：一个多序列比对文件（FASTA/PHYLIP/NEXUS），`-o` 指定输出目录，`-p` 指定输出文件前缀（所有结果文件都以此前缀命名）。

## 零基础概念速览 | Concepts in plain words { #concepts }

不熟悉进化分析术语的话，先看这张表，后面的参数说明都会用到：

| 术语 | 通俗理解 |
|------|----------|
| 多序列比对 | 把多条同源序列「对齐排好」的文本，像把同一段话的不同版本逐字对齐 |
| 系统发育树 | 物种/序列的「家族谱系图」，亲缘近的靠在一起 |
| 最大似然(ML) | 在「哪种树最能让现有序列出现」里挑概率最大的，本工具的建树方法 |
| 进化模型 | 描述「碱基/氨基酸怎么变」的数学规则；不懂就选自动，程序自己挑 |
| Bootstrap/UFBoot | 把数据重抽样很多遍再看「这个分支还稳不稳」，得分越高越可信（满分 100） |
| 外群 | 已知「最不像大家」的一条，用来给树定方向（哪边是祖先） |
| 分区分析 | 把比对按基因/位点分成几段，各用各的模型，多基因数据更准 |
| 一致性因子(gCF) | 用一堆基因树看「主树里每个分支有多少基因支持」，判断分支可靠程度 |
| 祖先状态重建 | 反推「老祖宗当年每个位点是什么碱基/氨基酸」 |

## 输入 | Input { #input }

一个多序列比对文件，支持 FASTA / PHYLIP / NEXUS 格式（IQ-TREE 自动识别）。序列数至少 3 条，序列名应唯一、避免特殊字符（否则下游可视化可能出问题）。

示例（FASTA）：

```text
>species_A
ATCGATCGATCG
>species_B
ATCGATCGATCC
>species_C
ATCGATCGATCG
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 三个必须给的：输入比对、输出目录、输出前缀。前缀决定了所有结果文件名（如 `my_tree.treefile`），换项目换前缀即可避免混淆。

### 模型与线程 | Model & threads

**通俗理解|In plain words:** 模型参数管「用哪套数学规则建树」。**不懂就完全不写**——不指定时自动用 ModelFinder（MFP）挑最好的，几乎总是最优解；只有你明确要复现某篇文献的模型时才手写。`--sequence-type` 管「把序列当 DNA 还是当蛋白」——默认程序自动嗅探并显式告知 IQ-TREE，一般不用动；只有 NEXUS 等少见格式或形态/二进制数据才需要手填。线程数管并行加速，一般不用动。

### Bootstrap 支持值 | Bootstrap

**通俗理解|In plain words:** 管「给每个分支打几分可信度」。重复次数越大越稳、越慢；默认 1000 够用。UFBoot 比标准 Bootstrap 快得多、结果同样可靠，默认即 UFBoot，一般不用动；`--save-boot-trees` 仅在你想拿到全部 bootstrap 树做后续分析时才开。

### 外群与约束 | Outgroup & constraint

**通俗理解|In plain words:** 外群用来「给树定方向」（哪个分支代表祖先），用逗号分隔多个名称。约束树用来「强制结果符合某个已知拓扑」（比如你确信某两个物种必须在一起），不确信就别给。

### 分区分析 | Partition

**通俗理解|In plain words:** 多基因/多位点数据才需要。分区文件告诉程序「哪几列是一个基因、各用什么模型」。三种模式（edge-linked 默认 / edge-equal / edge-unlinked）区分「不同分区共享多少参数」，默认 edge-linked 最稳，一般不用动。

### 高级功能 | Advanced

**通俗理解|In plain words:** `--concordance` 需要另给一个基因树文件，用来算每个分支的基因一致性因子，做系统发育组学（phylogenomics）时才有意义；`--ancestral` 开启祖先状态重建。两者都不开也不影响建树。

### 其他 | Other

**通俗理解|In plain words:** `--seed` 固定随机种子（保证结果可复现）；`--runs` 独立跑多次取最好的树（默认 1，一般不用动）；`--redo` 强制重算（见 FAQ）；`--iqtree-path` 指定 IQ-TREE 程序路径（默认自动找，一般不用动）。

## 分析流程 | Pipeline { #pipeline }

```text
步骤1: 检查依赖（IQ-TREE 是否可用）
   |
   v
步骤2: 构建系统发育树
   |-- 未指定模型 -> ModelFinder 自动选模型
   |-- 有 bootstrap -> 计算分支支持值
   +-- 输出 {prefix}.treefile / .contree / .iqtree / .log
   |
   v
步骤3: 一致性因子分析（可选，需 --concordance）
   +-- 输出 {prefix}_concordance.cf.tree / .cf.branch / .cf.stat
   |
   v
步骤4: 祖先状态重建（可选，需 --ancestral）
   +-- 输出 {prefix}_ancestral.state / .treefile / _ancestral_sequences.fasta / _ancestral_summary.txt
```

## 输出 | Output { #output }

所有结果都写在你指定的输出目录里，文件名以 `-p` 的前缀开头：

```text
output_dir/
|-- {prefix}.treefile                    # 最佳 ML 树（有精确枝长，供进一步分析）
|-- {prefix}.contree                     # 共识树（带 bootstrap 支持值，供可视化/发表）
|-- {prefix}.iqtree                      # 主结果文件（模型参数、统计信息，文本可读）
|-- {prefix}.log                         # IQ-TREE 运行日志
|-- {prefix}.model.gz                    # 模型选择结果（仅自动选模型时）
|-- {prefix}.ufboot                      # UFBoot 全部树（仅 --save-boot-trees 时）
|-- {prefix}_concordance.cf.tree         # 一致性因子结果（仅 --concordance）
|-- {prefix}_concordance.cf.branch
|-- {prefix}_concordance.cf.stat
|-- {prefix}_ancestral.state             # 祖先状态重建原始结果（仅 --ancestral）
|-- {prefix}_ancestral.treefile          # 带祖先状态的树
|-- {prefix}_ancestral_sequences.fasta   # 重建出的祖先序列
|-- {prefix}_ancestral_summary.txt       # 祖先状态汇总报告
+-- iqtree_analysis.log                  # biopytools 流程日志
```


## 结果解读 | Interpreting Results { #interpreting }

- **`.treefile`**：最佳 ML 树，枝长是精确估计的，用于后续进化分析（如分子钟、选择分析）。新树里节点没有支持值时，用它看拓扑。
- **`.contree`**：共识树，节点上标的数字就是 bootstrap/UFBoot 支持值（0-100）。**数值越高分支越可信**：通常 >=95 很可靠，70-94 中等，<70 需要谨慎（该分组证据不足）。
- **`.iqtree`**：文本报告，含所用模型、参数、log-likelihood、树形统计。论文 Methods 部分「用了什么模型」就从这里抄。
- **`.model.gz`**：ModelFinder 挑出的最优模型（按 BIC 排序），可用来佐证「为什么用这个模型」。
- **一致性因子 `.cf.stat`**：每个分支的 gCF/sCF 值。gCF 高 = 大多数基因树支持该分支；低 = 基因间有冲突（可能是杂交/不完全谱系分选），说明该分支的「物种树 vs 基因树」不一致。
- **祖先状态 `.state` 与 `_ancestral_summary.txt`**：`_ancestral_summary.txt` 汇总了每个位点重建出的状态分布，`_ancestral_sequences.fasta` 给出整条祖先序列，供后续功能分析。

## 参数选择建议 | Parameter Guidance { #guidance }

- **完全新手**：只给 `-i -o -p`，其余全默认（自动选模型 + UFBoot 1000），结果已可直接用。
- **要发表、看支持值**：默认就有 `.contree`，无需额外参数。
- **多基因/多位点数据**：加 `--partition`，配合分区文件（默认 edge-linked）。
- **系统发育组学、基因树很多**：加 `--concordance 基因树文件`，看每个分支的一致性。
- **要祖先序列做功能重建**：加 `--ancestral`。
- **想要标准（非 UFBoot）支持值**：`--boot-type standard`。
- **结果不稳定、想固定随机性**：加 `--seed`；想多跑几次取最优用 `--runs`。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入比对文件｜Input alignment file |
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--prefix, -p` | 必填 |  | 输出文件前缀｜Output file prefix |
| `--model, -m` | — |  | 进化模型｜Evolutionary model |
| `--sequence-type` | — | DNA/AA/BIN/MORPH | 序列类型 (默认按比对字母表自动嗅探并显式传-st, 规避IQ-TREE 3对简并码富集比对的误判)｜Sequence type (auto-sniffed from alignment alphabet by default) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--bootstrap, -b` | `1000` | int | Bootstrap重复次数｜Bootstrap replicates |
| `--boot-type` | `ufboot` | ufboot/standard | Bootstrap类型｜Bootstrap type |
| `--save-boot-trees` | — |  | 保存所有bootstrap树｜Save all bootstrap trees |
| `--outgroup` | — |  | 外群名称(逗号分隔)｜Outgroup taxon names (comma-separated) |
| `--constraint` | — | Path | 约束树文件｜Constraint tree file |
| `--partition` | — | Path | 分区文件｜Partition file |
| `--partition-mode` | `edge-linked` | edge-linked/edge-equal/edge-unlinked | 分区模式｜Partition mode |
| `--concordance` | — | Path | 一致性因子基因树文件｜Concordance factor gene tree file |
| `--ancestral` | — |  | 启用祖先状态重建｜Enable ancestral state reconstruction |
| `--seed` | — | int | 随机种子｜Random seed |
| `--runs` | `1` | int | 独立运行次数｜Number of independent runs |
| `--redo` | — |  | 重新运行分析｜Redo analysis |
| `--iqtree-path` | `iqtree` |  | IQ-TREE软件路径｜IQ-TREE program path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入比对文件 (FASTA/PHYLIP/NEXUS)｜Input alignment file (FASTA/PHYLIP/NEXUS) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-p, --prefix` | 必填 |  | 输出文件前缀｜Output file prefix |
| `-m, --model` | — |  | 进化模型 (不指定则自动选择)｜Evolutionary model (auto-select if not specified) |
| `--sequence-type` | — | DNA/AA/BIN/MORPH | 序列类型 (不指定则按比对字母表嗅探; NEXUS等格式回退IQ-TREE自动检测)｜Sequence type (sniffed from alignment alphabet if omitted; NEXUS etc. fall back to IQ-TREE auto-detection) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `-b, --bootstrap` | `1000` | int | Bootstrap重复次数｜Bootstrap replicates (default: 1000) |
| `--boot-type` | `ufboot` | ufboot/standard | Bootstrap类型｜Bootstrap type (default: ufboot) |
| `--save-boot-trees` | — | store_true | 保存所有bootstrap树到文件｜Save all bootstrap trees to file |
| `--outgroup` | — |  | 外群名称 (多个用逗号分隔)｜Outgroup taxon names (comma-separated) |
| `--constraint` | — |  | 约束树文件｜Constraint tree file |
| `--partition` | — |  | 分区文件｜Partition file |
| `--partition-mode` | `edge-linked` | edge-linked/edge-equal/edge-unlinked | 分区模式｜Partition mode (default: edge-linked) |
| `--concordance` | — |  | 一致性因子分析：基因树文件｜Concordance factor: gene tree file |
| `--ancestral` | — | store_true | 启用祖先状态重建｜Enable ancestral state reconstruction |
| `--seed` | — | int | 随机种子｜Random seed |
| `--runs` | `1` | int | 独立运行次数｜Number of independent runs (default: 1) |
| `--redo` | — | store_true | 重新运行分析｜Redo analysis |
| `--iqtree-path` | — |  | IQ-TREE程序路径｜IQ-TREE program path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- IQ-TREE（默认路径 `~/miniforge3/envs/phylo/bin/iqtree`，conda 环境 `phylo`，可用 `IQTREE_PATH` 环境变量或 `--iqtree-path` 覆盖）
- Python 3（biopytools 运行环境）

## 常见问题 | FAQ { #faq }

**Q1：换参数重跑，结果没变？**
断点续传按输出文件存在性判断——只要 `{prefix}.treefile` 已存在就跳过建树。想真正重算，加 `--redo`（会传给 IQ-TREE 覆盖旧文件）；或删掉输出目录里的旧 `.treefile/.iqtree/.log` 再跑。

**Q2：跑一半中断了，重跑会从头再来吗？**
IQ-TREE 自己有断点机制，会复用已有的中间文件继续。但如果中断留下了不完整的 `.iqtree/.log`，可能报「文件已存在」错误，此时加 `--redo` 强制重来最省事。

**Q3：`-p` 到底是前缀还是分区文件？**
在 `biopytools iqtree` 里，`-p/--prefix` 是**输出前缀**；分区文件用的是 `--partition`（长参数，无 `-p` 短参）。注意这和 IQ-TREE 原生命令不同（原生 `-p` 是分区），别照搬原生用法。

**Q4：contree 里某个分支支持值很低，树还能用吗？**
低支持值只说明「这个分组」证据弱，不影响整棵树。可以用 `.treefile` 看完整拓扑，重点讨论高支持值分支；低支持分支在结论里要谨慎表述。

**Q5：报 `ERROR: Unknown sequence type` 退出？**
IQ-TREE 3 的序列类型自动检测对「简并码/缺失富集」的 DNA 比对会误判——全比对非 ACGT 字符（R/Y/N 等杂合简并码）占比 ≥10% 就直接报错退出（IQ-TREE 2 无此问题）。VCF 转来的群体 SNP 比对恰好杂合多、极易超线。模块现在默认按字母表嗅探并显式传 `-st DNA` 规避；若仍报错（如 NEXUS 输入未被嗅探），手动加 `--sequence-type DNA`。

**Q5：UFBoot 和标准 Bootstrap 哪个好？**
两者都给分支支持值，UFBoot 计算量小很多、结果通常与标准 Bootstrap 高度一致，默认 UFBoot 即可。个别审稿场景要求「标准 bootstrap」时再用 `--boot-type standard`。
