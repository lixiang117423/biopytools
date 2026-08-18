# 系统发育树构建（MAFFT + FastTree） | Phylogenetic Tree Builder (MAFFT + FastTree)

一句话理解：**把一堆未比对的序列一键「对齐」再用快速算法「画成一棵进化树」，是快速建树的标准流水线**。

输入一个 FASTA（蛋白或核酸），自动完成序列清理、MAFFT 多序列比对、FastTree 建树三步，输出比对文件和 Newick 树文件。

## 功能概述 | Overview { #overview }

- MAFFT 做多序列比对 + FastTree 构建最大似然树，两步串联一键完成
- 自动检测序列类型（核酸/蛋白），也支持手动指定
- 建树前自动清理序列：ID 里特殊字符换成下划线、重复 ID 加后缀，并生成新旧 ID 对照表
- FastTree 速度极快，适合几百到上万条序列的大规模建树
- 断点续传：无（每次运行从头跑，日志文件启动时会被覆盖）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools mafft-fasttree -i sequences.fa -o phylo_results
```

最小输入：一个 FASTA 文件（蛋白或核酸序列），`-o` 指定输出目录（默认 `./phylo_output`）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| FASTA | 存序列的标准文本格式，每行 `>` 开头是名字，下面是序列 |
| 多序列比对(MSA) | 把同源序列「对齐」到同一把尺子上，让同一位点排在同一列 |
| Newick 树 | 进化树的标准文本格式，一串括号，树软件都能读 |
| 核酸/蛋白序列 | DNA/RNA 是核酸（A/T/C/G），蛋白是氨基酸（20 种字母） |
| 枝长 | 树上每个分支的长度，代表该分支累积了多少替换 |

## 输入 | Input { #input }

一个 FASTA 文件，蛋白或核酸均可（自动检测）。序列名建议用字母数字和下划线，避免冒号、竖线、空格、括号等特殊字符（程序会自动替换，但影响下游对 ID 的识别）。

```text
>seq1
MSTNPKPQRKTKRN
>seq2
MSTNPKPQRKTKRN
>seq3
MSTNPKPQRKTKRS
```

## 参数说明 | Parameters { #parameters }

### 输入输出与序列类型 | I/O & seq type

**通俗理解|In plain words:** `--seq-type` 管「序列是核酸还是蛋白」。**一般不用写**，程序会自动检测；只有自动检测判错（比如某些氨基酸组合看起来像核酸）时才手动指定。`-o` 指定输出目录。

### 线程 | Threads

**通俗理解|In plain words:** 管并行加速（主要加速 MAFFT）。默认 12，一般不用动。

### MAFFT 与 FastTree 参数 | MAFFT & FastTree params

**通俗理解|In plain words:** 这两项是把额外参数原样「塞给」底层软件。`--mafft-params` 默认 `--auto`（让 MAFFT 自己挑策略），**一般不用动**；`--fasttree-params` 默认空（蛋白用默认模型，核酸自动加 `-nt`）。只有你对 MAFFT/FastTree 很熟、要改比对精度或模型时才填。

### 工具路径 | Tool paths

**通俗理解|In plain words:** 指定 MAFFT、FastTree 的可执行文件路径。默认 `mafft`/`fasttree`（即在系统 PATH 里找），一般不用动；只有软件装在特殊位置时才写全路径。

## 分析流程 | Pipeline { #pipeline }

```text
步骤1: 检查依赖（MAFFT、FastTree 是否可用）
   |
   v
步骤2: 检测序列类型（蛋白 / 核酸）
   |
   v
步骤3: 清理序列 + 处理 ID（特殊字符换下划线、重复 ID 加后缀）
   |     -> {base}.cleaned.fa 与 {base}.id_mapping.txt
   v
步骤4: MAFFT 多序列比对 -> {base}.mafft.fa
   |
   v
步骤5: FastTree 建树 -> {base}.fasttree.nwk
```

## 输出 | Output { #output }

所有文件都写在输出目录里，`{base}` 是输入文件去扩展名后的名字：

```text
output_dir/
|-- {base}.cleaned.fa        # 清理后的序列（ID 已规范化）
|-- {base}.mafft.fa          # MAFFT 比对结果（FASTA）
|-- {base}.fasttree.nwk      # FastTree 系统发育树（Newick 格式）
|-- {base}.id_mapping.txt    # 新旧 ID 对照表（含 Modified/Unchanged 状态）
+-- mafft_fasttree.log       # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting }

- **`{base}.fasttree.nwk`**：最终的系统发育树（Newick 文本）。用 iTOL、FigTree、ggtree 等打开即可看拓扑与枝长。**注意**：FastTree 默认不输出 bootstrap 支持值，只有拓扑和枝长；需要支持值请改用 `biopytools iqtree`。
- **`{base}.mafft.fa`**：比对结果，可留作下游分析（如再修剪、再建树）。
- **`{base}.cleaned.fa`** 与 **`{base}.id_mapping.txt`**：若建树后发现有序列名对不上原始样本，查 id_mapping.txt 就能还原「新 ID 对应原来的哪个 ID」。
- **好坏判据**：先看树上亲缘近的物种是否聚在一起（符合生物学预期）；若拓扑异常，多半是序列 ID 清理、或输入里混了反向/污染序列。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规快速建树**：只给 `-i -o`，全默认即可。
- **序列很多（成千上万条）**：这个组合就是为大规模设计的，默认参数直接跑。
- **蛋白序列被误判成核酸**：加 `--seq-type protein` 手动纠正。
- **想要更高的比对精度**：`--mafft-params "--localpair --maxiterate 1000"`（或 `--globalpair`）。
- **需要分支支持值**：本工具不提供，换 `biopytools iqtree` 或 `biopytools raxml`。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入序列文件(FASTA格式)｜Input sequence file (FASTA format) |
| `--output, -o` | `./phylo_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--seq-type` | — | protein/nucleotide | 序列类型(未指定时自动检测)｜Sequence type (auto-detect if not specified) |
| `--mafft-params` | `--auto` |  | MAFFT额外参数｜Additional MAFFT parameters |
| `--fasttree-params` | `` |  | FastTree额外参数｜Additional FastTree parameters |
| `--mafft-path` | `mafft` |  | MAFFT软件路径｜MAFFT software path |
| `--fasttree-path` | `fasttree` |  | FastTree软件路径｜FastTree software path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入序列文件 (FASTA格式)｜Input sequence file (FASTA format) |
| `-o, --output` | `./phylo_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `--seq-type` | — | protein/nucleotide | 序列类型 (不指定则自动检测)｜Sequence type (auto-detect if not specified) |
| `--mafft-params` | `--auto` |  | MAFFT额外参数｜Additional MAFFT parameters |
| `--fasttree-params` | `` |  | FastTree额外参数｜Additional FastTree parameters |
| `--mafft-path` | — |  | MAFFT软件路径(默认域环境自动解析)｜MAFFT software path (default: auto domain env) |
| `--fasttree-path` | — |  | FastTree软件路径(默认域环境自动解析)｜FastTree software path (default: auto domain env) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- MAFFT（默认 `mafft`，在 PATH 中查找；可用 `--mafft-path` 指定）
- FastTree（默认 `fasttree`，在 PATH 中查找；可用 `--fasttree-path` 指定）
- Python 3 + BioPython（biopytools 运行环境）

## 常见问题 | FAQ { #faq }

**Q1：我的序列名带冒号/空格，会不会出错？**
程序会把 ID 里的 `:`、`|`、空格、`;`、`,`、括号等换成下划线，重复 ID 自动加 `_1`、`_2` 后缀。对照关系都写在 `{base}.id_mapping.txt` 里，不会丢信息。但为了下游省心，最好输入前就把序列名规范好。

**Q2：树上的名字和我的样本名对不上？**
查 `{base}.id_mapping.txt`：里面一行一个 `新ID 原始ID 状态`，Modified 的就是被改动过的，用原始ID 去对应你的样本。

**Q3：能输出带 bootstrap 的树吗？**
不能，FastTree 默认只给拓扑和枝长。需要支持值请用 `biopytools iqtree`。

**Q4：蛋白序列被当成核酸建树了怎么办？**
用 `--seq-type protein` 手动指定即可（或检查序列是否只含 A/T/C/G 字母导致误判）。

**Q5：重跑会复用旧结果吗？**
不会，本模块没有断点续传，每次从头跑，且日志文件 `mafft_fasttree.log` 会被覆盖。
