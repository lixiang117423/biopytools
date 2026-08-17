# PanMAN 泛基因组网络构建与分析 | PanMAN Pangenome Network Construction and Analysis

一句话理解：把多个基因组的"泛基因组图"叠加上系统发育树，构建一个既能存序列、又能查变异、还能追溯演化关系的网络数据结构，之后可从中提取各种格式的序列/变异/比对。

## 功能概述 | Overview

- 从 PanGraph JSON / GFA / MSA + Newick 树构建 PanMAN 泛基因组网络
- 内置 generate-pangraph：从 FASTA 直接生成 PanGraph JSON + Newick 树
- 从 PanMAN 提取多种数据：摘要、FASTA、MSA、VCF、GFA、Newick、MAF、氨基酸翻译等
- 高级操作：子网络提取、节点注释、重新扎根、创建网络、打印突变、范围查询
- 三种运行后端：conda(默认)、docker、singularity

## 快速开始 | Quick Start

```bash
biopytools panman build -p input.json -n tree.nwk -o my_panman
```

`-p` 是 PanGraph JSON(也可用 `-g` GFA 或 `-m` MSA)，`-n` 是 Newick 树(必需)，`-o` 是输出前缀。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 泛基因组图(pangenome graph) | 把多个基因组的序列"拼"成一张有分叉的网络图，既能表示共有部分也能表示差异 |
| PanGraph | 一个具体的泛基因组图构建工具，输出 JSON 格式的图 |
| GFA | 图的通用文本格式，描述节点(片段)和边(连接关系) |
| MSA | 多序列比对，把同源序列对齐好排成列 |
| Newick 树 | 系统发育树的文本记法，用括号嵌套表示亲缘关系 |
| VCF | 记录基因组间/个体间"哪里不一样"的标准变异格式 |
| MAF | 多序列比对的对齐块格式 |
| 重新扎根(reroot) | 给树换一个"根"，改变祖先-后代的叙述方向 |
| ACR | 祖先状态重建，推断祖先节点上是什么序列/状态 |
| 后端(backend) | 用哪种方式运行软件：conda 环境 / docker / singularity 容器 |

## 输入 | Input

三个子命令的输入不同：

### build

- 至少一个图输入：`-p` PanGraph JSON、`-g` GFA、`-m` MSA 三选一或组合
- `-n` Newick 树文件(必需)

### extract

- `-i` 输入 PanMAN 文件(必需，由 build 生成)
- 至少一个提取选项(如 `--summary`、`--vcf` 等)

### generate-pangraph

- `-i` 输入 FASTA 文件(含多条序列)

## 参数说明 | Parameters

### 构建输入 | Build inputs

**通俗理解|In plain words:** build 时告诉程序"图"从哪来(`-p`/`-g`/`-m` 至少给一个)，"树"从哪来(`-n` 必需)。三者可组合使用。

### 提取选项 | Extraction options

**通俗理解|In plain words:** extract 时用一堆 `--xxx` 开关勾选"要导出哪些格式"。每个开关对应一种输出文件；可以一次勾多个。VCF/重新扎根需要 `-r` 参考序列名；子网络/注释/建网络需要 `--input-file`。

### 后端与环境 | Backend & environment

**通俗理解|In plain words:** `--backend` 选运行方式：默认 conda(用 `--conda-env` 指定环境名，默认 `panman_v.0.1.4`)；无 conda 时可用 singularity/docker(需给 `--sif-image`/`--singularity-path`)。**一般保持默认 conda 即可。**

### 高级 | Advanced

**通俗理解|In plain words:** `--acr`(祖先状态重建方法 fitch/mppa)、`--tree-id`(VCF 提取指定哪棵树)、`--range-*`(按坐标范围提取序列)等，仅在特定提取任务中需要，一般不用动。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-p, --pangraph` | — |  | PanGraph JSON文件｜PanGraph JSON file path |
| `-g, --gfa` | — |  | GFA文件路径｜GFA file path |
| `-m, --msa` | — |  | MSA文件路径 (FASTA格式)｜MSA file path (FASTA format) |
| `-n, --newick` | 必填 |  | Newick树文件路径｜Newick tree file path |
| `-o, --output-prefix` | `output` |  | 输出文件前缀｜Output file prefix |
| `--output-dir` | `./panman_output` |  | 输出目录｜Output directory |
| `--backend` | `conda` | conda/docker/singularity | 后端选择｜Backend selection |
| `--conda-env` | `panman_v.0.1.4` |  | Conda环境名称｜Conda environment name |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--sif-image` | — |  | PanMAN Singularity SIF镜像路径｜PanMAN Singularity SIF image path |
| `--singularity-path` | — |  | Singularity可执行文件路径｜Singularity executable path |
| `-i, --input-panman` | 必填 |  | PanMAN文件路径｜PanMAN file path |
| `-r, --reference` | — |  | 参考序列名称 (VCF提取/重新扎根需要)｜Reference sequence name (required for VCF/reroot) |
| `--summary` | — |  | 提取摘要统计｜Extract summary statistics |
| `--extract-fasta` | — |  | 提取FASTA序列｜Extract FASTA sequences |
| `--extract-msa` | — |  | 提取MSA比对｜Extract MSA alignment |
| `--vcf` | — |  | 提取VCF变异｜Extract VCF variants |
| `--extract-gfa` | — |  | 提取GFA格式｜Extract GFA format |
| `--extract-newick` | — |  | 提取Newick树｜Extract Newick tree |
| `--extended-newick` | — |  | 提取扩展Newick格式｜Extract extended Newick format |
| `--maf` | — |  | 提取MAF格式｜Extract MAF format |
| `--aa` | — |  | 提取氨基酸翻译｜Extract amino acid translations |
| `--subnet` | — |  | 提取子网络｜Extract subnet (requires --input-file) |
| `--annotate` | — |  | 注释节点｜Annotate nodes (requires --input-file) |
| `--reroot` | — |  | 重新扎根树｜Reroot tree (requires --reference) |
| `--create-network` | — |  | 创建网络｜Create network (requires --input-file) |
| `--print-mutations` | — |  | 打印突变信息｜Print mutations |
| `--acr` | `fitch` | fitch/mppa | ACR方法｜ACR method (fitch/mppa) |
| `--tree-id` | — |  | 树ID (VCF提取可选)｜Tree ID (optional for VCF extraction) |
| `--input-file` | — |  | 输入文件路径 (subnet/annotate/create-network需要)｜Input file path (required for subnet/annotate/create-network) |
| `--range-index` | — |  | 范围查询index参数｜Range query index parameter |
| `--range-start` | — | int | 范围查询起始坐标｜Range query start coordinate |
| `--range-end` | — | int | 范围查询结束坐标｜Range query end coordinate |
| `-i, --fasta` | 必填 |  | 输入FASTA文件路径｜Input FASTA file path |
| `--pangraph-path` | — |  | PanGraph可执行文件路径｜PanGraph executable path |
| `--pangraph-sif` | — |  | PanGraph Singularity SIF镜像路径｜PanGraph Singularity SIF image path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--build` | — | store_true | 构建模式｜Build mode: build PanMAN from alignments |
| `--extract` | — | store_true | 提取模式｜Extract mode: extract data from PanMAN |
| `--generate-pangraph` | — | store_true | PanGraph生成模式｜Generate PanGraph mode: generate PanGraph JSON from FASTA |
| `-i, --fasta` | — |  | 输入FASTA文件路径 (PanGraph生成模式必需)｜Input FASTA file path (required for generate-pangraph mode) |
| `-P, --pangraph` | — |  | PanGraph JSON文件路径｜PanGraph JSON file path |
| `-G, --gfa` | — |  | GFA文件路径｜GFA file path |
| `-M, --msa` | — |  | MSA文件路径 (FASTA格式)｜MSA file path (FASTA format) |
| `-N, --newick` | — |  | Newick树文件路径 (构建模式必需)｜Newick tree file path (required for build mode) |
| `-I, --input-panman` | — |  | PanMAN文件路径 (提取模式必需)｜PanMAN file path (required for extract mode) |
| `-o, --output-prefix` | `output` |  | 输出文件前缀｜Output file prefix |
| `--output-dir` | `./panman_output` |  | 输出目录｜Output directory |
| `-r, --reference` | — |  | 参考序列名称 (VCF提取需要)｜Reference sequence name (required for VCF extraction) |
| `--summary` | — | store_true | 提取摘要统计｜Extract summary statistics |
| `--extract-fasta` | — | store_true | 提取FASTA序列｜Extract FASTA sequences |
| `--extract-msa` | — | store_true | 提取MSA比对｜Extract MSA alignment |
| `--vcf` | — | store_true | 提取VCF变异｜Extract VCF variants |
| `--extract-gfa` | — | store_true | 提取GFA格式｜Extract GFA format |
| `--extract-newick` | — | store_true | 提取Newick树｜Extract Newick tree |
| `--extended-newick` | — | store_true | 提取扩展Newick格式｜Extract extended Newick format |
| `--maf` | — | store_true | 提取MAF格式｜Extract MAF format |
| `--aa` | — | store_true | 提取氨基酸翻译｜Extract amino acid translations |
| `--subnet` | — | store_true | 提取子网络｜Extract subnet |
| `--annotate` | — | store_true | 注释节点｜Annotate nodes |
| `--reroot` | — | store_true | 重新扎根树｜Reroot tree |
| `--create-network` | — | store_true | 创建网络｜Create network |
| `--print-mutations` | — | store_true | 打印突变信息｜Print mutations |
| `--backend` | `conda` | conda/docker/singularity | 后端选择｜Backend selection (default: conda) |
| `--conda-env` | `panman_v.0.1.4` |  | Conda环境名称｜Conda environment name (default: panman_v.0.1.4) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--pangraph-path` | — |  | PanGraph可执行文件路径｜PanGraph executable path |
| `--pangraph-sif` | — |  | PanGraph Singularity SIF镜像路径｜PanGraph Singularity SIF image path |
| `--sif-image` | — |  | Singularity SIF镜像路径｜Singularity SIF image path |
| `--singularity-path` | — |  | Singularity可执行文件路径｜Singularity executable path |
| `--acr` | `fitch` | fitch/mppa | ACR方法｜ACR method (fitch/mppa) |
| `--tree-id` | — |  | 树ID (VCF提取需要)｜Tree ID (required for VCF extraction) |
| `--input-file` | — |  | 输入文件路径 (用于subnet/annotate/create-network)｜Input file path (for subnet/annotate/create-network) |
| `--range-index` | — |  | 范围查询index参数｜Range query index parameter |
| `--range-start` | — | int | 范围查询起始坐标｜Range query start coordinate |
| `--range-end` | — | int | 范围查询结束坐标｜Range query end coordinate |

<!-- END PARAMS:auto -->

## 分析流程 | Pipeline

**通俗理解|In plain words:** 要么"从 FASTA 建图"(generate-pangraph)，要么"拿现成图+树建网络"(build)，建好后"按需导出"(extract)。

```text
build:        PanGraph/GFA/MSA + Newick 树 → panmanUtils → {prefix}.panman
extract:      {prefix}.panman → panmanUtils → 各格式输出文件
generate-pangraph: FASTA → pangraph build → {prefix}.json + {prefix}.nwk
```

## 输出 | Output

### build 输出

```text
panman_output/
└── {prefix}.panman      # PanMAN 网络文件(后续 extract 的输入)
```

### extract 输出(按勾选的选项生成)

```text
panman_output/
├── {prefix}_summary.txt          # --summary 摘要统计
├── {prefix}.fa                   # --extract-fasta 序列
├── {prefix}_msa.fa               # --extract-msa 多序列比对
├── {prefix}.vcf                  # --vcf 变异(需 -r 参考)
├── {prefix}.gfa                  # --extract-gfa 图格式
├── {prefix}.nwk                  # --extract-newick 树
├── {prefix}_extended.nwk         # --extended-newick 扩展树
├── {prefix}.maf                  # --maf 对齐块
├── {prefix}_aa.tsv               # --aa 氨基酸翻译
├── {prefix}_subnet.panman        # --subnet 子网络(需 --input-file)
├── {prefix}_annotated.panman     # --annotate 节点注释(需 --input-file)
├── {prefix}_rerooted.panman      # --reroot 重新扎根(需 -r)
├── {prefix}_network.panman       # --create-network(需 --input-file)
├── {prefix}_mutations.txt        # --print-mutations 突变信息
└── {prefix}_range_x_y.fa         # 范围查询(需 --range-*)
```

### generate-pangraph 输出

```text
panman_output/
├── {prefix}.json         # PanGraph JSON
├── {prefix}.nwk          # 从日志提取的 Newick 树
└── {prefix}_pangraph.log # PanGraph 运行日志
```

## 结果解读 | Interpreting Results

### 1. PanMAN 文件（`{prefix}.panman`）

**通俗理解|In plain words:** 这是"把序列、图、树绑在一起"的最终产物，是后续一切提取的基础。

- 文件本身是二进制/结构化数据，一般用 extract 子命令读取，不直接打开
- 有了它就能反复导出 VCF/FASTA/比对等，无需重跑构建

### 2. VCF（`{prefix}.vcf`）

**通俗理解|In plain words:** 逐位点列出"相对参考序列，各基因组哪里不一样"。

- 需要 `-r` 指定参考序列名；变异坐标相对该参考
- 可交给下游变异分析工具进一步过滤/注释

### 3. 摘要（`{prefix}_summary.txt`）

- 网络规模、节点/边数量、序列长度等统计，快速了解数据概况

### 4. 氨基酸翻译（`{prefix}_aa.tsv`）

- 每条序列的蛋白翻译结果，用于看功能序列是否完整(是否提前终止等)

## 参数选择建议 | Parameter Guidance

- **后端选择**：有 conda 环境(默认 `panman_v.0.1.4`)直接默认；无 conda 用 `--backend singularity` 并给 `--sif-image`/`--singularity-path`
- **build 输入**：有 PanGraph JSON 用 `-p`(信息最全)；只有图用 `-g` GFA；只有比对用 `-m` MSA
- **extract**：先 `--summary` 看概况，再按需导出 VCF/FASTA/MSA 等；一次可勾多个选项
- **VCF/重新扎根**：必须先 `-r` 指定参考序列名，否则报错
- **generate-pangraph**：默认走 singularity(PanGraph 镜像)，也可用 `--pangraph-path` 指定本地可执行文件

## 依赖 | Dependencies

- panmanUtils（conda 环境 `panman_v.0.1.4`，或 docker/singularity 镜像）
- PanGraph（generate-pangraph 需要；本地可执行文件或 `~/software/singularity/pangraph_0.7.3.sif` 镜像）
- Singularity（`~/miniforge3/envs/singularity_v.3.8.7/bin/singularity`，仅 singularity 后端需要）
- Python 3

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
不支持。build/extract/generate-pangraph 每次都是完整执行，没有按步骤跳过的机制。

**Q2：conda 环境不存在怎么办？**
程序会报错并提示 `conda create -n panman_v.0.1.4 panman`。也可改用 `--backend singularity` 或 `--backend docker`。

**Q3：VCF 提取报"需要参考序列"？**
`--vcf` 和 `--reroot` 都需要 `-r/--reference` 指定参考序列名(图中某条序列的名字)。

**Q4：子网络/注释/建网络为什么报错？**
`--subnet`、`--annotate`、`--create-network` 都需要 `--input-file` 提供节点列表/注释/网络清单文件。

**Q5：build 的输入三者都要给吗？**
不需要。`-p`/`-g`/`-m` 至少给一个即可；`-n` Newick 树是必需项。

**Q6：generate-pangraph 用哪种方式跑 PanGraph？**
默认走 singularity 镜像(若 `~/software/singularity/pangraph_0.7.3.sif` 存在)；可用 `--pangraph-path` 指定本地可执行文件，或 `--pangraph-sif` 指定其它镜像。
