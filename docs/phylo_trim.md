# 系统发育树修剪整合流程（MAFFT + FastTree + trimAl） | Phylogenetic Tree Trimming Pipeline (MAFFT + FastTree + trimAl)

一句话理解：**把「比对 → 建树」和「修掉比对里的噪声列 → 再建一棵树」串成一条流水线，让你一次拿到修剪前、修剪后两棵树做对比**。

输入一个原始序列 FASTA，自动完成 MAFFT 比对、FastTree 建树（修剪前）、trimAl 修剪、FastTree 再建树（修剪后）。

## 功能概述 | Overview { #overview }

- 一条命令串起 MAFFT + FastTree + trimAl，输出「修剪前」和「修剪后」两棵系统发育树
- trimAl 修剪掉比对里 gap 太多、信息量低的列，让树更可靠
- 支持 trimAl 六种修剪策略（automated1 / gappyout / strict / strictplus / gt / cons）
- 可 `--skip-trimal` 只出修剪前的一棵树
- 断点续传：前流程、trimAl、后树三步各自按输出文件存在性跳过
- 复用 mafft_fasttree 与 trimal 两个模块（只 import，不改动它们）

## 快速开始 | Quick Start { #quick-start }

`@bash
biopytools phylo-trim -i seqs.fa -o out/
`@

最小输入：一个原始序列 FASTA（蛋白或核酸），`-o` 指定输出目录（默认 `./phylo_trim_output`）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 多序列比对(MSA) | 把同源序列「对齐」，同一位点排同一列 |
| gap | 对齐里用来补空的 `-`，代表这条序列在这里缺位 |
| 修剪(trim) | 删掉比对里「gap 太多、没信息」的列，让剩下的列更可信 |
| 修剪前/后树 | 用「没修剪」和「修剪过」的比对各建一棵树，对比差异 |
| 保守度(cons) | 一列上「大家有多一致」，越高越可信 |

## 输入 | Input { #input }

一个原始序列 FASTA（蛋白或核酸均可，序列类型自动检测，也可 `--seq-type` 指定）。

`@text
>seq1
MSTNPKPQRKTKRN
>seq2
MSTNPKPQRKTKRN
>seq3
MSTNPKPQRKTKRS
`@

## 参数说明 | Parameters { #parameters }

### 输入输出与序列类型 | I/O & seq type

**通俗理解|In plain words:** `-i` 给原始 FASTA，`-o` 给输出目录。`--seq-type` 一般不用写（自动检测）。`--sample-name` 管「输出文件叫什么名字」，默认取输入文件去扩展名的名字，**建议保持默认**（见 FAQ）。

### trimAl 修剪 | trimAl trimming

**通俗理解|In plain words:** `--trimal-method` 选修剪策略——默认 `automated1`（全自动、兼顾信息量与噪声），**一般不用动**；`gappyout` 只按 gap 比例删、`strict`/`strictplus` 更激进。`gt-threshold`（默认 0.9）和 `cons-threshold`（默认 80）只在选 `gt`（按 gap 阈值删列）或 `cons`（按保守度删列）时才生效。`--trimal-format` 管修剪后比对写成什么格式，默认 keep（沿用输入格式），一般不用动。`--skip-trimal` 跳过修剪、只出修剪前树。

### MAFFT / FastTree 与性能 | MAFFT / FastTree & performance

**通俗理解|In plain words:** `-t`（默认 88）管 MAFFT 线程；`--mafft-params`（默认 `--auto`）、`--fasttree-params`（默认空）是透传给底层软件的参数，**一般不用动**；`--mafft-path`/`--fasttree-path` 指定软件路径，默认自动解析功能域环境，缺失时回退 PATH 查找。

## 分析流程 | Pipeline { #pipeline }

`@text
步骤1: 检测序列类型 -> 建目录
   |
   v
步骤2: 前流程 = MAFFT 比对 -> FastTree 建树（修剪前树）
   |     -> 01_mafft_fasttree/{sample}.mafft.fa / .fasttree.nwk
   v
步骤3: trimAl 修剪（可 --skip-trimal 跳过）
   |     -> 02_trimal/01_trimal/{sample}.trimmed{ext}
   v
步骤4: 后流程 = FastTree 在修剪后比对再建树（修剪后树）
   |     -> 03_fasttree_trimmed/{sample}.fasttree.trimmed.nwk
   v
步骤5: 写 software_versions.yml
`@

## 输出 | Output { #output }

`@text
output_dir/
|-- 00_pipeline_info/
|   +-- software_versions.yml           # 软件版本与参数记录
|-- 01_mafft_fasttree/                  # 前流程（修剪前）
|   |-- {sample}.cleaned.fa             # 清理后的序列
|   |-- {sample}.mafft.fa               # MAFFT 比对
|   |-- {sample}.fasttree.nwk           # 修剪前的树
|   +-- {sample}.id_mapping.txt         # 新旧 ID 对照
|-- 02_trimal/
|   +-- 01_trimal/
|       +-- {sample}.trimmed{ext}       # trimAl 修剪后的比对
|-- 03_fasttree_trimmed/
|   +-- {sample}.fasttree.trimmed.nwk   # 修剪后的树
+-- 99_logs/
    +-- phylo_trim.log                  # 运行日志
`@

## 结果解读 | Interpreting Results { #interpreting }

- **`01_mafft_fasttree/{sample}.fasttree.nwk`**：修剪前的树（用未修剪的完整比对建的）。
- **`03_fasttree_trimmed/{sample}.fasttree.trimmed.nwk`**：修剪后的树（用 trimAl 清理过的比对建的），一般更可靠。
- **`02_trimal/01_trimal/{sample}.trimmed{ext}`**：修剪后的比对，可直接看删掉了多少列（比对长度变短）。
- **两棵树对比**：若两棵树拓扑一致，说明结果稳健、噪声列没影响结论；若差异大，说明那些 gap 多的列在干扰信号，**应以修剪后树为准**。
- **`00_pipeline_info/software_versions.yml`**：记录了 mafft/fasttree/trimal 的版本和关键参数，写论文 Methods 时从这里抄。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规流程**：只给 `-i -o`，默认 automated1 + 出两棵树，直接跑。
- **数据 gap 特别多**：可试 `--trimal-method gappyout` 或 `strict`（删得更狠）。
- **想按自己的 gap 阈值删列**：`--trimal-method gt --gt-threshold 0.8`（0-1 之间）。
- **想按保守度删列**：`--trimal-method cons --cons-threshold 70`（0-100）。
- **只要修剪前那棵树、跳过修剪**：加 `--skip-trimal`。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入原始序列 FASTA｜Input raw sequence FASTA |
| `--output-dir, -o` | `./phylo_trim_output` | Path | 输出目录｜Output directory |
| `--skip-trimal` | — |  | 跳过 trimal,只出 trimal 前树｜Skip trimal, before-tree only |
| `--trimal-method` | `automated1` | automated1/gappyout/strict/strictplus/gt/cons | trimal 方法｜trimal method |
| `--gt-threshold` | `0.9` | float | trimal gap 阈值[0,1]｜trimal gap threshold [0,1] |
| `--cons-threshold` | `80` | int | trimal 保守度[0,100]｜trimal conservation [0,100] |
| `--trimal-format` | `keep` | keep/fasta/phylip/phylip_paml/clustal/nexus/nbrf/mega | trimal 输出格式｜trimal output format (keep=沿用输入｜input format) |
| `--seq-type` | — | protein/nucleotide | 序列类型(默认自动检测)｜Sequence type (auto-detect by default) |
| `--threads, -t` | `88` | int | MAFFT 线程数｜MAFFT threads |
| `--mafft-params` | `--auto` |  | MAFFT 额外参数｜Additional MAFFT parameters |
| `--fasttree-params` | `` |  | FastTree 额外参数｜Additional FastTree parameters |
| `--mafft-path` | `mafft` |  | MAFFT 路径｜MAFFT path |
| `--fasttree-path` | `fasttree` |  | FastTree 路径｜FastTree path |
| `--sample-name` | — |  | 输出文件名前缀(默认输入 basename)｜Output filename prefix (default: input basename) |
| `--log-file` | — | Path | 日志文件路径｜Log file path |
| `--verbose, -v` | — |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入原始序列 FASTA｜Input raw sequence FASTA |
| `-o, --output-dir` | `./phylo_trim_output` |  | 输出目录｜Output directory (default: ./phylo_trim_output) |
| `--skip-trimal` | — | store_true | 跳过 trimal,只出 trimal 前树｜Skip trimal, before-tree only |
| `--trimal-method` | `automated1` |  | trimal 方法｜trimal method (default: automated1) |
| `--gt-threshold` | `0.9` | float | trimal gap 阈值[0,1]｜trimal gap threshold [0,1] (default: 0.9) |
| `--cons-threshold` | `80` | int | trimal 保守度[0,100]｜trimal conservation [0,100] (default: 80) |
| `--trimal-format` | `keep` |  | trimal 输出格式｜trimal output format (default: keep) |
| `--seq-type` | — | protein/nucleotide | 序列类型(不指定自动检测)｜Sequence type (auto-detect if not specified) |
| `-t, --threads` | `88` | int | MAFFT 线程数｜MAFFT threads (default: 88) |
| `--mafft-params` | `--auto` |  | MAFFT 额外参数｜Additional MAFFT parameters (default: --auto) |
| `--fasttree-params` | `` |  | FastTree 额外参数｜Additional FastTree parameters |
| `--mafft-path` | — |  | MAFFT 路径(默认域环境自动解析)｜MAFFT path (default: auto domain env) |
| `--fasttree-path` | — |  | FastTree 路径(默认域环境自动解析)｜FastTree path (default: auto domain env) |
| `--sample-name` | — |  | 输出文件名前缀｜Output filename prefix (default: 输入 basename｜input basename) |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `-v, --verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- MAFFT（自动解析 phylo 域环境并经 conda run 调用；`--mafft-path` 或环境变量 MAFFT_PATH 覆盖；域环境缺失时回退 PATH 直接调用）
- FastTree（无对应功能域环境，默认 `fasttree`；`--fasttree-path` 或环境变量 FASTTREE_PATH 覆盖，回退 PATH 直接调用）
- trimAl（默认路径 `~/miniforge3/envs/phylo/bin/trimal`，自动解析 phylo 域环境并经 conda run 调用，可用 `TRIMAL_PATH` 覆盖）
- Python 3（biopytools 运行环境）

## 常见问题 | FAQ { #faq }

**Q1：换参数重跑，结果没变？**
本模块**有断点续传**：前树、trimAl、后树各自按输出文件存在性跳过。换 `--trimal-method` 或阈值重跑前，先删掉对应的旧产物（如 `02_trimal/` 和 `03_fasttree_trimmed/` 下的文件），否则会复用旧结果。

**Q2：`--sample-name` 能随便改吗？**
建议保持默认（取输入文件 basename）。内部 MAFFT/FastTree 的中间文件按输入 basename 命名，而 phylo-trim 自己的前/后树和 trimAl 文件按 `sample-name` 命名；两者不一致会导致「跳步判断」和「找中间文件」出错。要改请确保与输入 basename 一致。

**Q3：`--skip-trimal` 之后还有后树吗？**
没有了，只会出 `01_mafft_fasttree/` 下的修剪前树，`02_trimal/` 和 `03_fasttree_trimmed/` 目录不会生成。

**Q4：`gt` 和 `cons` 方法为什么要额外给阈值？**
`gt` 按「gap 比例超过 `--gt-threshold` 的列删掉」；`cons` 按「保守度低于 `--cons-threshold` 的列删掉」。默认 `automated1` 则自动判断、无需阈值，所以大多数情况用默认就好。

**Q5：跑一半中断了，重跑会接着来吗？**
会。已完成的步骤（前树/trimAl/后树）会按输出文件存在性跳过，从断掉的那一步继续。
