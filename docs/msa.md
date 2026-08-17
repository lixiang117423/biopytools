# 多序列比对 | Multiple Sequence Alignment (MSA)

一句话理解：**把一堆同源序列「对齐到同一把尺子上」，让相同位置排在同一列，是建树、找保守区、分析变异前的必要一步**。

输入一个未比对的 FASTA，可选 MAFFT / Clustal Omega / MUSCLE / T-Coffee 四种方法，输出比对文件 + 一份统计报告。

## 功能概述 | Overview { #overview }

- 支持四种主流比对工具：MAFFT（默认）、Clustal Omega、MUSCLE、T-Coffee
- 输出格式可选：FASTA（默认）、Clustal、PHYLIP、NEXUS
- 自动生成比对统计报告（序列数、比对长度、平均一致性）
- 断点续传：无（每次运行重新比对）

## 快速开始 | Quick Start { #quick-start }

`@bash
biopytools msa -i sequences.fa -o alignment
`@

最小输入：一个 FASTA 文件，`-o` 指定**输出文件前缀**（不是目录，见 FAQ）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 多序列比对(MSA) | 把多条序列「对齐」，同一位点排同一列，中间用 `-`（gap）补空 |
| 一致性(identity) | 某一列上「大家有多像」的平均值，越接近 100% 越保守 |
| 保守区 | 几乎所有序列都一样的区域，通常功能重要 |
| 比对格式 | 同一份比对结果的不同「写法」（FASTA/Clustal/PHYLIP/NEXUS），给不同下游软件用 |
| 迭代次数 | 比对算法「反复打磨」的轮数，越多越精细也越慢 |

## 输入 | Input { #input }

一个未比对的 FASTA 文件（蛋白或核酸均可，比对软件自动处理）。

`@text
>seq1
MSTNPKPQRKTKRN
>seq2
MSTNPKPQRKTKRN
>seq3
MSTNPKPQRKTKRS
`@

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 输入 FASTA + 输出前缀。注意 `-o` 是**前缀**不是目录：最终会生成 `<前缀>.<格式>` 比对文件、`<前缀>_stats.txt` 统计、`<前缀>.log` 日志三个文件。

### 比对方法与输出格式 | Method & format

**通俗理解|In plain words:** `-m` 选用哪个比对软件——**默认 MAFFT 即可**，速度快、质量好，绝大多数场景不用换；只有需要和其他研究对齐工具口径时才换。`-f` 选输出格式，取决于下游软件要什么（建树常用 FASTA，某些软件要 PHYLIP/NEXUS），默认 FASTA 最通用。

### 各工具的迭代参数 | Per-tool iterations

**通俗理解|In plain words:** 这些参数管「对应工具打磨多少遍」。`--mafft-maxiterate`（默认 1000）只在自动策略下生效，**一般不用动**；`--clustalo-iterations` 默认 0（不迭代）、`--muscle-maxiters` 默认 16，也基本不用调——调大只会更慢，质量提升有限。

### MAFFT 策略 | MAFFT strategy

**通俗理解|In plain words:** 只有用 MAFFT 时相关。默认 `auto` 让 MAFFT 自己挑策略，**一般不用动**；`linsi`/`ginsi`/`einsi` 分别对应局部/全局/结构敏感的精细比对，序列少且要最高精度时才用。

### 工具路径 | Tool paths

**通俗理解|In plain words:** 指定各比对软件的可执行文件路径，默认在系统 PATH 里找，**一般不用动**。

## 分析流程 | Pipeline { #pipeline }

`@text
步骤1: 检查所选比对工具是否可用
   |
   v
步骤2: 统计输入序列条数
   |
   v
步骤3: 调用所选工具执行比对 -> <前缀>.<格式>
   |
   v
步骤4: 计算比对统计 -> <前缀>_stats.txt
`@

## 输出 | Output { #output }

`@text
<前缀>.fasta          # 比对结果（格式由 -f 决定，默认 fasta）
<前缀>_stats.txt      # 比对统计报告（序列数/比对长度/平均一致性）
<前缀>.log            # 运行日志
`@

## 结果解读 | Interpreting Results { #interpreting }

- **`<前缀>.<格式>`**：比对结果本体，可直接用于下游建树（如 `biopytools iqtree`、`biopytools mafft-fasttree`）或可视化。序列现在等长、同一位点对齐。
- **`<前缀>_stats.txt`**：看三个数——序列数量、比对长度、**平均一致性**。一致性高（如 >90%）说明序列保守；低（<50%）说明差异大，比对可能含较多 gap，建树前或需修剪。
- **好坏判据**：用可视化工具（如 `biopytools msaviz`）检查比对，gap 分布合理、保守区清晰即可；若 gap 满天飞、明显错位，考虑换更精细的策略（如 `--mafft-strategy linsi`）或先修序列。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规比对**：只给 `-i -o`，默认 MAFFT + auto + FASTA 即可。
- **要建树**：保持 `-f fasta`（大多数建树工具接受 FASTA 比对）。
- **下游要 PHYLIP**：加 `-f phylip`。
- **序列少、要最高精度**：MAFFT 下加 `--mafft-strategy linsi`。
- **需要和文献口径一致**：按对方工具选 `-m`。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入序列文件(FASTA格式)｜Input sequence file (FASTA format) |
| `--output, -o` | 必填 |  | 输出文件前缀｜Output file prefix |
| `--method, -m` | `mafft` | mafft/clustalo/muscle/t_coffee | 比对方法｜Alignment method |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--format, -f` | `fasta` | fasta/clustal/phylip/nexus | 输出格式｜Output format |
| `--mafft-strategy` | `auto` | auto/linsi/ginsi/einsi/fftns/fftnsi | MAFFT比对策略｜MAFFT alignment strategy |
| `--mafft-maxiterate` | `1000` | int | MAFFT最大迭代次数｜MAFFT max iterations |
| `--clustalo-iterations` | `0` | int | Clustal Omega迭代次数｜Clustal Omega iterations |
| `--muscle-maxiters` | `16` | int | MUSCLE最大迭代次数｜MUSCLE max iterations |
| `--mafft-path` | `mafft` |  | MAFFT程序路径｜MAFFT program path |
| `--clustalo-path` | `clustalo` |  | Clustal Omega程序路径｜Clustal Omega program path |
| `--muscle-path` | `muscle` |  | MUSCLE程序路径｜MUSCLE program path |
| `--tcoffee-path` | `t_coffee` |  | T-Coffee程序路径｜T-Coffee program path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入序列文件 (FASTA格式)｜Input sequence file (FASTA format) |
| `-o, --output` | 必填 |  | 输出文件前缀｜Output file prefix |
| `-m, --method` | `mafft` | mafft/clustalo/muscle/t_coffee | 比对方法｜Alignment method |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `-f, --format` | `fasta` | fasta/clustal/phylip/nexus | 输出格式｜Output format |
| `--mafft-strategy` | `auto` | auto/linsi/ginsi/einsi/fftns/fftnsi | MAFFT策略｜MAFFT strategy |
| `--mafft-maxiterate` | `1000` | int | MAFFT最大迭代次数｜MAFFT max iterations |
| `--clustalo-iterations` | `0` | int | Clustal Omega迭代次数｜Clustal Omega iterations |
| `--muscle-maxiters` | `16` | int | MUSCLE最大迭代次数｜MUSCLE max iterations |
| `--mafft-path` | `mafft` |  | MAFFT路径｜MAFFT path |
| `--clustalo-path` | `clustalo` |  | Clustal Omega路径｜Clustal Omega path |
| `--muscle-path` | `muscle` |  | MUSCLE路径｜MUSCLE path |
| `--tcoffee-path` | `t_coffee` |  | T-Coffee路径｜T-Coffee path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- MAFFT（默认 `mafft`）、Clustal Omega（`clustalo`）、MUSCLE（`muscle`）、T-Coffee（`t_coffee`）——按所选方法至少需要一个，均在 PATH 中查找
- Python 3 + BioPython（用于生成统计报告；缺失时统计步骤会跳过但不影响比对结果）

## 常见问题 | FAQ { #faq }

**Q1：`-o` 是目录还是文件名？**
是**输出前缀**（不是目录）。比如 `-o alignment` 会生成 `alignment.fasta`、`alignment_stats.txt`、`alignment.log`。想放某目录就写 `-o /path/to/alignment`。

**Q2：选 MAFFT 还是 MUSCLE/Clustal Omega？**
默认 MAFFT 即可，速度快、质量好。MUSCLE 注意其 v5 的命令行格式与 v3/v4 不同，本模块按 v3/v4 语法调用，若你的 MUSCLE 是 v5 可能不兼容；有特殊需求再换。

**Q3：统计报告没生成？**
统计步骤依赖 BioPython，没装的话会跳过（比对结果不受影响）。需要统计就装 BioPython。

**Q4：重跑会复用旧结果吗？**
不会，每次从头比对，且日志 `<前缀>.log` 会被覆盖。

**Q5：比对结果里有好多 `-`（gap）正常吗？**
正常，gap 代表「这条序列在这里缺位」。但如果 gap 特别多、分布很乱，说明序列之间差异大或方向不对（混了反向互补序列），建议先检查输入。
