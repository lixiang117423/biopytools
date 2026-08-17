# wgsim 测序数据模拟 | wgsim Sequencing Simulation

一句话理解：**给一条参考基因组序列，人工「造」出符合设定错误率、读长、插入片段大小的 Illumina 双端测序 reads**，用来测试下游流程或补足缺少的真实测序数据。

## 功能概述 | Overview

- 基于 wgsim 从基因组 FASTA 模拟 Illumina 双端 reads
- 支持单个基因组文件或整个目录批量模拟
- 可控读长、reads 数、错误率、突变率、插入片段距离、随机种子
- 自动修复质量行（替换为 Q40 字符），避免下游 QC 工具误删全部 reads
- 输出 gzip 压缩的 `.fq.gz`，断点续传（已存在则跳过，见 FAQ）

## 快速开始 | Quick Start

```bash
biopytools wgsim -i genome.fna -o output_dir
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 双端测序 | 从一条 DNA 片段的两端各测一次，得到 `_1` 和 `_2` 两个文件，中间隔一段没测到的「插入片段」 |
| 读长 | 每次测序读多少碱基，如 150bp（bp = 碱基对，即字母个数） |
| 插入片段距离 | 两端读段之间那段没测到的 DNA 大概多长；outer 是两端最外侧距离，inner 是两端内侧间距 |
| 错误率 | 模拟测序时故意引入的碱基错误比例，越接近真实测序的噪声水平越好 |
| 突变率 | 模拟 reads 相对参考基因组的变异比例，用来模拟个体差异 |
| 随机种子 | 随机数生成器的「起点」，同一种子每次结果一模一样，方便复现 |
| Phred 质量 | 每个碱基的「可信度打分」，Q40 表示很高的可信度 |

## 输入 | Input

- 基因组 FASTA 文件（`.fna` / `.fa` / `.fasta`），或含这些文件的目录（批量处理）。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 告诉工具「用哪个基因组、写到哪」。

- `-i, --input`：输入基因组文件或目录。
- `-o, --output-dir`：输出目录。

### 测序量 | Read volume

**通俗理解|In plain words:** 决定「造多少数据」。reads 数越大、读长越长，数据越多、覆盖度越高，但也越慢、文件越大。覆盖度约等于 `reads数 × 读长 × 2 ÷ 基因组大小`。

- `-N, --num-reads`：模拟 reads 数量（默认 50000000，即 5000 万）。
- `-1, --read-length`：reads 长度（默认 150bp）。

### 模拟质量 | Error & mutation

**通俗理解|In plain words:** 决定造出来的 reads「像不像真实测序」。错误率模拟测序噪声，突变率模拟个体差异。**一般用默认值即可。**

- `-e, --error-rate`：测序错误率（默认 0.020）。
- `-r, --mutation-rate`：突变率（默认 0.001）。

### 插入片段 | Insert size

**通俗理解|In plain words:** 决定双端 reads 之间的间隔。默认值对应常见 500bp 左右插入片段文库；要模拟特定文库时再改。

- `-d, --outer-distance`：外部距离（默认 500）。
- `-D, --inner-distance`：内部距离（默认 0）。

### 随机性 | Randomness

**通俗理解|In plain words:** 控制结果能否复现。要复现结果就设一个固定正整数；默认 0 表示用时间做种子，每次跑结果都不同。

- `-s, --seed`：随机种子（默认 0）。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先跑 wgsim 生成原始 FASTQ，再修复质量行，最后压缩成 `.fq.gz`。

1. wgsim 从基因组模拟出 `{base}_1.fq` / `{base}_2.fq` 双端 reads
2. 修复质量行：把 wgsim 恒定的低质量字符替换为 `I`（Q40）
3. gzip 压缩为 `{base}_1.fq.gz` / `{base}_2.fq.gz`

## 输出 | Output

```text
output_dir/
├── genome_1.fq.gz          # 双端 Read1（压缩）
├── genome_2.fq.gz          # 双端 Read2（压缩）
└── wgsim.log               # 运行日志
```

- 文件名以基因组文件名（去扩展名）为前缀，后缀 `_1.fq.gz` / `_2.fq.gz`。
- 中间未压缩的 `.fq` 在压缩成功后保留为 `.fq.gz`（原 `.fq` 由 gzip 就地处理）。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 拿到 `.fq.gz` 后可用 fastp/fastqc 等再验证一遍；主要确认 reads 数、读长、配对是否正常。

- 文件大小应与预期 reads 数 × 读长 × 2 大致匹配（gzip 后会更小）。
- 质量行统一为 Q40，属预期行为（工具刻意修复，见 FAQ），不代表真实质量分布。
- 若某个基因组模拟失败，日志会记录，且对应输出文件会被清理。

## 参数选择建议 | Parameter Guidance

- **覆盖度**：先算目标覆盖度，再反推 `-N`。公式 `N = 覆盖度 × 基因组大小 ÷ (读长 × 2)`。例如 500Mb 基因组、150bp 读长、要 30x，则 `N = 30 × 5e8 ÷ 300 = 5e7`（约 5000 万，正是默认值）。
- **`-s`**：需要可复现时设固定正整数（如 42），否则默认 0 每次结果不同。
- **`-e` / `-r`**：默认值贴近常规 Illumina 数据；想模拟更高/更低噪声再调。
- **`-d` / `-D`**：模拟特定插入片段文库（如 350bp 或 1kb）时按文库参数调整，否则默认。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入基因组文件或目录｜Input genome file or directory |
| `-o, --output-dir` | 必填 | Path | 输出目录｜Output directory |
| `-N, --num-reads` | `50000000` | int | 模拟reads数量｜Number of reads to simulate |
| `-1, --read-length` | `150` | int | reads长度｜Read length |
| `-s, --seed` | `0` | int | 随机种子｜Random seed |
| `-e, --error-rate` | `0.02` | float | 测序错误率｜Sequencing error rate |
| `-r, --mutation-rate` | `0.001` | float | 突变率｜Mutation rate |
| `-d, --outer-distance` | `500` | int | 外部距离｜Outer distance |
| `-D, --inner-distance` | `0` | int | 内部距离｜Inner distance |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | [FILE/DIR] 输入基因组文件或目录｜Input genome file or directory |
| `-o, --output-dir` | 必填 |  | [DIR] 输出目录｜Output directory |
| `-N, --num-reads` | `50000000` | int | [INT] 模拟reads数量，默认50000000｜Number of reads to simulate, default 50000000 |
| `-1, --read-length` | `150` | int | [INT] reads长度，默认150｜Read length, default 150 |
| `-s, --seed` | `0` | int | [INT] 随机种子，默认0｜Random seed, default 0 |
| `-e, --error-rate` | `0.02` | float | [FLOAT] 测序错误率，默认0.020｜Sequencing error rate, default 0.020 |
| `-r, --mutation-rate` | `0.001` | float | [FLOAT] 突变率，默认0.001｜Mutation rate, default 0.001 |
| `-d, --outer-distance` | `500` | int | [INT] 外部距离，默认500｜Outer distance, default 500 |
| `-D, --inner-distance` | `0` | int | [INT] 内部距离，默认0｜Inner distance, default 0 |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- wgsim（默认路径 `~/miniforge3/envs/align/bin/wgsim`，可用环境变量 `WGSIM_PATH` 或配置文件覆盖）
- gzip（系统自带，用于压缩输出）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持。若某基因组的 `{base}_1.fq.gz` 和 `{base}_2.fq.gz` 都已存在，则跳过该样本不重算。

**Q2：为什么质量行全是 `I`（Q40）？**
wgsim 会按错误率生成恒定的低质量值，下游 fastp 等工具可能因质量阈值过高而丢弃全部 reads。工具因此把质量行统一替换为 `I`（Q40），保证 reads 能被正常保留；这属预期行为。

**Q3：默认 5000 万 reads 会不会太大？**
对大基因组是合理覆盖度，但很耗时间和磁盘。小基因组或快速测试时请把 `-N` 调小（如 100 万）。

**Q4：为什么两次运行结果不一样？**
默认 `--seed 0` 表示用时间做随机种子，每次不同；设一个固定正整数即可复现。

**Q5：输出为什么是 `.fq.gz` 压缩格式？**
模拟结果体积大，默认 gzip 压缩省空间，下游工具（fastp、bwa 等）都能直接读 `.fq.gz`。
