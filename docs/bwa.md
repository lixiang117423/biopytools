# bwa - BWA 全基因组比对 | BWA Whole Genome Alignment

一句话理解：**把一批双端 FASTQ 测序数据用 BWA-MEM 比对到参考基因组，自动完成建索引、比对、排序、(可选)去重、统计和覆盖度分析，输出每个样本的 BAM 及质控指标。**

## 功能概述 | Overview

- 封装 BWA-MEM 全流程：基因组索引 → 双端比对 → SAM 转 BAM → 排序 →（可选）标记/移除重复 → 建 BAM 索引
- 每个样本自动生成比对统计（`samtools flagstat` / `stats`）、单碱基覆盖度、滑窗覆盖度
- 自动查找成对 FASTQ（按 `-p` 模式匹配，默认 `_1.clean.fq.gz`）
- 暴露 BWA-MEM 的算法/打分/输出参数（`--bwa-*`），高级用户可调
- 支持 `--resume` 断点续传（跳过已完成样本）；基因组索引自动复用/构建
- 生成整体汇总报告 `alignment_summary.txt`

## 快速开始 | Quick Start

~~~bash
biopytools bwa -g genome.fa -i fastq_dir
~~~

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 比对(mapping) | 把每条读段"贴回"参考基因组上它最可能来自的位置 |
| BWA-MEM | 目前最常用的短读长比对算法，本工具的核心引擎 |
| 双端(paired-end) | 一个片段两头各测一次，得到一对读段，成对信息让比对更准 |
| 种子长度(k) | BWA 先找"短片段种子"再扩展，种子越短越敏感但越慢 |
| 重复标记(markdup) | 把 PCR 产生的完全一样读段标记出来，避免它们干扰后续分析 |
| 覆盖度(深度) | 每个碱基平均被测多少次；滑窗覆盖度是按固定窗口统计的"区域性"深度 |

## 输入 | Input

- 参考基因组 FASTA（`-g`，必需）
- FASTQ 目录（`-i`，必需），内含成对 FASTQ，read1 文件名含 `-p` 指定的模式（默认 `_1.clean.fq.gz`），read2 由 `_1.` 换成 `_2.` 自动推出

~~~text
fastq_dir/
├── sampleA_1.clean.fq.gz
├── sampleA_2.clean.fq.gz
├── sampleB_1.clean.fq.gz
└── sampleB_2.clean.fq.gz
~~~

- 找不到配对文件（read2）的 read1 会被跳过并警告

## 参数说明 | Parameters

### 基本与性能 | Basic & performance

**通俗理解|In plain words:** `-g` 参考基因组、`-i` FASTQ 目录必填；`-p` 是 read1 的匹配模式（不同项目的后缀不同时才改）；`-o` 输出目录；`-t` 线程数（默认 12，核多可调大）。

### BWA-MEM 算法参数 | BWA-MEM algorithm

**通俗理解|In plain words:** 这一组直接透传 BWA-MEM 的底层参数（种子长度 `--bwa-k`、带宽 `--bwa-w`、X-dropoff `--bwa-d`、种子因子 `--bwa-r`、种子出现次数阈值 `--bwa-c`、链丢弃比例 `--bwa-D`、配对拯救轮数 `--bwa-m` 等）。**除非你明确知道在调什么，否则全部保持默认**——这些是 BWA 官方默认值，绝大多数数据效果最好。

### BWA-MEM 打分参数 | BWA-MEM scoring

**通俗理解|In plain words:** 匹配得分 `--bwa-A`、错配罚分 `--bwa-B`、gap 开放/延伸罚分 `--bwa-O`/`--bwa-E`、末端剪切罚分 `--bwa-L`、未配对罚分 `--bwa-U`。这些决定"给一段比对打多少分、低于多少分不报"。**同样一般不用动**；做特殊应用（如长读、特定变异）时才按 BWA 手册调。

### 后处理 | Post-processing

**通俗理解|In plain words:** `--markdup` 标记重复读段，`--remove-dup` 直接移除（需要先 markdup 语义）；`--keep-sam` 保留中间 SAM 文件（默认删）。`--markdup` 后最终 BAM 名为 `样本.markdup.bam`，否则为 `样本.bam`。做变异检测前一般建议 `--markdup`。

### 覆盖度与窗口 | Coverage & window

**通俗理解|In plain words:** `--min-base-quality` / `--min-mapping-quality` / `--max-depth` 传给 `samtools depth` 做过滤（默认 0 即不过滤）；`--window-size`（默认 1Mb）与 `--step-size`（默认 100kb）决定滑窗覆盖度的窗口大小与步长。想画全基因组深度图时按需调窗口。

## 分析流程 | Pipeline

~~~text
参考基因组(-g) + FASTQ目录(-i)
    │
    ▼
步骤1: 检查依赖(bwa + samtools)
步骤2: 检查/构建 BWA 基因组索引(.amb/.ann/.bwt/.pac/.sa)
步骤3: 按模式查找成对 FASTQ(样本, read1, read2)
步骤4: 逐样本比对:
   bwa mem → SAM → BAM(view) → 排序(sort) → (可选markdup) → 索引(index)
步骤5: 生成比对统计(flagstat + stats)
步骤6: 覆盖度分析(samtools depth + 滑窗)
步骤7: 生成汇总报告 alignment_summary.txt
~~~

## 输出 | Output

~~~text
bwa_output/
├── bam/
│   ├── sampleA.bam               # 最终BAM(或 sampleA.markdup.bam)
│   └── sampleA.bam.bai           # BAM 索引
├── coverage/
│   └── sampleA.depth.txt         # 单碱基深度(chrom pos depth)
├── windows/
│   └── sampleA.window.bed        # 滑窗平均深度(chrom start end avg_depth)
├── stats/
│   ├── sampleA.flagstat.txt      # samtools flagstat 比对统计
│   └── sampleA.stats.txt         # samtools stats 详细统计
├── logs/
│   └── bwament.log               # 运行日志
└── alignment_summary.txt         # 汇总报告(参数+已处理样本清单)
~~~

- 中间文件（SAM、raw BAM、sorted BAM）默认清理；`--keep-sam` 保留 SAM

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 每个样本先看 `stats/*.flagstat.txt` 的比对率，再看 `windows/*.window.bed` 的深度是否均匀。

- `flagstat` 中 `mapped` 比例越高越好（通常 >95% 为正常）；`properly paired` 比例高说明双端比对正常
- `coverage/*.depth.txt`：每行 `染色体 位置 深度`，可画深度曲线；`windows/*.window.bed` 是滑窗平均，适合看大尺度覆盖
- `stats/*.stats.txt`：含 GC、插入片段、错误率等更细的质控数据
- `alignment_summary.txt`：本次运行的参数与成功处理的样本列表，用于留档

## 参数选择建议 | Parameter Guidance

- **常规重测序比对**：默认参数即可，`biopytools bwa -g ref.fa -i fastq_dir`
- **下游做变异检测**：加 `--markdup`
- **中途断掉重跑**：加 `--resume`（自动跳过已产出最终 BAM 的样本）
- **要保留 SAM 排查**：加 `--keep-sam`
- **画全基因组深度图**：按需设 `--window-size / --step-size`（如 100kb 窗口、10kb 步长）

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 参考基因组文件｜Reference genome file |
| `--input, -i` | 必填 |  | 输入FASTQ目录｜Input FASTQ directory |
| `--pattern, -p` | `_1.clean.fq.gz` | str | FASTQ文件匹配模式｜FASTQ file pattern |
| `--output-dir, -o` | `./bwa_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--bwa-k` | `19` | int | 最小种子长度｜Minimum seed length |
| `--bwa-w` | `100` | int | 带宽｜Band width |
| `--bwa-d` | `100` | int | X-dropoff｜Off-diagonal X-dropoff |
| `--bwa-r` | `1.5` | float | 内部种子因子｜Internal seed factor |
| `--bwa-c` | `500` | int | 种子出现次数阈值｜Seed occurrence threshold |
| `--bwa-D` | `0.5` | float | 短链丢弃比例｜Short chain drop fraction |
| `--bwa-W` | `0` | int | 最小链长｜Minimum chain length |
| `--bwa-m` | `50` | int | 配对拯救轮数｜Mate rescue rounds |
| `--bwa-S` | — |  | 跳过配对拯救｜Skip mate rescue |
| `--bwa-P` | — |  | 跳过配对｜Skip pairing |
| `--bwa-A` | `1` | int | 匹配得分｜Match score |
| `--bwa-B` | `4` | int | 错配罚分｜Mismatch penalty |
| `--bwa-O` | `6,6` |  | Gap开放罚分｜Gap open penalty |
| `--bwa-E` | `1,1` |  | Gap延伸罚分｜Gap extension penalty |
| `--bwa-L` | `5,5` |  | 末端剪切罚分｜Clipping penalty |
| `--bwa-U` | `17` | int | 未配对罚分｜Unpaired penalty |
| `--bwa-M` | — |  | 标记次要比对｜Mark shorter split hits as secondary |
| `--bwa-T` | `30` | int | 最小输出得分｜Minimum score to output |
| `--bwa-a` | — |  | 输出所有比对｜Output all alignments |
| `--bwa-C` | — |  | 附加FASTQ注释｜Append FASTA/FASTQ comment |
| `--bwa-V` | — |  | 输出参考序列头｜Output reference FASTA header |
| `--bwa-Y` | — |  | 软剪切补充比对｜Soft clipping for supplementary alignments |
| `--markdup` | — |  | 标记重复序列｜Mark duplicate reads |
| `--remove-dup` | — |  | 移除重复序列｜Remove duplicate reads |
| `--min-base-quality` | `0` | int | 最小碱基质量｜Minimum base quality |
| `--min-mapping-quality` | `0` | int | 最小比对质量｜Minimum mapping quality |
| `--max-depth` | `0` | int | 最大深度限制｜Max depth limit |
| `--window-size` | `1000000` | int | 窗口大小｜Window size in bp |
| `--step-size` | `100000` | int | 步长｜Step size in bp |
| `--resume` | — |  | 断点续传｜Resume skip completed samples |
| `--keep-sam` | — |  | 保留SAM文件｜Keep SAM files |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 参考基因组文件｜Reference genome file |
| `-i, --input` | 必填 |  | 输入FASTQ文件目录｜Input FASTQ directory |
| `-p, --pattern` | `_1.clean.fq.gz` |  | FASTQ文件匹配模式｜FASTQ file pattern (default: _1.clean.fq.gz) |
| `-o, --output-dir` | `./bwa_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `--bwa-k` | `19` | int | 最小种子长度｜Minimum seed length |
| `--bwa-w` | `100` | int | 带宽｜Band width |
| `--bwa-d` | `100` | int | X-dropoff｜Off-diagonal X-dropoff |
| `--bwa-r` | `1.5` | float | 内部种子因子｜Internal seed factor |
| `--bwa-c` | `500` | int | 种子出现次数阈值｜Seed occurrence threshold |
| `--bwa-D` | `0.5` | float | 短链丢弃比例｜Short chain drop fraction |
| `--bwa-W` | `0` | int | 最小链长｜Minimum chain length |
| `--bwa-m` | `50` | int | 配对拯救轮数｜Mate rescue rounds |
| `--bwa-S` | — | store_true | 跳过配对拯救｜Skip mate rescue |
| `--bwa-P` | — | store_true | 跳过配对｜Skip pairing |
| `--bwa-A` | `1` | int | 匹配得分｜Match score |
| `--bwa-B` | `4` | int | 错配罚分｜Mismatch penalty |
| `--bwa-O` | `6,6` |  | Gap开放罚分｜Gap open penalty |
| `--bwa-E` | `1,1` |  | Gap延伸罚分｜Gap extension penalty |
| `--bwa-L` | `5,5` |  | 末端剪切罚分｜Clipping penalty |
| `--bwa-U` | `17` | int | 未配对罚分｜Unpaired penalty |
| `--bwa-M` | — | store_true | 标记次要比对｜Mark shorter split hits as secondary |
| `--bwa-T` | `30` | int | 最小输出得分｜Minimum score to output |
| `--bwa-a` | — | store_true | 输出所有比对｜Output all alignments |
| `--bwa-C` | — | store_true | 附加FASTQ注释｜Append FASTA/FASTQ comment |
| `--bwa-V` | — | store_true | 输出参考序列头｜Output reference FASTA header |
| `--bwa-Y` | — | store_true | 软剪切补充比对｜Soft clipping for supplementary alignments |
| `--markdup` | — | store_true | 标记重复序列｜Mark duplicate reads |
| `--remove-dup` | — | store_true | 移除重复序列｜Remove duplicate reads |
| `--min-base-quality` | `0` | int | 最小碱基质量｜Minimum base quality |
| `--min-mapping-quality` | `0` | int | 最小比对质量｜Minimum mapping quality |
| `--max-depth` | `0` | int | 最大深度限制｜Max depth limit |
| `--window-size` | `1000000` | int | 窗口大小｜Window size in bp |
| `--step-size` | `100000` | int | 步长｜Step size in bp |
| `--resume` | — | store_true | 断点续传｜Resume skip completed samples |
| `--keep-sam` | — | store_true | 保留SAM文件｜Keep SAM files |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- `bwa`、`samtools`（自动解析 align 域环境并经 conda run 调用；可用环境变量 BWA_PATH / SAMTOOLS_PATH 覆盖；域环境缺失时回退 PATH 直接调用）

## 常见问题 | FAQ

**Q1：会不会断点续传？**
基因组索引会自动复用（缺了才重建）。样本级断点续传需要显式加 `--resume`：已存在最终 BAM 的样本会跳过。不加 `--resume` 则每次全部重跑。

**Q2：提示"未找到匹配的FASTQ文件对"？**
read1 文件名必须包含 `-p` 指定的模式（默认 `_1.clean.fq.gz`），且同目录有对应的 `_2.clean.fq.gz`。后缀不同就改 `-p`。

**Q3：markdup 之后 BAM 名字变了？**
是的。`--markdup` 时最终文件是 `样本.markdup.bam`，否则是 `样本.bam`。`--remove-dup` 会直接删掉重复读段（比标记更激进）。

**Q4：线程数默认到底是多少？**
通过 `biopytools bwa` 入口时默认 12 线程（模块直调默认 88）。核多的机器建议显式 `-t 32` 或更高。

**Q5：`--bwa-*` 参数要调吗？**
绝大多数不用。这些是 BWA-MEM 的官方默认值，只有明确知道在做什么（如特殊读长、特定应用）时才按 BWA 手册调整。
