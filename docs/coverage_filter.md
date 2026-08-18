# 覆盖度过滤 | Coverage Filter

一句话理解：**根据测序 reads 在每段序列上的比对覆盖度，把基因组里的序列分成「高质量 / 中等质量 / 低质量」三档，并把它们分别抽出来**，帮你快速剔除可疑的低覆盖序列、留下值得信任的序列继续下游分析。

## 功能概述 | Overview

- 以 BAM 比对文件的覆盖度为依据，对基因组 FASTA 中的每条序列做质量分级
- 分三档：高质量（覆盖度 >= 高阈值）、中等质量（介于两阈值之间）、低质量（低于低阈值）
- 每档各输出一份 FASTA 文件和一份序列 ID 清单，方便单独取用
- 附带覆盖度明细表和一份汇总报告，可追溯每条序列的覆盖度数值
- 流程单次跑完、无断点续传（重跑会重新生成并覆盖旧产物）

## 快速开始 | Quick Start

```bash
biopytools coverage-filter -i sample.bam -f genome.fa -o filtered
```

最小输入：一个 BAM 文件 + 一个基因组 FASTA 文件 + 一个输出前缀（结果默认写到当前目录）。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>Plain meaning |
|------|----------|
| BAM | 记录「每条 read 比对到基因组哪个位置」的文件，像一份查岗记录表 |
| 覆盖度 coverage | 某段序列被 reads 覆盖到的比例（0-100%），相当于「这段序列被多少条 read 覆盖到」 |
| FASTA | 存放序列（DNA 字母串）的标准文本格式，每条序列有一个 `>` 开头的名字 |
| contig / 序列 | 基因组拼接后得到的一条条序列片段，本工具逐条给它们打分 |
| seqtk / seqkit | 两个常用的序列处理小工具，分别负责「按名单抽序列」和「做统计」 |

## 输入 | Input

### BAM 文件

标准 BAM 比对文件，需能正常用于 `samtools coverage`（一般比对流程的最终产物即可）。

### 基因组 FASTA 文件

与 BAM 对应的参考基因组 FASTA。序列名须与 BAM 中的参考序列名一致，这是「按名抽序列」能对上号的前提。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 这三样缺一不可：BAM 提供「每段序列被覆盖了多少」的证据，FASTA 提供「要抽取的序列本体」，输出前缀决定所有结果文件叫什么名字。注意 `-o` 不是输出目录，而是文件名前缀（填 `filtered` 会得到 `filtered_high_quality.fa` 等），这点容易和别的工具混淆。

### 输出位置 | Output location

**通俗理解|In plain words:** `-d` 才是输出目录，默认 `.`（当前目录）。一般不用动；只有想把结果集中放到某个文件夹时才指定。

### 线程与性能 | Threads

**通俗理解|In plain words:** `-t` 只影响最后一步 seqkit 统计的速度，不参与覆盖度计算。默认 12，一般不用动；机器核少时调小避免抢占资源。

### 覆盖度阈值 | Coverage thresholds

**通俗理解|In plain words:** 这两个数决定「好 / 中 / 差」的分界线。`--high-cov`（默认 90）是高质量的门槛，覆盖度达到 90% 以上的序列才算高质量；`--medium-cov-min`（默认 50）是中等质量的下限，覆盖度在 50% 到 90% 之间算中等，低于 50% 算低质量。调大 `--high-cov` 会收紧「高质量」标准（留下的高质量序列更少但更可靠）；调大 `--medium-cov-min` 会让更多序列掉进「低质量」。绝大多数项目用默认 90 / 50 即可，只有对覆盖度有特殊要求时才需要动。

## 分析流程 | Pipeline

```text
输入 BAM + FASTA
  |
  v
1. samtools coverage 计算每条序列覆盖度 -> <prefix>_coverage.txt
  |
  v
2. 按阈值分成高/中/低三档，各生成一份序列 ID 清单 (.list)
  |
  v
3. seqtk subseq 按清单从 FASTA 抽取序列 -> 三份 .fa
  |
  v
4. seqkit stats 统计 + 生成汇总报告 -> <prefix>_report.txt
```

## 输出 | Output

```text
输出目录/
├── filtered_coverage.txt          # 每条序列的覆盖度明细(samtools coverage 原样输出)
├── filtered_high_quality.list     # 高质量序列 ID 清单
├── filtered_medium_quality.list   # 中等质量序列 ID 清单
├── filtered_low_quality.list      # 低质量序列 ID 清单
├── filtered_high_quality.fa       # 高质量序列 FASTA(推荐下游使用)
├── filtered_medium_quality.fa     # 中等质量序列 FASTA
├── filtered_low_quality.fa        # 低质量序列 FASTA(可疑,建议复查)
├── filtered_report.txt            # 汇总报告
└── coverage_filter.log            # 运行日志(位于输出前缀的父目录)
```

（以上文件名中的 `filtered` 即 `-o` 传入的前缀。）

- `filtered_coverage.txt`：`samtools coverage` 的原始输出，其中的 `coverage` 列（第 6 列）即本工具用来分级的覆盖度百分比
- `filtered_*_quality.list`：各类序列的 ID 清单，每行一个序列名
- `filtered_*_quality.fa`：按清单抽出的对应序列 FASTA
- `filtered_report.txt`：配置参数、三档序列数量、seqkit 统计结果、输出文件说明的总览

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 打开 `filtered_report.txt` 看「三档各有多少条序列」就能判断数据质量：高质量占比高说明数据干净，低质量占比高说明比对或组装可能有问题。

- **三档数量**：报告里的「高质量序列 / 中等质量序列 / 低质量序列」行。高质量越多越好；低质量序列过多（比如过半）说明大量序列没被 reads 覆盖到，需检查比对步骤或测序深度
- **覆盖度明细**：`filtered_coverage.txt` 的 `coverage` 列是每条序列的实际覆盖度，可定位具体哪条序列偏低
- **好坏判据**：一般希望「高质量」是最大的一档；若「低质量」最多，先确认 BAM 与 FASTA 是否匹配（序列名是否一致）

## 参数选择建议 | Parameter Guidance

- **默认 90 / 50 即可**：大多数组装项目的覆盖度分级用默认值就够
- **想要更严的高质量**：把 `--high-cov` 提到 95，只保留几乎完全覆盖的序列
- **想让中等质量更宽**：把 `--medium-cov-min` 调到 30，让更多序列落入中等档
- **只想留两档**：可把 `--medium-cov-min` 设为 0，让序列要么高质量、要么低质量

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--bam-file, -i` | 必填 |  | BAM文件｜BAM file path |
| `--fasta-file, -f` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--output-prefix, -o` | 必填 |  | 输出文件前缀｜Output file prefix |
| `--output-dir, -d` | `.` |  | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--high-cov` | `90.0` | float | 高质量覆盖度阈值｜High quality coverage threshold |
| `--medium-cov-min` | `50.0` | float | 中等质量最小覆盖度｜Medium quality minimum coverage |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --bam-file` | 必填 |  | BAM文件路径｜BAM file path |
| `-f, --fasta-file` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output-prefix` | 必填 |  | 输出文件前缀｜Output file prefix |
| `-d, --output-dir` | `.` |  | 输出目录｜Output directory (default: current directory) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--high-cov` | `90.0` | float | 高质量覆盖度阈值｜High quality coverage threshold (default: 90.0) |
| `--medium-cov-min` | `50.0` | float | 中等质量最小覆盖度｜Medium quality minimum coverage (default: 50.0) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- samtools（`samtools coverage`，自动解析 align 域环境）
- seqkit（`seqkit stats`，自动解析 misc 域环境）
- seqtk（`seqtk subseq`，无对应功能域环境）

三个软件经 conda run 自动检测包装；可用环境变量 SAMTOOLS_PATH / SEQKIT_PATH / SEQTK_PATH 覆盖；域环境缺失时回退 PATH 直接调用。

## 常见问题 | FAQ

**Q1：换阈值重跑，为什么结果没变？**
本工具单次跑完、没有断点续传，但每次都会重新生成覆盖度文件。若结果没变，确认是否真的改了参数、以及是否写进了正确的输出目录（注意 `-o` 是前缀、`-d` 才是目录）。

**Q2：日志文件 coverage_filter.log 在哪？**
日志写在「输出前缀的父目录」下。若 `-o filtered`（不带路径），日志落在当前目录；若 `-o sub/filtered`，日志落在 `sub/` 下。这与结果文件（写在 `-d` 指定目录）可能不是同一个位置。

**Q3：提示缺少依赖怎么办？**
程序会检查 samtools、seqtk、seqkit 三个命令是否在 PATH 中，缺哪个装哪个（可用 conda 分别安装，或确认它们所在目录已加入 PATH）。

**Q4：三档序列加起来不等于总数？**
理论上三档（高 >=90、中 50~90、低 <50）互斥且覆盖全部序列。若对不上，检查 `filtered_coverage.txt` 中是否有异常行（如序列名缺失或覆盖度列为空）。
