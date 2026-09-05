# GTX WGS 全流程分析 | GTX WGS Pipeline

一句话理解：**把一堆质控好的 FASTQ 测序数据(多个样本)一次性跑完「比对→变异检测→群体联合」，得到每个样本的 BAM 和 gVCF，以及可选的多样本联合 VCF**，解决「从原始测序数据到最终变异结果」的端到端问题。

## 功能概述 | Overview { #overview }

- 输入一个含 clean FASTQ 的目录，自动按 R1/R2 文件名模式配对样本
- 自动检查/构建参考基因组索引(bwa index + samtools faidx)，缺失时补建
- 每个样本用 GTX 的 `wgs` 子命令跑变异检测(产出 BAM + gVCF)
- 可选 `--enable-joint` 做多样本联合变异检测(gtx joint)，产出合并 VCF
- 断点续传：已完成的样本跳过，联合结果已存在也跳过
- 生成分析总结报告 `gtx_analysis_summary.txt`

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools gtx -i input_dir -o results/ -r ref.fa --enable-joint
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| WGS | 全基因组测序，把整个基因组都测一遍 |
| clean FASTQ | 质控后干净的双端测序数据，R1/R2 是一对「前后两段读序」 |
| read group | 比对时给每个样本打的「身份标签」，告诉下游这段读序属于哪个样本 |
| gVCF(genomic VCF) | 每个位点都有记录(包括没变异的位点)的 VCF，信息更全，适合多样本联合 |
| joint calling | 把多个样本的变异放在一起重新判读，能提高低频变异的检出率 |
| faketime | 一个伪装系统时间的工具，用于绕过 GTX 的 license 时间校验 |
| 参考索引 | 参考基因组为快速比对预生成的一堆辅助文件(像书的目录) |

## 输入 | Input { #input }

### FASTQ 目录

`-i` 指向一个目录，里面放每个样本的一对 clean FASTQ。默认按 `*_1_clean.fq.gz` / `*_2_clean.fq.gz` 配对(可用 `--read1-pattern` / `--read2-pattern` 改)。

```text
input_dir/
├── sample1_1_clean.fq.gz
├── sample1_2_clean.fq.gz
├── sample2_1_clean.fq.gz
└── sample2_2_clean.fq.gz
```

### 参考基因组

`-r` 指定参考基因组 FASTA。索引若缺失会自动用 bwa index + samtools faidx 构建。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** `-i` 是放 clean FASTQ 的目录，`-o` 是结果目录，`-r` 是参考基因组。三个都要给。

### 运行参数 | Runtime

**通俗理解|In plain words:** `-t` 控制每个样本比对/变异检测用的线程数，默认 12，调大更快但更吃 CPU。`--gtx-path` 是 GTX 可执行文件位置(默认 `~/software/gtx/GTX.CAT_2.2.1/bin/gtx`)，`--tmp-dir` 是临时目录(默认输出目录下的 `tmp/`)。**一般不用动。**

### Joint Calling 参数 | Joint calling

**通俗理解|In plain words:** `--enable-joint` 是开关，加了才会在单样本跑完后做多样本联合(默认关闭)。`--joint-output`(默认 `merged_gtx.vcf.gz`)是合并结果文件名，必须以 `.vcf.gz` 结尾；`--joint-threads`(默认 88)是联合这步单独用的线程数。**只有多样本、且想要统一的群体变异结果时才开。**

### 质控参数 | Quality control

**通俗理解|In plain words:** 这四个决定「什么样的变异算数」。`--min-confidence`(默认 30)是最低置信度，调高更严格、假阳性少但可能漏真变异；`--min-base-quality`(默认 20)是最低碱基质量；`--ploidy`(默认 2)是倍性，一般物种都是二倍体不用动；`--pcr-indel-model`(默认 CONSERVATIVE)是 PCR 引起的 indel 处理方式。**绝大多数项目用默认值即可，一般不用动。**

### 文件模式参数 | File patterns

**通俗理解|In plain words:** 告诉程序「R1/R2 文件怎么命名」。只有你的文件命名和默认 `*_1_clean.fq.gz` / `*_2_clean.fq.gz` 不一样时才需要改，且必须含 `*` 通配符。

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 先检查环境和参考索引，再逐个样本跑 GTX，最后(可选)做联合，收尾写总结报告。

```text
输入 clean FASTQ 目录 + 参考基因组
    |
    ▼
步骤0: 依赖检查(GTX + faketime)
    |
    ▼
步骤1: 参考索引检查/构建(bwa index + samtools faidx)
    |
    ▼
步骤2: 单样本分析(gtx wgs -g → bam/<样本>_sorted.bam + vcf/<样本>.vcf.gz)
    |
    ▼
步骤3: Joint Calling(可选, gtx joint → joint/merged_gtx.vcf.gz)
    |
    ▼
步骤4: 生成总结报告 gtx_analysis_summary.txt
```

## 输出 | Output { #output }

```text
results/
├── bam/                        # 每个样本的比对结果
│   └── sample1_sorted.bam
├── vcf/                        # 每个样本的 gVCF
│   └── sample1.vcf.gz
├── joint/                      # 仅 --enable-joint 时生成
│   ├── sample_map.txt          # 样本名到 VCF 的映射文件
│   └── merged_gtx.vcf.gz       # 多样本联合 VCF(核心群体结果)
├── tmp/                        # 临时文件目录
├── gtx_analysis_summary.txt    # 分析总结报告
└── gtx_analysis.log            # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. vcf/<样本>.vcf.gz(单样本 gVCF)

**通俗理解|In plain words:** 每个样本自己的变异记录(gVCF 格式，含未变异位点)。适合存档或喂给下游联合分析。

### 2. joint/merged_gtx.vcf.gz(联合 VCF)

**通俗理解|In plain words:** 把所有样本放在一起重新判读后的群体变异结果，是多样本项目最常用的下游输入(做 GWAS、群体遗传等)。

- `sample_map.txt` 记录了哪些样本的 VCF 参与了联合
- 联合需要至少 2 个成功处理的样本，样本不足时会跳过并提示

### 3. gtx_analysis_summary.txt(总结报告)

**通俗理解|In plain words:** 一张「跑了哪些样本、各文件多大」的收尾清单，含分析配置、质控参数、文件统计和每个样本的 BAM/VCF 大小。

- 检查 `分析完成状态` 里的单样本完成数和 Joint calling 状态，确认没有样本失败

## 参数选择建议 | Parameter Guidance { #guidance }

- **单样本只想要 gVCF**：不加 `--enable-joint`，最快
- **多样本要做群体分析**：加 `--enable-joint`，一次拿到个体 gVCF + 联合 VCF
- **样本多、联合这步慢**：调大 `--joint-threads`(默认 88)
- **数据质量差、假阳性多**：调高 `--min-confidence`(如 40)；想更保守可同时调高 `--min-base-quality`
- **非二倍体物种**：改 `--ploidy`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input-dir, -i` | 必填 |  | 输入目录(包含clean FASTQ文件)｜Input directory path (containing clean FASTQ files) |
| `--output-dir, -o` | 必填 | Path | 输出目录路径｜Output directory path |
| `--reference, -r` | 必填 |  | 参考基因组文件路径｜Reference genome file path |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--gtx-path` | `~/software/gtx/GTX.CAT_2.2.1/bin/gtx` | Path | GTX程序路径｜GTX program path |
| `--tmp-dir` | — | Path | 临时目录路径｜Temporary directory path |
| `--enable-joint` | — |  | 启用Joint Calling｜Enable Joint Calling for multi-sample variant detection |
| `--joint-output` | `merged_gtx.vcf.gz` |  | Joint calling输出VCF文件名｜Joint calling output VCF filename |
| `--joint-threads` | `88` | int | Joint calling使用的线程数｜Number of threads for joint calling |
| `--min-confidence` | `30` | int | 最小置信度阈值｜Minimum confidence threshold |
| `--min-base-quality` | `20` | int | 最小碱基质量阈值｜Minimum base quality threshold |
| `--ploidy` | `2` | int | 倍性｜Ploidy |
| `--pcr-indel-model` | `CONSERVATIVE` | CONSERVATIVE/AGGRESSIVE | PCR indel模型｜PCR indel model |
| `--read1-pattern` | `*_1_clean.fq.gz` | str | R1文件匹配模式｜R1 file pattern |
| `--read2-pattern` | `*_2_clean.fq.gz` | str | R2文件匹配模式｜R2 file pattern |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-dir` | 必填 |  | 输入目录路径(包含clean FASTQ文件)｜Input directory path (containing clean FASTQ files) |
| `-o, --output-dir` | 必填 |  | 输出目录路径｜Output directory path |
| `-r, --reference` | 必填 |  | 参考基因组文件路径｜Reference genome file path |
| `--force-rebuild-index` | — | store_true | 强制重新构建参考基因组索引｜Force rebuild reference genome index |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `--gtx-path` | — |  | GTX程序路径｜GTX program path |
| `--tmp-dir` | — |  | 临时目录路径｜Temporary directory path |
| `--enable-joint` | — | store_true | 启用Joint Calling｜Enable Joint Calling |
| `--joint-output` | `merged_gtx.vcf.gz` |  | Joint calling输出文件名｜Joint calling output filename |
| `--joint-threads` | `88` | int | Joint calling使用的线程数｜Number of threads for joint calling |
| `--min-confidence` | `30` | int | 最小置信度阈值｜Minimum confidence threshold |
| `--min-base-quality` | `20` | int | 最小碱基质量阈值｜Minimum base quality threshold |
| `--ploidy` | `2` | int | 倍性｜Ploidy |
| `--pcr-indel-model` | `CONSERVATIVE` |  | PCR indel模型｜PCR indel model |
| `--read1-pattern` | `*_1.fq.gz` |  | R1文件匹配模式｜R1 file pattern |
| `--read2-pattern` | `*_2.fq.gz` |  | R2文件匹配模式｜R2 file pattern |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- GTX(核心变异检测工具，默认路径 `~/software/gtx/GTX.CAT_2.2.1/bin/gtx`，可用 --gtx-path 覆盖；无对应功能域环境)
- bwa(参考基因组索引构建，自动解析 align 域环境并经 conda run 调用，可用环境变量 BWA_PATH 覆盖；域环境缺失时回退 PATH 直接调用)
- samtools(参考基因组索引构建，自动解析 align 域环境并经 conda run 调用，可用环境变量 SAMTOOLS_PATH 覆盖；域环境缺失时回退 PATH 直接调用)
- faketime(必须能在 PATH 中找到，GTX 命令依赖它绕过 license 时间校验)

## 常见问题 | FAQ { #faq }

**Q1：为什么报「faketime 未安装」？**
GTX 命令会用 `faketime '2020-10-20 00:00:00'` 伪装系统时间来绕过 license 时间校验，所以 faketime 必须能通过 `which` 找到。请先安装 libfaketime。

**Q2：支持断点续传吗？**
支持。单个样本按 `bam/<样本>_sorted.bam` 和 `vcf/<样本>.vcf.gz` 是否都存在判断跳过；Joint Calling 结果 `joint/merged_gtx.vcf.gz` 已存在也会跳过；参考索引缺失时自动补建。换参数重跑前需先删除旧产物。

**Q3：Joint Calling 没跑起来？**
需要满足两点：加了 `--enable-joint`，且至少 2 个样本成功处理。样本不足会在日志里给出 WARNING。另外 `--joint-output` 的文件名必须以 `.vcf.gz` 结尾。

**Q4：线程数到底是多少？**
通过 `biopytools gtx` 调用时 `-t` 默认 12；直接调用 `python -m biopytools.gtx` 时默认 88。Joint Calling 那步单独用 `--joint-threads`(默认 88)。

**Q5：`--force-rebuild-index` 是什么？**
它是模块直调(main.py)的参数，用于强制重建参考索引；click 包装器未暴露，用默认行为(索引存在则跳过)即可。
