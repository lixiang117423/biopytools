
# Fastq2VCF (GTX) - 基于 GTX 的全流程变异检测 | Fastq to VCF (GTX) pipeline

一句话理解：**从原始 FASTQ 出发，经质控、GTX 比对、联合变异检测、过滤，一键得到干净的 SNP/INDEL VCF**，解决「大批量重测序样本从 reads 到变异」的端到端问题。

## 功能概述 | Overview { #overview }

- 五步流水线：fastp 质控 → GTX 基因组索引 → GTX WGS 比对 → 联合变异检测 → SNP/INDEL 过滤
- GTX 一体化完成比对+变异检测；样本数小于 200 用单机模式，大于等于 200 自动切集群模式（生成操作指南手动投递）
- 支持 --input（原始 reads）或 --clean-fastq-dir（已清洗 reads，自动跳过质控）二选一
- 默认启用断点续传，已完成的步骤/索引自动跳过
- 运行前系统预检查：必需工具(samtools/bcftools/biopytools/faketime)、磁盘(>=200GB)、内存(>=64GB)

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools fastq2vcf-gtx -i /path/to/raw_fastq -g /path/to/genome.fa -o /path/to/output
```

最小输入：一个放原始 FASTQ 的目录 + 一个参考基因组 FASTA（输出目录默认当前目录 .）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| GTX | 本流程的比对+变异检测引擎，一体化完成 |
| fastp | 质控工具，剪掉低质量碱基、去接头 |
| 基因组索引 | 给参考基因组建「目录」，让比对更快 |
| gVCF | 记录全基因组每个位点证据的中间文件 |
| 联合变异检测(joint calling) | 把所有样本的 gVCF 合起来一起判定变异 |
| 分块窗口(window size) | 把基因组切成多大一块一块处理（默认 20Mb） |
| 单机模式/集群模式 | 样本少在单台机器跑，样本多自动拆到集群 |
| 断点续传 | 已完成的步骤重跑时自动跳过 |

## 输入 | Input { #input }

- **原始 FASTQ 目录**（`-i`）与 **已清洗 FASTQ 目录**（`--clean-fastq-dir`）**二选一**：提供 --clean-fastq-dir 时自动跳过质控
- **参考基因组**（`-g`）：FASTA，必需

R1/R2 配对模式默认 _1.fq.gz / _2.fq.gz（可用 --read1-pattern-fastp/--read2-pattern-fastp 改）。

## 参数说明 | Parameters { #parameters }

### 必需与路径 | Required & paths

**通俗理解|In plain words:** -i（原始 reads）和 --clean-fastq-dir（清洁 reads）二选一，-g 参考基因组必填。--mapping-dir/--gvcf-dir/--bam-dir/--joint-dir/--filter-dir 可自定义各步骤输出位置，**不填会自动在输出目录下建 01_fastp 到 05_filtered_snp_indel**，一般不用动。

### 过滤参数 | Filtering

**通俗理解|In plain words:** --min-depth（默认 5）和 --quality（默认 30）是 SNP/INDEL 共同的「深度/质量」门槛；--indel-min-dp/--indel-min-qual 可给 INDEL 单独设更严/更松的标准。调大=更严格、变异更少但更可信。**默认值适合常规重测序，一般不用动。**

### GTX 参数 | GTX options

**通俗理解|In plain words:** 这是 GTX 引擎的旋钮。--gtx-single-threshold（默认 200）决定「多少样本以下在单机跑」；--gtx-window-size（默认 20Mb）是分块大小；--gtx-pcr-indel-model、--gtx-min-confidence、--gtx-min-base-qual、--gtx-ploidy 影响 INDEL 判定与置信度。**默认值经实践验证，一般不用动；非二倍体物种才需改 --gtx-ploidy。**

### 执行控制 | Execution

**通俗理解|In plain words:** --step 1~5 单独跑某一步；--no-checkpoint 关闭断点续传；--force 强制覆盖；--dry-run 只打印命令不执行；--keep-intermediate 保留中间文件。**日常用默认即可。**

## 分析流程 | Pipeline { #pipeline }

```text
原始 FASTQ 目录 + 参考基因组
    │
    ▼
步骤1: fastp 质控 (→ 01_fastp/)
    │
    ▼
步骤2: GTX 基因组索引 (→ 02_genome_index/)
    │
    ▼
步骤3: GTX WGS 比对 (→ 03_mapping/bam + 03_mapping/vcf)
    │
    ▼
步骤4: 联合变异检测 (→ 04_joint_calling/)
    │   (样本>=200 时切换集群模式,生成操作指南手动投递)
    ▼
步骤5: SNP/INDEL 过滤 (→ 05_filtered_snp_indel/)
    │
    ▼
最终报告
```

## 输出 | Output { #output }

```text
output/
├── 01_fastp/                      # fastp 清洗后的 FASTQ
├── 02_genome_index/               # GTX 基因组索引
├── 03_mapping/
│   ├── bam/                       # 各样本比对 BAM
│   └── vcf/                       # 各样本 gVCF
├── 04_joint_calling/              # 联合检测结果 (最终 raw VCF)
├── 05_filtered_snp_indel/         # 过滤后的 SNP/INDEL VCF (最终结果)
└── 99_logs/                       # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. 05_filtered_snp_indel/（最终结果）

**通俗理解|In plain words:** 这是最终交付的、过滤后的 SNP/INDEL VCF（.vcf.gz），直接用于下游分析（GWAS、群体遗传等）。文件按变异类型/染色体拆分。

### 2. 04_joint_calling/（联合检测原始结果）

**通俗理解|In plain words:** 过滤前的完整变异集合，想用不同阈值重新过滤时从这里出发。

### 3. 99_logs/（日志）

**通俗理解|In plain words:** 每一步的完整命令与进度，排查问题时先看这里。

## 参数选择建议 | Parameter Guidance { #guidance }

- **--gtx-single-threshold**：样本小于 200 保持默认（单机）；大于等于 200 会自动切集群模式，无需手动改。
- **--min-depth / --quality**：低深度数据可适当下调（如 depth 3）；想要更保守的结果可上调。
- **--gtx-ploidy**：非二倍体（如多倍体作物）改成对应倍性。
- **--clean-fastq-dir**：已经做过质控的 reads，直接传清洁目录，省去质控步骤。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | — | Path | 原始FASTQ目录(与 --clean-fastq-dir 二选一)｜Raw FASTQ directory (mutually exclusive with --clean-fastq-dir) |
| `--genome, -g` | 必填 | Path | 参考基因组｜Reference genome file path |
| `--output-dir, -o` | `.` | Path | 输出目录｜Output directory path |
| `--clean-fastq-dir` | — | Path | 清洗后FASTQ目录｜Clean FASTQ files directory path |
| `--mapping-dir` | — | Path | 比对结果目录｜Mapping results directory path |
| `--gvcf-dir` | — | Path | gVCF文件目录｜gVCF files directory path |
| `--bam-dir` | — | Path | BAM文件目录｜BAM files directory path |
| `--joint-dir` | — | Path | 联合检测目录｜Joint calling results directory path |
| `--filter-dir` | — | Path | 过滤结果目录｜Filtering results directory path |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--min-depth, --snp-min-dp` | `5` | int | 最小测序深度｜Minimum depth for SNP/InDel |
| `--quality, -q, --min-qual, --snp-min-qual` | `30` | int | 最小质量值｜Minimum quality for SNP/InDel |
| `--indel-min-dp` | `5` | int | InDel最小深度｜InDel minimum depth |
| `--indel-min-qual` | `30` | int | InDel最小质量｜InDel minimum quality |
| `--gtx-single-threshold` | `200` | int | GTX单机样本阈值｜GTX single machine sample count threshold |
| `--gtx-window-size` | `20000000` | int | GTX窗口大小｜GTX chunk window size in bp |
| `--gtx-bin` | `~/software/gtx/bin/gtx` | Path | GTX可执行文件路径｜GTX executable path |
| `--use-gtx-wgs` | `True` |  | 使用GTX WGS模式｜Use GTX WGS |
| `--gtx-pcr-indel-model` | `CONSERVATIVE` |  | GTX PCR InDel模型｜GTX PCR InDel model |
| `--gtx-min-confidence` | `30` | int | GTX最小置信度｜GTX minimum confidence |
| `--gtx-min-base-qual` | `20` | int | GTX最小碱基质量｜GTX minimum base quality |
| `--gtx-ploidy` | `2` | int | GTX倍性｜GTX ploidy |
| `--gtx-cmd-gen-script` | `${HOME}/software/scripts/51.生成GTX按染色体合并gVCF的脚本.sh` |  | GTX命令生成脚本｜GTX command generation script path |
| `--step, -s` | — | 1/2/3/4/5 | 运行指定步骤｜Run only specified step |
| `--no-checkpoint` | — |  | 禁用检查点恢复｜Disable checkpoint resume |
| `--dry-run` | — |  | 试运行｜Dry run |
| `--force, -f` | — |  | 强制覆盖｜Force overwrite |
| `--keep-intermediate` | — |  | 保留中间文件｜Keep intermediate files |
| `--verbose, -v` | — |  | 详细输出｜Verbose output |
| `--quiet` | — |  | 仅输出错误｜Quiet mode |
| `--log-file` | — | Path | 日志文件｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--skip-qc` | — |  | 跳过质控｜Skip QC |
| `--skip-mapping` | — |  | 跳过比对｜Skip mapping |
| `--read1-pattern-fastp` | `_1.fq.gz` |  | R1文件模式｜QC R1 file pattern |
| `--read2-pattern-fastp` | `_2.fq.gz` |  | R2文件模式｜QC R2 file pattern |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  | 输入reads文件目录路径｜Input reads directory path |
| `-g, --genome` | 必填 |  | 参考基因组文件路径｜Reference genome file path |
| `-o, --output-dir` | `.` |  | 输出目录路径｜Output directory path |
| `--clean-fastq-dir` | — |  | 清洁FASTQ文件目录路径｜Clean FASTQ files directory path |
| `--mapping-dir` | — |  | 比对结果目录路径｜Mapping results directory path |
| `--gvcf-dir` | — |  | gVCF文件目录路径｜gVCF files directory path |
| `--bam-dir` | — |  | BAM文件目录路径｜BAM files directory path |
| `--joint-dir` | — |  | 联合检测结果目录路径｜Joint calling results directory path |
| `--filter-dir` | — |  | 过滤结果目录路径｜Filtering results directory path |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--min-depth, --snp-min-dp` | `5` | int | SNP/InDel最小深度｜Minimum depth for SNP/InDel |
| `-q, --quality, --min-qual, --snp-min-qual` | `30` | int | SNP/InDel最小质量｜Minimum quality for SNP/InDel |
| `--indel-min-dp` | `5` | int | InDel最小深度｜InDel minimum depth |
| `--indel-min-qual` | `30` | int | InDel最小质量｜InDel minimum quality |
| `--gtx-single-threshold` | `200` | int | GTX单机模式样本数阈值｜GTX single machine sample count threshold |
| `--gtx-window-size` | `20000000` | int | GTX分块窗口大小｜GTX chunk window size in bp |
| `--gtx-bin` | `~/software/gtx/bin/gtx` |  | GTX可执行文件路径｜GTX executable path |
| `--gtx-cmd-gen-script` | `${HOME}/software/scripts/51.生成GTX按染色体合并gVCF的脚本.sh` |  | GTX命令生成脚本路径｜GTX command generation script path |
| `--use-gtx-wgs` | `True` | store_true | 使用GTX WGS｜Use GTX WGS |
| `--no-gtx-wgs` | — | store_false | 禁用GTX WGS｜Disable GTX WGS |
| `--gtx-pcr-indel-model` | `CONSERVATIVE` |  | GTX PCR InDel模型｜GTX PCR InDel model |
| `--gtx-min-confidence` | `30` | int | GTX最小置信度｜GTX minimum confidence |
| `--gtx-min-base-qual` | `20` | int | GTX最小碱基质量｜GTX minimum base quality |
| `--gtx-ploidy` | `2` | int | GTX倍性｜GTX ploidy |
| `-f, --force` | — | store_true | 强制覆盖已存在文件｜Force overwrite existing files |
| `--dry-run` | — | store_true | 模拟运行｜Dry run |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `--step` | — | 1/2/3/4/5 | 只运行指定步骤｜Run only specified step |
| `--no-checkpoint` | — | store_true | 禁用断点续传｜Disable checkpoint resume |
| `-v, --verbose` | `0` | count | 详细输出模式｜Verbose mode |
| `--quiet` | — | store_true | 静默模式｜Quiet mode |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--skip-qc` | — | store_true | 跳过质控步骤｜Skip QC step |
| `--skip-mapping` | — | store_true | 跳过比对步骤｜Skip mapping step |
| `--read1-pattern-fastp` | `_1.fq.gz` |  | 质控R1文件匹配模式｜QC R1 file pattern |
| `--read2-pattern-fastp` | `_2.fq.gz` |  | 质控R2文件匹配模式｜QC R2 file pattern |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- GTX（默认 ~/software/gtx/bin/gtx，可用 --gtx-bin 覆盖）
- fastp（质控）
- samtools、bcftools（必需，预检查会校验）
- biopytools（本工具自身，预检查会校验）
- faketime（**必需**：GTX 所有命令靠它绕过 license 时间校验）
- Python 3

## 常见问题 | FAQ { #faq }

**Q1：为什么要 faketime？**
GTX 的所有命令（index/wgs/joint）依赖 faketime 绕过 license 时间校验；缺少它会直接报「缺少必需工具」并退出。

**Q2：样本很多时怎么跑？**
样本数达到 --gtx-single-threshold（默认 200）时，联合检测步骤会切换集群模式并生成操作指南，需按指南手动投递，这不算是失败。

**Q3：磁盘/内存要求？**
预检查要求输出目录所在盘至少 200GB 空闲、内存 >=64GB，不够会直接退出。

**Q4：断点续传怎么工作？**
默认开启（--no-checkpoint 关闭）。已生成的索引、质控结果、比对结果存在即跳过；换参数重跑前用 --force 强制覆盖。

**Q5：已做过质控的 reads 怎么用？**
传 --clean-fastq-dir 而不是 -i，程序会自动跳过质控步骤。
