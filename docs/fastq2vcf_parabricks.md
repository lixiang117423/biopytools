
# Fastq2VCF (Parabricks) - 基于 GPU 的全流程变异检测 | Fastq to VCF (Parabricks) pipeline

一句话理解：**与 GTX 版相同的「reads 到 VCF」全流程，但比对步骤用 GPU 加速的 Parabricks 完成**，解决「大批量样本重测序变异检测耗时太长」的问题。

## 功能概述 | Overview { #overview }

- 五步流水线：fastp 质控 → 基因组索引 → Parabricks(GPU) 比对 → 联合变异检测(GTX) → SNP/INDEL 过滤
- 按样本数自动选引擎：小于 4 样本用 GATK、4~200 用 GTX 单机、大于等于 200 切集群模式
- 项目目录式管理：所有输出统一挂在 -p project_base 下，结构清晰
- 默认启用断点续传；运行前系统预检查（bwa/samtools/bcftools/gatk/biopytools + 磁盘/内存）
- 需要 NVIDIA GPU（Parabricks 依赖）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools fastq2vcf-parabricks -i /path/to/raw_fastq -g /path/to/genome.fa -p /path/to/project
```

最小输入：原始 FASTQ 目录 + 参考基因组 FASTA + 一个项目根目录（-p，所有结果都放在这里）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| Parabricks | NVIDIA 的 GPU 版生信工具集，比对/变异检测比 CPU 快一个数量级 |
| GPU | 图形处理器，天生擅长大规模并行计算 |
| GTX | 本流程的联合变异检测引擎 |
| 项目根目录(project-base) | 一个文件夹，装下本项目的所有中间与最终结果 |
| 样本数阈值 | 程序按样本多少自动选择用哪种引擎跑 |

## 输入 | Input { #input }

- **原始 FASTQ 目录**（`-i`，必需）
- **参考基因组**（`-g`，必需）
- **项目根目录**（`-p`，必需）：所有输出挂在这个目录下

## 参数说明 | Parameters { #parameters }

### 必需与路径 | Required & paths

**通俗理解|In plain words:** -i（reads）、-g（参考基因组）、-p（项目目录）三个必填。--clean-fastq-dir/--mapping-dir/--joint-dir/--filter-dir 等可自定义，不填会自动在项目目录下建 01_data/clean、02_mapping、03_joint_calling、04_filtered_snp_indel。**一般不用动。**

### 计算资源 | Computing resources

**通俗理解|In plain words:** --threads-mapping/--threads-gtx/--threads-filter 分别控制比对、GTX、过滤三阶段的线程数（CLI 界面默认各 12）。比对主要靠 GPU，CPU 线程是辅助；按机器配置调。

### 过滤参数 | Filtering

**通俗理解|In plain words:** --snp-min-dp（默认 5）/--snp-min-qual（默认 30）是 SNP 的深度/质量门槛，--indel-min-dp/--indel-min-qual 给 INDEL。调大=更严格。**默认适合常规重测序，一般不用动。**

### 引擎选择阈值 | Engine thresholds

**通俗理解|In plain words:** --gatk-threshold（默认 4）决定「样本小于 4 时改用 GATK」；--gtx-single-threshold（默认 200）决定「样本小于 200 用 GTX 单机、否则切集群」。**基本不用改，程序会自动选择。**

### 执行控制 | Execution

**通俗理解|In plain words:** --step 1~5 单独跑某一步；--no-checkpoint 关断点续传；--force 强制覆盖；--dry-run 只打印不执行。**日常用默认。**

## 分析流程 | Pipeline { #pipeline }

```text
原始 FASTQ + 参考基因组 + 项目目录
    │
    ▼
步骤1: fastp 质控 (→ 01_data/clean)
    │
    ▼
步骤2: 基因组索引 (→ 02_mapping)
    │
    ▼
步骤3: Parabricks(GPU) 比对 (→ 02_mapping/bam + vcf)
    │
    ▼
步骤4: 联合变异检测 (→ 03_joint_calling；按样本数选 GATK/GTX单机/集群)
    │
    ▼
步骤5: SNP/INDEL 过滤 (→ 04_filtered_snp_indel)
    │
    ▼
最终报告 (ANALYSIS_REPORT.txt)
```

## 输出 | Output { #output }

```text
project/
├── 01_data/clean/                 # fastp 清洗后的 FASTQ
├── 02_mapping/
│   ├── bam/                       # 各样本比对 BAM
│   └── vcf/                       # 各样本 gVCF
├── 03_joint_calling/              # 联合检测结果
├── 04_filtered_snp_indel/         # 过滤后的 SNP/INDEL VCF (最终结果)
├── 99_logs/                       # 运行日志
└── ANALYSIS_REPORT.txt            # 最终分析报告
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. 04_filtered_snp_indel/（最终结果）

**通俗理解|In plain words:** 最终过滤后的 SNP/INDEL VCF，直接用于下游分析。

### 2. ANALYSIS_REPORT.txt（最终报告）

**通俗理解|In plain words:** 一页纸汇总本次运行的输入、参数、样本数、输出文件与大小，相当于回执单。

### 3. 99_logs/pipeline.log（日志）

**通俗理解|In plain words:** 每一步完整命令与进度，排查问题先看这里。

## 参数选择建议 | Parameter Guidance { #guidance }

- **GPU**：本流程核心依赖 NVIDIA GPU（Parabricks），没有 GPU 请改用 biopytools fastq2vcf-gtx。
- **--threads-mapping**：GPU 机器可适当调高；CPU 受限时对加速有限。
- **--snp-min-dp/--snp-min-qual**：低深度数据可下调，想保守可上调。
- **样本数阈值**：一般无需改，程序自动按样本数选择 GATK/GTX单机/集群。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--version, -V` | — |  | 显示版本信息｜Show version information |
| `--input, -i` | 必填 | Path | 原始FASTQ文件目录路径｜Raw FASTQ files directory path |
| `--genome, -g` | 必填 | Path | 参考基因组文件路径｜Reference genome file path |
| `--project-base, -p` | 必填 | Path | 项目根目录路径｜Project base directory path |
| `--clean-fastq-dir` | — | Path | 清洁FASTQ文件目录路径｜Clean FASTQ files directory path |
| `--mapping-dir` | — | Path | 比对结果目录路径｜Mapping results directory path |
| `--joint-dir` | — | Path | 联合检测结果目录路径｜Joint calling results directory path |
| `--filter-dir` | — | Path | 过滤结果目录路径｜Filtering results directory path |
| `--output-dir, -o` | — | Path | 输出目录路径｜Output directory path |
| `--threads-mapping` | `12` | int | 比对线程数｜Number of mapping threads |
| `--threads-gtx` | `12` | int | GTX线程数｜Number of GTX threads |
| `--threads-filter` | `12` | int | 过滤线程数｜Number of filtering threads |
| `--snp-min-dp` | `5` | int | SNP最小深度｜SNP minimum depth |
| `--snp-min-qual` | `30` | int | SNP最小质量｜SNP minimum quality |
| `--indel-min-dp` | `5` | int | InDel最小深度｜InDel minimum depth |
| `--indel-min-qual` | `30` | int | InDel最小质量｜InDel minimum quality |
| `--gatk-threshold` | `4` | int | GATK模式样本数阈值｜GATK sample count threshold |
| `--gtx-single-threshold` | `200` | int | GTX单机模式样本数阈值｜GTX single machine sample count threshold |
| `--gtx-window-size` | `20000000` | int | GTX分块窗口大小｜GTX chunk window size in bp |
| `--gtx-bin` | `~/software/gtx/bin/gtx` | Path | GTX可执行文件路径｜GTX executable path |
| `--step, -s` | — | 1/2/3/4/5 | 只运行指定步骤｜Run only specified step |
| `--no-checkpoint` | — |  | 禁用断点续传｜Disable checkpoint resume |
| `--dry-run` | — |  | 测试模式｜Test mode |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose output mode |
| `--quiet` | — |  | 静默模式｜Quiet mode |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--log-file` | — | Path | 日志文件路径｜Log file path |
| `--force, -f` | — |  | 强制覆盖已存在的结果｜Force overwrite existing results |
| `--skip-qc` | — |  | 跳过质控步骤｜Skip QC step |
| `--skip-mapping` | — |  | 跳过比对步骤｜Skip mapping step |
| `--read1-pattern-fastp` | `_1.fq.gz` |  | 质控R1文件匹配模式｜QC R1 file pattern |
| `--read2-pattern-fastp` | `_2.fq.gz` |  | 质控R2文件匹配模式｜QC R2 file pattern |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-V, --version` | — | version |  |
| `-i, --input` | 必填 |  | 原始FASTQ文件目录路径｜Raw FASTQ files directory path |
| `-g, --genome` | 必填 |  | 参考基因组文件路径｜Reference genome file path |
| `-p, --project-base` | 必填 |  | 项目根目录路径｜Project base directory path |
| `--clean-fastq-dir` | — |  | 清洁FASTQ文件目录路径｜Clean FASTQ files directory path |
| `--mapping-dir` | — |  | 比对结果目录路径｜Mapping results directory path |
| `--gvcf-dir` | — |  | gVCF文件目录路径｜gVCF files directory path |
| `--bam-dir` | — |  | BAM文件目录路径｜BAM files directory path |
| `--joint-dir` | — |  | 联合检测结果目录路径｜Joint calling results directory path |
| `--filter-dir` | — |  | 过滤结果目录路径｜Filtering results directory path |
| `-o, --output-dir` | — |  | 输出目录路径｜Output directory path |
| `--threads-mapping` | `88` | int | 比对线程数｜Number of mapping threads |
| `--threads-gtx` | `88` | int | GTX线程数｜Number of GTX threads |
| `--threads-filter` | `88` | int | 过滤线程数｜Number of filtering threads |
| `--snp-min-dp` | `5` | int | SNP最小深度｜SNP minimum depth |
| `--snp-min-qual` | `30` | int | SNP最小质量｜SNP minimum quality |
| `--indel-min-dp` | `5` | int | InDel最小深度｜InDel minimum depth |
| `--indel-min-qual` | `30` | int | InDel最小质量｜InDel minimum quality |
| `--gatk-threshold` | `4` | int | GATK模式样本数阈值｜GATK sample count threshold |
| `--gtx-single-threshold` | `200` | int | GTX单机模式样本数阈值｜GTX single machine sample count threshold |
| `--gtx-window-size` | `20000000` | int | GTX分块窗口大小｜GTX chunk window size in bp |
| `--gtx-bin` | `~/software/gtx/bin/gtx` |  | GTX可执行文件路径｜GTX executable path |
| `--step` | — | 1/2/3/4/5 | 只运行指定步骤｜Run only specified step |
| `--no-checkpoint` | — | store_true | 禁用断点续传｜Disable checkpoint resume |
| `--dry-run` | — | store_true | 测试模式｜Test mode |
| `-v, --verbose` | `0` | count | 详细输出模式｜Verbose output mode |
| `--quiet` | — | store_true | 静默模式｜Quiet mode |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR/CRITICAL | 日志级别｜Log level |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `-f, --force` | — | store_true | 强制覆盖已存在的结果｜Force overwrite existing results |
| `--skip-qc` | — | store_true | 跳过质控步骤｜Skip QC step |
| `--skip-mapping` | — | store_true | 跳过比对步骤｜Skip mapping step |
| `--read1-pattern-fastp` | `_1.fq.gz` |  | 质控R1文件匹配模式｜QC R1 file pattern |
| `--read2-pattern-fastp` | `_2.fq.gz` |  | 质控R2文件匹配模式｜QC R2 file pattern |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Parabricks（GPU 比对/变异检测，需 NVIDIA GPU）
- GTX（默认 ~/software/gtx/bin/gtx，联合检测；可用 --gtx-bin 覆盖）
- BWA、SAMtools、bcftools、GATK（预检查会校验；自动解析 align 域环境并经 conda run 调用，可用环境变量 BWA_PATH / SAMTOOLS_PATH / BCFTOOLS_PATH / GATK_PATH 覆盖，域环境缺失时回退 PATH 直接调用）
- fastp（质控）
- biopytools（自身，预检查会校验）
- Python 3

## 常见问题 | FAQ { #faq }

**Q1：没有 GPU 能跑吗？**
不能。Parabricks 依赖 NVIDIA GPU；无 GPU 请用 biopytools fastq2vcf-gtx（纯 CPU 版）。

**Q2：样本很少（小于 4）会怎样？**
程序会自动改用 GATK 做变异检测（--gatk-threshold 默认 4），无需手动干预。

**Q3：样本很多（大于等于 200）会怎样？**
联合检测步骤切集群模式，生成操作指南需手动投递，不算失败。

**Q4：断点续传怎么工作？**
默认开启（--no-checkpoint 关闭）。已完成的索引/质控/比对结果存在即跳过；换参数重跑用 --force 强制覆盖。

**Q5：磁盘/内存要求？**
预检查要求项目目录所在盘 >=200GB、内存 >=64GB。
