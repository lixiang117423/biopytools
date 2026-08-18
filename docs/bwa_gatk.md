
# BWA-GATK - 全基因组比对与变异检测 | Whole-genome Alignment & Variant Calling

一句话理解：**把原始测序 reads（FASTQ）一路加工成「哪里和参考基因组不一样」的变异清单（VCF）**，解决「从测序数据到可信变异集合」的标准流程问题。

## 功能概述 | Overview { #overview }

- 完整流程：fastp 质控 → BWA 比对 → BAM 排序/去重 → GATK HaplotypeCaller 生成 GVCF → 合并 GVCF → 联合分型 → 硬/软过滤 → 统计报告
- 自动检测并构建参考基因组索引（BWA + SAMtools + GATK dict）
- 支持多样本（自动从 FASTQ 目录识别样本、支持双端/单端）
- 硬过滤（直接删除不合格位点）与软过滤（只标记不删除）双份输出
- 全程断点续传：每个中间文件存在即跳过，--force-restart 强制重跑

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools bwa-gatk -i fastq_dir/ -g ref.fasta -o output/
```

最小输入：一个放 FASTQ 的目录（含双端 _R1/_R2 或自定义后缀）+ 一个参考基因组 FASTA。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| FASTQ | 测序仪吐出的原始 reads（含序列+每个碱基的质量分） |
| 比对(alignment) | 把每条 read 贴回参考基因组最像的位置 |
| BAM/SAM | 比对结果的文件格式（SAM 是文本，BAM 是压缩二进制） |
| 去重(dedup) | PCR 扩增会造成同一片段被读很多遍，去掉这些「复印本」避免假信号 |
| GVCF | 记录「全基因组每个位点」证据的变异文件，比只记变异位点的 VCF 更全 |
| 联合分型(joint genotyping) | 把多个样本的 GVCF 合起来一起判定变异 |
| 硬过滤 | 不合格的变异直接删除 |
| 软过滤 | 不合格的变异只打标记(FILTER 列)，不删除，留给用户决定 |
| 倍性(ploidy) | 每个位点有几份染色体拷贝（二倍体=2） |

## 输入 | Input { #input }

- **FASTQ 目录**（`-i`）：原始或已清洗的 FASTQ；双端文件用 _R1/_R2（或自定义 --qc-read1-suffix/--qc-read2-suffix）配对，单端用 --qc-single-end
- **参考基因组**（`-g`）：FASTA；首次运行会自动构建 BWA/SAMtools/GATK 索引

格式示例（目录内文件名）：

```text
sample1_R1.fastq.gz    sample1_R2.fastq.gz
sample2_R1.fastq.gz    sample2_R2.fastq.gz
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 两个必填：-i 是放 reads 的目录，-g 是参考基因组。输出目录 -o 不填也行，默认 ./bwa_gatk_output。

### 计算资源 | Computing resources

**通俗理解|In plain words:** --threads（线程数，CLI 界面默认 12）越大比对越快、越吃 CPU；内存方面，GATK 这类 Java 程序吃内存，模块直调可用 -m/--mem-per-thread（默认每线程 10GB）控制。**按机器配置调，一般不用动。**

### 质控参数 | QC parameters

**通俗理解|In plain words:** 默认用 fastp 做质控。--qc-quality-threshold（默认 20）是「碱基平均质量低于多少就剪掉」，--qc-min-length（默认 50）是「剪完后太短的 read 直接丢弃」。输入已经是清洁数据时用 --skip-qc 跳过。**一般不用动。**

### 过滤参数 | Filtering

**通俗理解|In plain words:** 变异过滤按 GATK 最佳实践自动套用硬过滤表达式（QD/FS/SOR/MQ/MQRankSum/ReadPosRankSum 等，SNP 与 INDEL 阈值不同），--snp-qd/--snp-fs/--indel-qd 等可微调。调严=删更多、更保守。**默认即 GATK 推荐值，一般不用动。**

### 运行控制 | Run control

**通俗理解|In plain words:** --force-restart 强制从头重跑（忽略已有文件）；--dry-run 只打印命令不执行；-v 详细日志。断点续传默认开启（中间文件存在即跳过）。

### 工具路径 | Tool paths

**通俗理解|In plain words:** --bwa-path/--gatk-path/--samtools-path/--fastp-path 指定软件路径，默认用系统 PATH 里的 bwa/gatk/samtools/fastp。软件装在 conda 环境里时程序会自动包装。

## 分析流程 | Pipeline { #pipeline }

```text
输入 FASTQ 目录 + 参考基因组
    │
    ▼
步骤1: fastp 质控 (→ 01_qc/)
    │
    ▼
步骤2: 参考基因组索引 (BWA + SAMtools faidx + GATK dict)
    │
    ▼
步骤3: 样本识别 (配对双端/单端)
    │
    ▼
步骤4: BWA mem 比对 → SAM→BAM→排序→MarkDuplicates 去重 (→ 02_alignment/)
    │
    ▼
步骤5: HaplotypeCaller 逐样本生成 GVCF (→ 03_gvcf/)
    │
    ▼
步骤6: CombineGVCFs 合并 → GenotypeGVCFs 联合分型 (→ 04_filter/raw_variants)
    │
    ▼
步骤7: 分离 SNP/INDEL → 硬过滤 + 软过滤 (→ 04_filter/)
    │
    ▼
步骤8: 统计报告 (→ 05_stats/pipeline_statistics.txt)
```

## 输出 | Output { #output }

```text
output/
├── 00_pipeline_info/               # 流程元数据
├── 01_qc/                          # fastp 清洗后的 FASTQ
├── 02_alignment/                   # 各样本 {sample}.dedup.bam (+.bai)
├── 03_gvcf/                        # 各样本 {sample}.g.vcf.gz + combined.g.vcf.gz
├── 04_filter/                      # 变异结果
│   ├── raw_variants.vcf.gz         # 联合分型原始变异
│   ├── raw_snp.vcf.gz / raw_indel.vcf.gz      # 分离后的 SNP/INDEL
│   ├── filtered.hard.SNP.vcf.gz / filtered.hard.INDEL.vcf.gz   # 硬过滤(仅标记)
│   ├── all_samples.filtered.SNP.vcf.gz / all_samples.filtered.INDEL.vcf.gz  # 硬过滤(删除后 PASS)
│   └── filtered.soft.SNP.vcf.gz / filtered.soft.INDEL.vcf.gz   # 软过滤(仅标记)
├── 05_stats/
│   ├── pipeline_statistics.txt     # 比对/变异统计总报告
│   └── {sample}.dedup_metrics.txt  # 各样本去重指标
├── 99_logs/                        # 运行日志
└── temp/                           # 中间文件(运行中)
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. 统计报告 05_stats/pipeline_statistics.txt

**通俗理解|In plain words:** 一页纸看懂本次结果：每个样本比对率、去重率、原始/过滤后变异数、过滤掉了多少。

- 比对率(mapping rate)正常大于 90%；去重率(duplicate rate)过高(如大于 30%)提示文库复杂度低或 PCR 过度
- SNP/INDEL 过滤率异常高(如大于 50%)提示阈值过严或数据质量差

### 2. 硬过滤 vs 软过滤

**通俗理解|In plain words:** 硬过滤结果（all_samples.filtered.*）是「删干净了、直接能用」的最终清单；软过滤结果（filtered.soft.*）保留了所有位点、只在 FILTER 列打标，适合你想自己控制宽松程度时用。

### 3. GVCF vs VCF

**通俗理解|In plain words:** 03_gvcf/ 里是每个样本的全位点证据（中间产物，供联合分型用）；04_filter/ 里才是最终「哪些位置有变异」的清单。

## 参数选择建议 | Parameter Guidance { #guidance }

- **线程与内存**：大基因组/多样本时调高 --threads（如 24~88）；GATK 步骤吃内存，模块直调时把 --mem-per-thread 给足。
- **--skip-qc**：输入已是 clean reads 时跳过，省时间。
- **--force-restart**：换参数或数据更新后，想彻底重跑时用。
- **--dry-run**：第一次跑或换新数据前，先干跑看一遍命令再正式投递。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTQ文件目录（包含原始或清洁的FASTQ文件）｜Input FASTQ directory (containing raw or clean FASTQ files) |
| `--genome, -g` | 必填 |  | 参考基因组FASTA｜Reference genome FASTA file |
| `--output-dir, -o` | `./bwa_gatk_output` | Path | 输出目录｜Output directory |
| `--sample-name` | — |  | 样本名｜Sample name |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--memory` | `300G` |  | 最大内存使用｜Maximum memory usage |
| `--bwa-algorithm` | `mem` | mem/aln | BWA比对算法｜BWA alignment algorithm |
| `--min-seed-length` | `19` | int | 最小种子长度｜Minimum seed length |
| `--band-width` | `100` | int | 带宽｜Band width |
| `--min-base-quality` | `20` | int | 最小碱基质量｜Minimum base quality |
| `--min-mapping-quality` | `20` | int | 最小比对质量｜Minimum mapping quality |
| `--stand-call-conf` | `30.0` | float | 调用置信度阈值｜Stand call confidence threshold |
| `--filter-expression` | — |  | 变异过滤表达式｜Variant filter expression |
| `--min-depth` | `10` | int | 最小测序深度｜Minimum sequencing depth |
| `--max-depth` | — | int | 最大测序深度｜Maximum sequencing depth |
| `--skip-bwa` | — |  | 跳过BWA比对｜Skip BWA alignment step |
| `--skip-gatk` | — |  | 跳过GATK变异检测｜Skip GATK variant calling step |
| `--only-align` | — |  | 仅执行比对｜Only perform alignment step |
| `--only-call` | — |  | 仅执行变异检测｜Only perform variant calling step |
| `--resume` | `True` |  | 启用断点恢复｜Enable resume from checkpoint |
| `--force` | — |  | 强制重新运行所有步骤｜Force rerun all steps |
| `--known-sites` | — | Path | 已知变异位点VCF｜Known variant sites VCF file |
| `--bqsr` | — |  | 启用碱基质量分数重校正｜Enable Base Quality Score Recalibration |
| `--remove-duplicates` | `True` |  | 移除PCR重复｜Remove PCR duplicates |
| `--emit-ref-confidence` | `NONE` | NONE/BP_RESOLUTION/GVCF | 参考置信度模式｜Reference confidence mode |
| `--bwa-path` | — |  | BWA工具路径｜BWA tool path |
| `--gatk-path` | — |  | GATK工具路径｜GATK tool path |
| `--samtools-path` | — |  | Samtools工具路径｜Samtools tool path |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose output mode |
| `--quiet, -q` | — |  | 静默模式｜Quiet mode |
| `--keep-intermediate` | — |  | 保留中间文件｜Keep intermediate files |
| `--skip-qc` | — |  | 跳过质控步骤｜Skip quality control step |
| `--fastp-path` | `fastp` |  | fastp可执行文件路径｜fastp executable path |
| `--qc-threads` | `12` | int | 质控线程数｜QC threads |
| `--qc-quality-threshold` | `20` | int | 质控质量阈值｜QC quality threshold |
| `--qc-min-length` | `50` | int | 质控最小长度｜QC minimum length |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `-i, --input` | 必填 |  | 输入FASTQ文件目录（包含原始或清洁的FASTQ文件）｜Input FASTQ directory (containing raw or clean FASTQ files) |
| `-o, --output-dir` | `./bwa_gatk_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `-p, --ploidy` | `2` | int | 倍性｜Ploidy |
| `-m, --mem-per-thread` | `10` | int | 每线程内存(GB)｜Memory per thread in GB |
| `--intervals` | — |  | 区间限定BED文件｜Intervals BED file for targeted analysis |
| `--snp-qd` | `2.0` | float | SNP QD阈值｜SNP QD threshold |
| `--snp-fs` | `60.0` | float | SNP FS阈值｜SNP FS threshold |
| `--indel-qd` | `2.0` | float | InDel QD阈值｜InDel QD threshold |
| `--force-restart` | — | store_true | 强制从头开始，忽略已有文件｜Force restart, ignore existing files |
| `--dry-run` | — | store_true | 干运行模式，仅显示命令不执行｜Dry-run mode, show commands without execution |
| `-v, --verbose` | — | store_true | 详细日志模式｜Verbose logging mode |
| `--skip-qc` | — | store_true | 跳过质控步骤，直接使用输入数据｜Skip QC step, use input data directly |
| `--fastp-path` | — |  | fastp可执行文件路径｜fastp executable path |
| `--qc-threads` | `12` | int | 质控线程数｜QC threads |
| `--qc-quality-threshold` | `20` | int | 质控质量阈值｜QC quality threshold |
| `--qc-min-length` | `50` | int | 质控最小长度｜QC minimum length |
| `--qc-unqualified-percent` | `40` | int | 不合格碱基百分比阈值｜Unqualified base percentage threshold |
| `--qc-n-base-limit` | `10` | int | N碱基数量限制｜N base count limit |
| `--qc-read1-suffix` | — |  | R1文件后缀｜R1 file suffix |
| `--qc-read2-suffix` | — |  | R2文件后缀｜R2 file suffix |
| `--qc-single-end` | — | store_true | 单末端测序模式｜Single-end sequencing mode |
| `--bwa-path` | — |  | BWA路径｜BWA path |
| `--samtools-path` | — |  | SAMtools路径｜SAMtools path |
| `--gatk-path` | — |  | GATK路径｜GATK path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- BWA（bwa，比对）
- SAMtools（samtools，格式转换/排序/索引/统计）
- GATK（gatk，HaplotypeCaller/MarkDuplicates/VariantFiltration 等）
- fastp（通过 biopytools fastp 子命令调用做质控）
- Python 3

以上软件默认自动解析功能域环境（bwa/samtools/gatk → align，fastp → rna）并经 conda run 调用；可通过 --*-path 参数或环境变量（BWA_PATH / SAMTOOLS_PATH / GATK_PATH / FASTP_PATH）覆盖；域环境缺失时回退 PATH 直接调用。

## 常见问题 | FAQ { #faq }

**Q1：断点续传怎么工作？**
每个中间文件（SAM/BAM/GVCF/VCF）存在即跳过对应步骤。中途失败后直接重跑同一命令，会从断点继续；用 --force-restart 强制从头重跑。

**Q2：为什么提示「参考基因组索引不完整」？**
首次运行会自动构建 BWA（.amb/.ann/.bwt/.pac/.sa）、SAMtools（.fai）、GATK（.dict）索引；若参考基因组文件没有写权限，索引构建会失败。

**Q3：双端 reads 怎么配对？**
默认按 _R1/_R2 后缀配对；如果你的命名不同（如 _1/_2），用 --qc-read1-suffix/--qc-read2-suffix 指定；单端数据加 --qc-single-end。

**Q4：CLI 界面上有些参数（如 --sample-name、--memory、--bwa-algorithm、--min-depth、--only-align 等）是什么？**
这些是 CLI 包装器（click）暴露的高级选项，其中一部分当前未接入后端 argparse（python -m biopytools.bwa_gatk），直接使用可能报「无法识别的参数」。请以参数速查表中「模块直调参数」实际支持的选项为准。

**Q5：硬过滤和软过滤有什么区别？**
硬过滤把不合格位点删除（最终结果）；软过滤只标记不删除，保留全部位点供你自己判断。两者都会输出，按需选用。
