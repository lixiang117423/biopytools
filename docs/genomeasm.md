# 基因组组装流程 | Genome Assembly Pipeline

一句话理解：**把一堆三代测序数据（HiFi 必需，可加 Hi-C/ONT/NGS）自动跑成一版染色体级别的基因组——从质控、组装、挂载到质量评估全流程自动完成**。

## 功能概述 | Overview

- 以 PacBio HiFi 为核心组装数据，自动整合 ONT、Hi-C、Illumina NGS 等辅助数据
- 自动检测输入目录里的数据类型，并据此选择组装策略
- 六阶段流程：环境检查 → 数据质控 → hifiasm 组装 → Hi-C 挂载 → 质量评估 → 结果整理
- 三种 Hi-C 挂载策略：complete_juicer（Juicer+3D-DNA）、standard_3ddna、simplified_salsa2（SALSA2）
- 质量评估输出 N50、BUSCO 完整性、gap 等指标，并生成最终报告与 README

## 快速开始 | Quick Start

```bash
biopytools genomeasm -i raw_data/ -o results/
```

把 HiFi reads（文件名含 hifi/pacbio/ccs 等关键字）放进输入目录，运行即可；有 Hi-C 就再放 Hi-C 的 R1/R2 文件。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| HiFi | PacBio 的高准确长读段，本流程的「主力组装数据」 |
| ONT | Oxford Nanopore 长读段，超长但错误率偏高，可辅助提升连续性 |
| NGS | Illumina 短读段，短而准，常用于精修/纠错 |
| Hi-C | 能测「染色体空间上谁挨着谁」的数据，用来把 contig 排成染色体（挂载） |
| 组装(assembly) | 把读段拼回基因组序列的过程 |
| 挂载(scaffolding) | 借助 Hi-C 把零散 contig 排序、定向、连成整条染色体 |
| contig / scaffold / 染色体 | 三级递进：contig 是连续无空洞的片段，scaffold 是排好序但中间可能有 N 的片段，染色体是最完整形态 |
| N50 | 把所有序列按长度从大到小排，累加到总长一半时那条序列的长度；越大越连续 |
| BUSCO | 一套「单拷贝保守基因」数据库，用来测基因组的完整性（缺了多少基因） |
| 单倍型/倍性 | 二倍体(diploid)有两套染色体，单倍型(haplotype)指其中一套；hifiasm 可分别组装两套 |

## 输入 | Input

输入是一个**目录**，程序按文件名关键字自动识别数据类型（HiFi 是必需的）：

| 数据类型 | 文件名关键字 | 是否必需 |
|----------|-------------|----------|
| HiFi | hifi / pacbio / ccs | 必需（检测不到会报错） |
| Hi-C | hic R1/R2 | 可选 |
| ONT | ont / nanopore / long | 可选 |
| NGS | ngs / illumina / short | 可选 |

示例目录结构：

```text
raw_data/
├── sample.hifi.fastq.gz        # HiFi（必需）
├── sample_hic_R1.fastq.gz      # Hi-C R1（可选）
└── sample_hic_R2.fastq.gz      # Hi-C R2（可选）
```

## 参数说明 | Parameters

### 基本参数 | Basic

**通俗理解|In plain words:** -i 是放测序数据的目录；-o 是结果输出目录；-n 是项目名（会出现在输出文件名里）；-t 是线程数。**注意 -t 在 CLI 帮助里显示默认 12，但实际不指定时模块用的是 88（见 FAQ）。**

### 组装参数 | Assembly

**通俗理解|In plain words:** --genome-size 是预估基因组大小（用于覆盖度估算和 QC，如 3g=30 亿、500m=5 亿），尽量给准；--species-type 是倍性，二倍体默认 diploid；--purge-level/--purge-max 控制「去掉冗余/重复序列」的力度（越大删得越狠，默认 1 较保守）；--telomere-motif 是端粒序列（动物默认 CCCTAA，植物换成 TTTAGGG）。**不熟悉就用默认。**

### Hi-C 参数 | Hi-C

**通俗理解|In plain words:** 只有目录里有 Hi-C 数据才用到。--hic-strategy 选挂载方案：complete_juicer 最完整质量最高、最慢；simplified_salsa2 最简化最快；standard_3ddna 折中。--restriction-enzyme 是建库用的酶（默认 MboI），必须和实际建库一致。

### 质控参数 | Quality control

**通俗理解|In plain words:** --skip-fastqc 默认跳过 FastQC 省时间；--min-hifi-coverage/--min-hic-coverage 是「覆盖度低于多少算不合格」的告警线，不达线只告警不中断；--min-mapping-rate 是 Hi-C 比对率下限；--busco-lineage 是 BUSCO 谱系（auto 自动选）。**一般不用动。**

### 工具路径 | Tool paths

**通俗理解|In plain words:** 各软件的路径，默认自动解析功能域环境，缺失时回退 PATH 里的命令名。只有软件装在不常见位置时才需要显式指定。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 六步走：先查环境 → 查数据质量 → 用 hifiasm 组装 → 有 Hi-C 就挂载成染色体 → 评估质量 → 整理成报告。

```text
输入目录（自动检测数据类型）
    │
    ▼
Phase 1: 环境检查 + 依赖验证 + 输入文件校验
    │
    ▼
Phase 2: 数据质量控制（覆盖度/读长/FastQC）
    │
    ▼
Phase 3: hifiasm 组装 → GFA 转 FASTA → 组装统计
    │
    ▼
Phase 4: Hi-C 挂载（有 Hi-C 数据时；失败不终止流程）
    ├─ complete_juicer: Juicer + 3D-DNA
    ├─ standard_3ddna: BWA + pairtools + 3D-DNA
    └─ simplified_salsa2: BWA + SALSA2
    │
    ▼
Phase 5: 质量评估（N50/L50、BUSCO、gap）
    │
    ▼
Phase 6: 结果整理 + 最终报告 + README + JSON 摘要
```

## 输出 | Output

```text
results/
├── logs/                         # 运行日志 + 命令记录
├── qc/                           # 数据质控报告
├── assembly/                     # hifiasm 组装产物（GFA/FASTA + 组装统计）
├── hic_processing/               # Hi-C 挂载中间产物
├── quality_assessment/           # 综合质量报告 + BUSCO 结果
├── temp/                         # 临时文件
└── results/                      # 最终整理结果（重点看这里）
    ├── FINAL_ASSEMBLY_REPORT.txt # 最终组装报告
    ├── assembly_results_summary.json  # 机器可读的结果摘要
    ├── README.md                 # 结果说明
    ├── final_assemblies/         # 最终组装序列 FASTA（+ .fai 索引）
    ├── scaffolds/                # 染色体挂载后的 FASTA
    ├── for_juicebox/             # Juicebox 手动校正文件（.hic/.assembly）
    └── statistics/               # 统计文件
```

## 结果解读 | Interpreting Results

### 1. 最终报告（results/FINAL_ASSEMBLY_REPORT.txt）

**通俗理解|In plain words:** 一张「组装成绩单」，先看它。

- 项目信息：组装策略、Hi-C 策略、数据类型、预估大小；
- 质量摘要：总体评级、最佳组装 N50、最佳挂载 N50、最高 BUSCO 完整性；
- 主要输出文件清单 + 改进建议。

### 2. 组装序列（results/final_assemblies/）

**通俗理解|In plain words:** hifiasm 组出来的最终序列。

- primary：主要组装（优先用）；alternate：备选组装；haplotype_1/2：分相的两套单倍型；
- 二倍体物种优先取 primary 版本。

### 3. 挂载序列（results/scaffolds/）

**通俗理解|In plain words:** 有 Hi-C 时，contig 被排成染色体的版本，命名形如 项目名_primary_xxx_scaffolded.fa。

### 4. 质量指标

**通俗理解|In plain words:** 看三个数判断好不好。

- N50 越大越好（染色体级通常几十 Mb）；
- BUSCO 完整性越接近 100% 越好（>95% 优秀，<90% 需关注）；
- 大小比例（组装大小 / 预估大小）在 0.9-1.1 之间算合理。

## 参数选择建议 | Parameter Guidance

- **--hic-strategy**：有 Hi-C 且追求高质量用默认 complete_juicer；赶时间或 Hi-C 质量一般用 simplified_salsa2
- **--genome-size**：尽量给准（如 3g、500m），影响覆盖度估算和 QC 判定
- **--species-type**：正常二倍体用默认 diploid；单倍体（如雄性膜翅目）用 haploid；多倍体用 polyploid
- **-t/--threads**：核多就调大，hifiasm/BWA 都能明显提速
- **--telomere-motif**：植物换成 TTTAGGG，动物默认 CCCTAA

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-dir` | 必填 | Path | () |
| `-o, --output-dir` | `./assembly_output` | Path |  |
| `-n, --project-name` | `genome_assembly` |  |  |
| `-t, --threads` | `12` | int |  |
| `--hic-strategy` | `complete_juicer` | complete_juicer/standard_3ddna/simplified_salsa2 | Hi-C |
| `--restriction-enzyme` | `MboI` | MboI/DpnII/HindIII/EcoRI |  |
| `--min-contig-size` | `15000` | int | contig |
| `--edit-rounds` | `2` | int | 3D-DNA |
| `--genome-size` | `3g` |  | (e.g., 3g, 500m) |
| `--species-type` | `diploid` | diploid/haploid/polyploid |  |
| `--telomere-motif` | `CCCTAA` |  | motif |
| `--purge-level` | `1` | 0/1/2/3 | Purging |
| `--purge-max` | `80` | int | Purging |
| `--similarity-threshold` | `0.75` | float |  |
| `--n-haplotypes` | `2` | int |  |
| `--skip-fastqc` | `True` |  | FastQC () |
| `--min-hifi-coverage` | `30` | int | HiFi |
| `--min-hic-coverage` | `50` | int | Hi-C |
| `--min-mapping-rate` | `0.7` | float |  |
| `--busco-lineage` | `auto` |  | BUSCO |
| `--hifiasm-path` | `hifiasm` |  | Hifiasm |
| `--bwa-path` | `bwa` |  | BWA |
| `--samtools-path` | `samtools` |  | Samtools |
| `--juicer-path` | `juicer.sh` |  | Juicer |
| `--pipeline-3ddna` | `3d-dna/run-asm-pipeline.sh` |  | 3D-DNA pipeline |
| `--juicer-tools` | `juicer_tools.jar` |  | Juicer tools JAR |
| `--salsa2-path` | `run_pipeline.py` |  | SALSA2 |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-dir` | 必填 |  | 输入数据目录 (自动检测文件类型)｜Input data directory (auto-detect file types) |
| `-o, --output-dir` | `./assembly_output` |  | 输出目录｜Output directory |
| `-n, --project-name` | `genome_assembly` |  | 项目名称｜Project name |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `--hic-strategy` | `complete_juicer` | complete_juicer/standard_3ddna/simplified_salsa2 | Hi-C处理策略｜Hi-C processing strategy |
| `--restriction-enzyme` | `MboI` | MboI/DpnII/HindIII/EcoRI | 限制性酶类型｜Restriction enzyme type |
| `--min-contig-size` | `15000` | int | 最小contig大小阈值｜Minimum contig size threshold |
| `--edit-rounds` | `2` | int | 3D-DNA编辑轮数｜3D-DNA editing rounds |
| `--genome-size` | `3g` |  | 预估基因组大小｜Estimated genome size (e.g., 3g, 500m) |
| `--species-type` | `diploid` | diploid/haploid/polyploid | 物种倍性｜Species ploidy |
| `--telomere-motif` | `CCCTAA` |  | 端粒序列motif｜Telomere sequence motif |
| `--purge-level` | `1` | 0/1/2/3 | Purging级别｜Purging level |
| `--purge-max` | `80` | int | Purging覆盖度上限｜Purging coverage upper limit |
| `--similarity-threshold` | `0.75` | float | 相似度阈值｜Similarity threshold |
| `--n-haplotypes` | `2` | int | 单倍型数量｜Number of haplotypes |
| `--skip-fastqc` | `True` | store_true | 跳过FastQC质量检查 (默认跳过，节省时间)｜Skip FastQC quality check (default: skip to save time) |
| `--run-fastqc` | — | store_true | 运行FastQC质量检查｜Run FastQC quality check |
| `--min-hifi-coverage` | `30` | int | 最小HiFi覆盖度｜Minimum HiFi coverage |
| `--min-hic-coverage` | `50` | int | 最小Hi-C覆盖度｜Minimum Hi-C coverage |
| `--min-mapping-rate` | `0.7` | float | 最小映射率｜Minimum mapping rate |
| `--busco-lineage` | `auto` |  | BUSCO谱系数据库｜BUSCO lineage database |
| `--hifiasm-path` | — |  | Hifiasm程序路径｜Hifiasm program path |
| `--bwa-path` | — |  | BWA程序路径｜BWA program path |
| `--samtools-path` | — |  | Samtools程序路径｜Samtools program path |
| `--juicer-path` | — |  | Juicer脚本路径｜Juicer script path |
| `--pipeline-3ddna` | — |  | 3D-DNA pipeline路径｜3D-DNA pipeline path |
| `--juicer-tools` | — |  | Juicer tools JAR路径｜Juicer tools JAR path |
| `--salsa2-path` | — |  | SALSA2脚本路径｜SALSA2 script path |
| `--skip-dependency-check` | — | store_true | 跳过依赖软件检查｜Skip dependency check |

<!-- END PARAMS:auto -->
## 依赖 | Dependencies

- 基础：hifiasm（asm 域）、bwa、samtools（align 域）、seqkit（misc 域）、fastqc
- Hi-C（按策略）：complete_juicer 需 juicer + juicer_tools + 3d-dna；standard_3ddna 需 bwa + pairtools（hic 域）+ 3d-dna + juicer_tools；simplified_salsa2 需 bwa + samtools + SALSA2（run_pipeline.py）
- 质量评估：BUSCO（busco 域，可选，装不上会自动跳过 BUSCO 步骤）
- Java（用于 juicer_tools）

conda 工具自动解析功能域环境并经 conda run 调用（环境变量 HIFIASM_PATH / BWA_PATH / SAMTOOLS_PATH / SEQKIT_PATH / PAIRTOOLS_PATH / BUSCO_PATH 或 --*-path 参数覆盖）；域环境缺失时回退 PATH 直接调用。fastqc、juicer、3d-dna、juicer_tools、SALSA2 无对应功能域环境，保持旧默认值（PATH 直接调用，环境变量 FASTQC_PATH / JUICER_PATH / PIPELINE_3DDNA_PATH / JUICER_TOOLS_PATH / SALSA2_PATH 覆盖）。

## 常见问题 | FAQ

**Q1：命令到底是 biopytools genomeasm 还是 biopytools assemble？**
是 **biopytools genomeasm**。CLI 注册表里命令名为 genomeasm（模块目录也叫 genomeasm）；内部 click 命令对象名叫 assemble，只影响 --help 的显示，不影响调用。

**Q2：为什么报「未检测到 HiFi 数据」？**
HiFi 是必需的。确认输入目录里 HiFi 文件名含 hifi/pacbio/ccs 关键字（不区分大小写的通配匹配），否则程序无法识别。

**Q3：-t 默认到底是 12 还是 88？**
CLI 帮助显示 12（click 层默认值），但 click 层不把默认值透传给模块，模块 argparse 实际默认 88，所以不指定 -t 时实际用 88 线程。想要指定就显式传 -t。

**Q4：Hi-C 挂载失败了会怎样？**
不会终止整个流程。挂载失败只跳过挂载结果，继续做质量评估和结果整理（组装结果仍会保留）。

**Q5：支持断点续传吗？**
不支持。每次运行都会从头走完整六阶段流程。中断后重跑会重新开始。

**Q6：为什么 --skip-fastqc 默认就是跳过？**
为节省时间默认跳过 FastQC（模块里 skip_fastqc 默认 True）。想跑 FastQC 需在模块直调时用 --run-fastqc；click 包装层未暴露 --run-fastqc。

