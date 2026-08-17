# 三代基因组组装流程 | Genome Assembly Pipeline

**多数据类型整合的三代基因组端到端组装流程 | End-to-end third-generation genome assembly pipeline with multi-data integration**

## 功能概述 | Overview

`genomeasm` 是一个面向参考基因组构建的一站式组装流程，以 PacBio HiFi 长读长为核心组装数据，自动整合 ONT、Hi-C、Illumina NGS 等辅助数据。流程会根据输入目录自动检测可用数据类型并选择组装策略。

流程划分为六个阶段：环境检查、数据质量控制（FastQC 等）、hifiasm 主组装、Hi-C 辅助挂载（Juicer + 3D-DNA 或 SALSA2）、质量评估（BUSCO/QUAST）、结果整理报告。Hi-C 挂载策略提供 `complete_juicer`、`standard_3ddna`、`simplified_salsa2` 三种可选方案。适用于植物、动物大基因组的染色体级别组装项目。

## 快速开始 | Quick Start

```bash
# 基本用法：HiFi 组装
biopytools assemble -i raw_data/ -o assembly_results/

# HiFi + Hi-C 染色体级别挂载
biopytools assemble -i data/ -o results/ --hic-strategy complete_juicer

# 指定项目名和基因组大小
biopytools assemble -i input/ -o output/ -n my_genome --genome-size 3g -t 64

# 使用 SALSA2 简化挂载策略
biopytools assemble -i data/ -o results/ --hic-strategy simplified_salsa2
```

注：CLI 命令名为 `assemble`（`biopytools assemble`），对应模块目录为 `genomeasm`。

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input-dir` | 输入目录（包含 HiFi/Hi-C/ONT/NGS 等数据文件）|

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./assembly_output` | 输出目录 |
| `-n, --project-name` | `genome_assembly` | 项目名称 |
| `-t, --threads` | `12` | 线程数 |
| `--genome-size` | `3g` | 基因组大小估计（如 `3g`、`500m`）|
| `--hic-strategy` | `complete_juicer` | Hi-C 挂载策略（complete_juicer / standard_3ddna / simplified_salsa2）|
| `--restriction-enzyme` | `MboI` | 限制性内切酶（MboI/DpnII/HindIII/EcoRI）|
| `--min-contig-size` | `15000` | 进入挂载的最小 contig 长度 |
| `--edit-rounds` | `2` | 3D-DNA 编辑轮数 |
| `--species-type` | `diploid` | 物种倍性（diploid/haploid/polyploid）|
| `--telomere-motif` | `CCCTAA` | 端粒 motif |
| `--purge-level` | `1` | Purge 级别（0-3）|
| `--purge-max` | `80` | Purge 最大覆盖度 |
| `--n-haplotypes` | `2` | 单倍型数 |
| `--skip-fastqc` | 开 | 跳过 FastQC（默认跳过以节省时间）|
| `--min-hifi-coverage` | `30` | HiFi 最低覆盖度阈值 |
| `--min-hic-coverage` | `50` | Hi-C 最低覆盖度阈值 |
| `--min-mapping-rate` | `0.7` | 最低比对率阈值 |
| `--busco-lineage` | `auto` | BUSCO 谱系 |
| `--hifiasm-path` | `hifiasm` | hifiasm 路径 |
| `--juicer-path` | `juicer.sh` | Juicer 脚本路径 |
| `--pipeline-3ddna` | `3d-dna/run-asm-pipeline.sh` | 3D-DNA pipeline 路径 |
| `--salsa2-path` | `run_pipeline.py` | SALSA2 脚本路径 |

（运行 `biopytools assemble -h` 查看完整参数列表）

## 输出 | Output

```
assembly_output/
├── 01_qc/              # 数据质量检查
├── 02_assembly/        # hifiasm 组装结果
├── 03_hic/             # Hi-C 挂载结果（Juicer/3D-DNA/SALSA2）
├── 04_qa/              # BUSCO/QUAST 质量评估
├── final_results/      # 最终染色体级别基因组
└── logs/               # 流程日志
```

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
| `--hifiasm-path` | `hifiasm` |  | Hifiasm程序路径｜Hifiasm program path |
| `--bwa-path` | `bwa` |  | BWA程序路径｜BWA program path |
| `--samtools-path` | `samtools` |  | Samtools程序路径｜Samtools program path |
| `--juicer-path` | `juicer.sh` |  | Juicer脚本路径｜Juicer script path |
| `--pipeline-3ddna` | `3d-dna/run-asm-pipeline.sh` |  | 3D-DNA pipeline路径｜3D-DNA pipeline path |
| `--juicer-tools` | `juicer_tools.jar` |  | Juicer tools JAR路径｜Juicer tools JAR path |
| `--salsa2-path` | `run_pipeline.py` |  | SALSA2脚本路径｜SALSA2 script path |
| `--skip-dependency-check` | — | store_true | 跳过依赖软件检查｜Skip dependency check |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- hifiasm（主组装）
- bwa、samtools（Hi-C 比对）
- Juicer + 3D-DNA（complete_juicer / standard_3ddna 策略）
- SALSA2（simplified_salsa2 策略）
- BUSCO、QUAST（质量评估）
- FastQC（可选）

## 引用 | Citation

- Cheng, H. et al. HiFiasm. *Nature Methods* 18, 170-175 (2021).
- Dudchenko, O. et al. De novo assembly of the Aedes aegypti genome using Hi-C yields chromosome-length scaffolds. *Science* 356, 92-95 (2017).
- Ghurye, J. et al. Integrating Hi-C links with assembly graphs for chromosome-scale assembly. *PLOS Computational Biology* (2019).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
