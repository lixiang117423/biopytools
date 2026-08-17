
# GATK Joint Genotyping - 多样本联合分型 | GATK Joint Genotyping

一句话理解：**把每个样本各自生成的 GVCF（记录全基因组每个位点的证据）合并起来，一次性「共同判定」所有样本的变异**，解决「单样本变异检测不一致、多样本如何统一」的问题。

## 功能概述 | Overview { #overview }

- 自动识别输入目录中的 GVCF（优先）或 VCF 文件
- 完整流程：GenomicsDB 导入 → GenotypeGVCFs 联合分型 → 提取/过滤 SNP → 提取/过滤 INDEL → 合并
- 自动从参考基因组生成全基因组区间文件（无需手动提供）
- 按 GATK 最佳实践分别过滤 SNP 与 INDEL（QD/FS/MQ/MQRankSum/ReadPosRankSum/SOR）
- 输出原始 VCF、分类型过滤 VCF、合并 VCF 与分析汇总报告

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools gatk-joint -i gvcf_folder/ -g ref.fasta -o results/
```

最小输入：一个含 GVCF/VCF 的目录 + 参考基因组 FASTA。输出目录默认 ./joint_genotyping_output。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| GVCF | 记录「全基因组每个位点证据」的文件，比只记变异位点的 VCF 更全 |
| 联合分型(joint genotyping) | 把所有样本的 GVCF 合起来，统一判定变异与基因型 |
| GenomicsDB | GATK 的一种高效数据仓库，先导入 GVCF 再统一分型 |
| sample map | 样本名与文件路径的对照表 |
| intervals(区间) | 分析哪些染色体/区段 |

## 输入 | Input { #input }

- **GVCF/VCF 目录**（`-i`）：优先识别 *.g.vcf.gz / *.g.vcf / *.gvcf.gz / *.gvcf；没有 GVCF 时识别 *.vcf.gz / *.vcf
- **参考基因组**（`-g`）：FASTA（.fasta/.fa），需有对应的 .dict/.fai（缺的会自动创建）

## 参数说明 | Parameters { #parameters }

### 必需与输出 | Required & output

**通俗理解|In plain words:** -i 是放 GVCF/VCF 的目录，-g 参考基因组，-o 输出目录。两个必填，一般用默认输出目录即可。

### 计算资源 | Computing resources

**通俗理解|In plain words:** --threads（默认 12）控制 GenomicsDB 读取线程；--memory（默认 100g）是给 GATK Java 程序的内存。样本多、基因组大时调高。**按机器配置调。**

### 区间设置 | Intervals

**通俗理解|In plain words:** --intervals（-L）限定只分析某些染色体/区段。不填时程序自动从参考基因组提取全基因组区间。**一般不用动，想只跑部分染色体时再用。**

### SNP 过滤 | SNP filtering

**通俗理解|In plain words:** QD>=2、FS<=60、MQ>=40、MQRankSum>=-12.5、ReadPosRankSum>=-8、SOR<=3，达不到即标记为低质量。**GATK 最佳实践默认值，一般不用动。**

### INDEL 过滤 | INDEL filtering

**通俗理解|In plain words:** INDEL 用更宽松标准（QD>=2、FS<=200、ReadPosRankSum>=-20、SOR<=10）。**默认即合理。**

### 工具路径 | Tool paths

**通俗理解|In plain words:** --gatk-path（默认 gatk）、--bcftools-path（默认 bcftools）指定软件路径，conda 环境中的软件会自动处理。**一般不用动。**

## 分析流程 | Pipeline { #pipeline }

```text
GVCF/VCF 目录 + 参考基因组
    │
    ▼
步骤1: 检测文件类型(GVCF优先) → 生成 sample_map.txt
    │
    ▼
步骤2: 生成全基因组区间(自动) → GenomicsDBImport 导入
    │
    ▼
步骤3: GenotypeGVCFs 联合分型 (→ joint_genotyping_raw.vcf.gz)
    │
    ▼
步骤4: 提取并过滤 SNP (SelectVariants + VariantFiltration)
    │
    ▼
步骤5: 提取并过滤 INDEL
    │
    ▼
步骤6: 合并 SNP+INDEL (bcftools concat)
    │
    ▼
汇总报告 (analysis_summary.txt)
```

## 输出 | Output { #output }

```text
output/
├── sample_map.txt                              # 样本名→文件路径映射
├── intervals.list                              # 自动生成的全基因组区间
├── genomicsdb_workspace/                       # GenomicsDB 数据仓库(中间产物)
├── joint_genotyping_raw.vcf.gz                 # 联合分型原始 VCF
├── joint_genotyping_snps_raw.vcf.gz            # 分离出的 SNP
├── joint_genotyping_indels_raw.vcf.gz          # 分离出的 INDEL
├── joint_genotyping_snps_filtered.vcf.gz       # 过滤后的 SNP
├── joint_genotyping_indels_filtered.vcf.gz     # 过滤后的 INDEL
├── joint_genotyping_merged_filtered.vcf.gz     # 合并后的最终 VCF (+.tbi)
└── analysis_summary.txt                        # 分析汇总报告
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. joint_genotyping_merged_filtered.vcf.gz（最终结果）

**通俗理解|In plain words:** 合并了所有样本、过滤后的 SNP+INDEL 最终 VCF，直接用于下游群体分析。

### 2. 分类型 VCF

**通俗理解|In plain words:** *_snps_filtered.vcf.gz 和 *_indels_filtered.vcf.gz 是分开的 SNP/INDEL 结果，需要单独用某类变异时取对应文件。

### 3. analysis_summary.txt（汇总报告）

**通俗理解|In plain words:** 记录输入信息、过滤参数、各阶段变异数量（原始/过滤后 SNP/INDEL/合并），一页纸看懂结果。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规多样本**：全部默认即可。
- **样本很多**：调高 --threads（GenomicsDB 读取）与 --memory（如 100g~300g）。
- **只分析部分染色体**：用 --intervals chr1（或区间文件）。
- **输入是普通 VCF（非 GVCF）**：程序会自动识别并按 VCF 处理，无需额外参数。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入目录(包含VCF/GVCF文件)｜Input directory (containing VCF/GVCF files) |
| `--genome, -g` | 必填 |  | 参考基因组文件(.fasta/.fa)｜Reference genome file (.fasta/.fa) |
| `--output-dir, -o` | `./joint_genotyping_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--memory, -m` | `100g` |  | Java内存设置｜Java memory setting |
| `--intervals, -L` | — |  | 分析区间(染色体或区间文件)｜Analysis intervals (chromosome or interval file) |
| `--snp-qd` | `2.0` | float | SNP QD阈值｜SNP QD threshold |
| `--snp-fs` | `60.0` | float | SNP FS阈值｜SNP FS threshold |
| `--snp-mq` | `40.0` | float | SNP MQ阈值｜SNP MQ threshold |
| `--snp-mqrs` | `-12.5` | float | SNP MQRankSum阈值｜SNP MQRankSum threshold |
| `--snp-rprs` | `-8.0` | float | SNP ReadPosRankSum阈值｜SNP ReadPosRankSum threshold |
| `--snp-sor` | `3.0` | float | SNP SOR阈值｜SNP SOR threshold |
| `--indel-qd` | `2.0` | float | INDEL QD阈值｜INDEL QD threshold |
| `--indel-fs` | `200.0` | float | INDEL FS阈值｜INDEL FS threshold |
| `--indel-rprs` | `-20.0` | float | INDEL ReadPosRankSum阈值｜INDEL ReadPosRankSum threshold |
| `--indel-sor` | `10.0` | float | INDEL SOR阈值｜INDEL SOR threshold |
| `--gatk-path` | `gatk` |  | GATK软件路径｜GATK software path |
| `--bcftools-path` | `bcftools` |  | BCFtools软件路径｜BCFtools software path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入目录(包含VCF/GVCF文件)｜Input directory (containing VCF/GVCF files) |
| `-g, --genome` | 必填 |  | 参考基因组文件｜Reference genome file (.fasta/.fa) |
| `-o, --output-dir` | `./joint_genotyping_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `-m, --memory` | `100g` |  | Java内存设置｜Java memory setting |
| `-L, --intervals` | — |  | 分析区间(染色体或区间文件)｜Analysis intervals (chromosome or interval file) |
| `--snp-qd` | `2.0` | float | SNP QD阈值｜SNP QD threshold |
| `--snp-fs` | `60.0` | float | SNP FS阈值｜SNP FS threshold |
| `--snp-mq` | `40.0` | float | SNP MQ阈值｜SNP MQ threshold |
| `--snp-mqrs` | `-12.5` | float | SNP MQRankSum阈值｜SNP MQRankSum threshold |
| `--snp-rprs` | `-8.0` | float | SNP ReadPosRankSum阈值｜SNP ReadPosRankSum threshold |
| `--snp-sor` | `3.0` | float | SNP SOR阈值｜SNP SOR threshold |
| `--indel-qd` | `2.0` | float | INDEL QD阈值｜INDEL QD threshold |
| `--indel-fs` | `200.0` | float | INDEL FS阈值｜INDEL FS threshold |
| `--indel-rprs` | `-20.0` | float | INDEL ReadPosRankSum阈值｜INDEL ReadPosRankSum threshold |
| `--indel-sor` | `10.0` | float | INDEL SOR阈值｜INDEL SOR threshold |
| `--gatk-path` | `gatk` |  | GATK软件路径｜GATK software path |
| `--bcftools-path` | `bcftools` |  | BCFtools软件路径｜BCFtools software path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- GATK（gatk，GenomicsDBImport/GenotypeGVCFs/SelectVariants/VariantFiltration/CreateSequenceDictionary）
- bcftools（bcftools，合并/索引/统计）
- samtools（自动创建 .fai 索引时用到）
- Python 3

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
不支持完整断点续传。每次运行都会重新导入 GenomicsDB（会先删除旧的 genomicsdb_workspace/ 再重建），并从联合分型开始重新跑。

**Q2：为什么自动生成 intervals.list？**
GenomicsDBImport 需要按区间导入；程序会自动从参考基因组的 .dict/.fai 提取染色体，生成 intervals.list，无需手动准备。

**Q3：参考基因组没有 .dict 怎么办？**
程序会自动调用 GATK CreateSequenceDictionary 创建；若失败会退而用 samtools faidx 的 .fai。

**Q4：GVCF 和 VCF 混在一起会怎样？**
程序优先识别 GVCF（*.g.vcf.gz 等），找到 GVCF 就按 GVCF 流程走；没有 GVCF 才按普通 VCF 处理。

**Q5：输出文件名为什么都是 joint_genotyping 开头？**
模块内部 base_name 固定为 joint_genotyping，所以输出前缀统一，无法通过 CLI 修改。
