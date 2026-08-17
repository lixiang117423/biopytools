# QIIME2 - 微生物组多样性分析 | QIIME2 Microbiome Diversity Analysis

一句话理解：**把双端 16S/ITS 扩增子测序数据自动跑完 QIIME2 全流程**——从原始 reads 得到 ASV/OTU 丰度表、物种注释、系统发育树、多样性指标，并导出可直接分析的表格。

## 功能概述 | Overview { #overview }

- 自动识别目录里的双端 FASTQ 样本（默认后缀 _1.clean.fq.gz / _2.clean.fq.gz）
- 全流程：导入→切引物(cutadapt)→去噪(DADA2 得 ASV)→可选 OTU 聚类(vsearch)→分类注释(sklearn)→系统发育(mafft+fasttree)→多样性+抽平→导出
- 支持 16S（默认）与 ITS 两种扩增子；ITS 自动跳过系统发育
- 分类器三级获取：用户指定 → 缓存复用 → 自动训练（SILVA/UNITE）
- 抽样深度自动取每样本 reads 的第 10 百分位，也可手动指定
- 逐步骤断点续传；--force 清空重跑

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools qiime2 -i raw_reads/ -o qiime2_output
```

最小输入：一个装双端 FASTQ 的目录。首次运行会自动训练分类器（约 30–60 分钟，之后缓存复用）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 扩增子(amplicon) | 用 PCR 扩增出的标记基因片段（16S 或 ITS） |
| ASV | 「精确序列型」，按单碱基差异分辨，去噪后的最小单元 |
| OTU | 「操作分类单元」，按相似度（默认 97%）把相似序列归成一簇 |
| DADA2 | 去噪算法，纠正测序错误并产出 ASV |
| vsearch | 聚类算法，把序列按相似度聚成 OTU |
| 去噪(denoise) | 把测序错误纠正掉，保留真实序列 |
| 抽平(rarefaction) | 把所有样本抽到相同 reads 数，消除测序深度差异后再比较 |
| alpha/beta 多样性 | alpha=一个样本内部的多样性；beta=样本之间的差异 |
| .qza / .qzv | QIIME2 的中间产物（.qza 数据 / .qzv 可交互可视化） |
| 分类器(classifier) | 预训练的「物种鉴定模型」，用于给序列贴分类标签 |

## 输入 | Input { #input }

- **双端 FASTQ 目录**（`-i`）：R1/R2 按后缀配对，默认 `_1.clean.fq.gz` / `_2.clean.fq.gz`（可用 `--r1-suffix` 与 `--r2-suffix` 修改）。
- **样品元数据 TSV**（`--metadata-file`，可选）：用于多样性分组的表；不提供则自动生成一个「全部归为 all 组」的最小元数据。

示例目录：

```text
raw_reads/
├── S1_1.clean.fq.gz
├── S1_2.clean.fq.gz
├── S2_1.clean.fq.gz
└── S2_2.clean.fq.gz
```

## 参数说明 | Parameters { #parameters }

### 必需与输出 | Required & output

**通俗理解|In plain words:** `-i` 输入目录、`-o` 输出目录（默认 ./qiime2_output）、`-t` 线程数（默认 12）。

### 分析类型 | Analysis type

**通俗理解|In plain words:** `--amplicon` 选 16s（默认）或 its；`--method` 选 asv（DADA2，默认，精度高）或 otu（vsearch 聚类，经典做法）。**ITS 请务必设 `--amplicon its`，且程序会自动跳过系统发育。**

### 引物与截断 | Primers & truncation

**通俗理解|In plain words:** `--fwd-primer` / `--rev-primer` 默认是 V3-V4 的 341F/806R；`--skip-cutadapt` 用于已经去掉引物的数据。`--trunc-len-f` / `--trunc-len-r` 决定保留的读长（默认 0=不截断），`--trim-left-f` / `--trim-left-r` 从左端裁掉若干碱基。**截断长度按你的读长定**（见参数选择建议），不设则整条用。

### 抽样深度 | Sampling depth

**通俗理解|In plain words:** `--sampling-depth` 决定多样性分析时每个样本抽多少 reads。默认 0=自动取第 10 百分位（稳妥，兼顾保留样本数和深度）；样本深度差异大时可手动指定一个大家都够得着的值。

### 聚类与分类 | Clustering & classification

**通俗理解|In plain words:** `--perc-identity` 是 OTU 聚类相似度（默认 0.97，仅 otu 模式）；`--confidence` 是分类置信度门槛（默认 0.7）；`--min-length` / `--max-length` 是训练分类器时从参考库提取目标区域的长度限制。**一般不用动。**

### 分类器与路径 | Classifier & paths

**通俗理解|In plain words:** `--classifier` 传现成的预训练分类器 .qza（省去训练）；`--database-dir` 指定 SILVA/UNITE 原始库目录（默认 `~/database/qiime2`）；`--qiime-path` 指定 qiime 可执行文件；`--classifier-cache-dir` 指定训练好的分类器缓存目录。**首次没训练过会自动训练并缓存，之后自动复用。**

### 命名与执行 | Naming & execution

**通俗理解|In plain words:** `--r1-suffix` / `--r2-suffix` 改配对后缀；`--skip-phylogeny` 跳过建树（ITS 自动）；`-f` 覆盖已有输出重跑；`-v` 详细日志。

## 分析流程 | Pipeline { #pipeline }

```text
检测双端样本 → 生成 manifest 与元数据
    │
    ▼
步骤1: 导入双端 FASTQ(demux.qza + 质量汇总)
    │
    ▼
步骤2: 切引物(cutadapt，可选，--skip-cutadapt 跳过)
    │
    ▼
步骤3: 去噪(DADA2 → ASV 表/代表序列/去噪统计)
    │
    ▼
步骤4: OTU 聚类(可选，--method otu 时用 vsearch)
    │
    ▼
步骤5: 分类注释(classify-sklearn)
    │
    ▼
步骤6: 系统发育(mafft+fasttree，ITS/--skip-phylogeny 跳过)
    │
    ▼
步骤7: 核心多样性 + alpha 稀释曲线(先定抽样深度)
    │
    ▼
步骤8: 导出(FASTA/BIOM/TSV/分类表/合并表)
```

## 输出 | Output { #output }

```text
qiime2_output/
├── 00_pipeline_info/                 # software_versions.yml、sample_manifest.tsv、sample_metadata.tsv
├── 01_import/                        # demux.qza、demux.qzv(质量图)
├── 02_trim/                          # demux_trimmed.qza
├── 03_denoise/                       # table.qza、rep_seqs.qza、denoising_stats.qza、table.qzv、rep_seqs.qzv
├── 04_otu/                           # otu_table.qza、otu_seqs.qza(仅 otu 模式)
├── 05_taxonomy/                      # taxonomy.qza、taxonomy.qzv
├── 06_phylogeny/                     # aligned/masked 序列、unrooted/rooted 树(.qza)
├── 07_diversity/                     # core_metrics(_phylogenetic)/、alpha_rarefaction.qzv
├── 08_export/
│   ├── asv_sequences.fasta           # 代表序列(otu 模式为 otu_sequences.fasta)
│   ├── feature_table.biom            # 丰度表 BIOM
│   ├── feature_table.tsv             # 丰度表 TSV
│   ├── taxonomy.tsv                  # 分类表
│   └── feature_table_w_taxonomy.tsv  # 合并分类的丰度表(核心交付物)
└── 99_logs/
    └── qiime2_pipeline.log
```

## 结果解读 | Interpreting Results { #interpreting-results }

### 1. 08_export/feature_table_w_taxonomy.tsv（核心交付物）

**通俗理解|In plain words:** 行是 ASV/OTU，列是样本，数值是丰度，最后一列 taxonomy 是分类注释。**这是下游分析最常用的表**，直接拿去画堆叠柱状图、做差异分析。

### 2. asv_sequences.fasta 与 feature_table.tsv

**通俗理解|In plain words:** 代表序列（可再做 BLAST/建树）与纯丰度表（无分类列）。feature_table.biom 是 QIIME2 生态可用的机器格式。

### 3. taxonomy.tsv

**通俗理解|In plain words:** 每个 ASV/OTU 的完整分类（门到种，分号分隔）。置信度低于阈值的位置会用 `Unassigned` 之类占位。

### 4. 多样性结果（07_diversity/）

**通俗理解|In plain words:** `core_metrics(_phylogenetic)/` 里是各 alpha/beta 指标与距离矩阵、PCoA；`alpha_rarefaction.qzv` 是 alpha 稀释曲线。用 `qiime tools view` 打开 .qzv 可交互查看。16S 用 phylogenetic 版（含 UniFrac/Faith），ITS 用非系统发育版。

### 5. 好坏判据

- **demux.qzv 质量图**：用于决定截断长度（质量在末端明显下跌就该截）。
- **denoising_stats.qza**：看每样本 reads 保留率，丢失过多提示截断/引物设置问题。
- **alpha 稀释曲线是否饱和**：抽平深度足够时曲线应趋于平台；未饱和说明深度不够或多样性高。
- **feature_table_w_taxonomy.tsv 中 Unassigned 比例**：过高可能是分类器与数据不匹配（引物区域不对）或置信度设太高。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- **16S V3-V4 双端 250bp**：`--trunc-len-f 220 --trunc-len-r 200`；**300bp** 用 `--trunc-len-f 230 --trunc-len-r 210`（依据 demux 质量图微调）。
- **ITS 数据**：`--amplicon its`（建议同时 `--skip-cutadapt` 或提供 ITS 专用引物）。
- **已去引物的数据**：加 `--skip-cutadapt`。
- **要 OTU 而非 ASV**：`--method otu`（默认 97% 相似度）。
- **想固定抽样深度**：`--sampling-depth <值>`；否则用默认自动（第 10 百分位）。
- **已有训练好的分类器**：`--classifier xx.qza` 跳过训练，快很多。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-dir` | 必填 |  | 双端FASTQ输入目录｜Input directory of paired-end FASTQ |
| `-o, --output-dir` | `./qiime2_output` | Path | 输出目录｜Output directory |
| `--amplicon` | `16s` | 16s/its | 扩增子类型｜Amplicon type |
| `--method` | `asv` | asv/otu | 聚类方法ASV(DADA2)或OTU(vsearch)｜Method: ASV or OTU |
| `--fwd-primer` | `CCTACGGGNGGCWGCAG` |  | 正向引物序列(IUPAC)｜Forward primer |
| `--rev-primer` | `GACTACHVGGGTATCTAATCC` |  | 反向引物序列(IUPAC)｜Reverse primer |
| `--trunc-len-f` | `0` | int | R1截断长度(0=不截断)｜R1 truncation length (0=none) |
| `--trunc-len-r` | `0` | int | R2截断长度(0=不截断)｜R2 truncation length (0=none) |
| `--trim-left-f` | `0` | int | R1左侧裁剪｜R1 trim left |
| `--trim-left-r` | `0` | int | R2左侧裁剪｜R2 trim left |
| `--sampling-depth` | `0` | int | 抽平深度(0=自动)｜Rarefaction depth (0=auto) |
| `--perc-identity` | `0.97` | float | OTU聚类相似度｜OTU identity |
| `--confidence` | `0.7` | float | 分类置信度｜Classification confidence |
| `--min-length` | `50` | int | extract-reads最小长度｜extract-reads min length |
| `--max-length` | `0` | int | extract-reads最大长度(0=不限)｜extract-reads max length (0=none) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--validate-level` | `min` | min/max | tools import校验级别｜import validate level |
| `--classifier` | — |  | 预训练分类器(.qza),省略则自动训练｜Pre-trained classifier (.qza) |
| `--database-dir` | — |  | 原始参考库目录(SILVA/UNITE)｜Raw reference DB directory |
| `--qiime-path` | — |  | qiime可执行文件路径｜qiime executable path |
| `--classifier-cache-dir` | — |  | 分类器缓存目录｜Classifier cache directory |
| `--r1-suffix` | `_1.clean.fq.gz` |  | R1文件后缀｜R1 file suffix |
| `--r2-suffix` | `_2.clean.fq.gz` |  | R2文件后缀｜R2 file suffix |
| `--skip-cutadapt` | — |  | 跳过引物切除｜Skip primer trimming |
| `--skip-phylogeny` | — |  | 跳过系统发育建树(ITS自动跳过)｜Skip phylogeny (auto for ITS) |
| `--metadata-file` | — |  | 样品元数据TSV(可选)｜Sample metadata TSV (optional) |
| `-f, --force` | — |  | 覆盖已有输出｜Overwrite existing outputs |
| `-v, --verbose` | — |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-dir` | 必填 |  | 双端FASTQ输入目录｜Input directory of paired-end FASTQ |
| `-o, --output-dir` | `./qiime2_output` |  | 输出目录｜Output directory (default: ./qiime2_output) |
| `--amplicon` | `16s` | 16s/its | 扩增子类型｜Amplicon type (default: 16s) |
| `--method` | `asv` | asv/otu | 聚类方法｜Method: ASV(DADA2) or OTU(vsearch) (default: asv) |
| `--fwd-primer` | `CCTACGGGNGGCWGCAG` |  | 正向引物序列(IUPAC)｜Forward primer (default: 341F) |
| `--rev-primer` | `GACTACHVGGGTATCTAATCC` |  | 反向引物序列(IUPAC)｜Reverse primer (default: 806R) |
| `--trunc-len-f` | `0` | int | R1截断长度(0=不截断)｜R1 truncation length (0=none) |
| `--trunc-len-r` | `0` | int | R2截断长度(0=不截断)｜R2 truncation length (0=none) |
| `--trim-left-f` | `0` | int | R1左侧裁剪｜R1 trim left |
| `--trim-left-r` | `0` | int | R2左侧裁剪｜R2 trim left |
| `--sampling-depth` | `0` | int | 抽平深度(0=自动取第10百分位)｜Rarefaction depth (0=auto) |
| `--perc-identity` | `0.97` | float | OTU聚类相似度｜OTU identity (default: 0.97) |
| `--confidence` | `0.7` | float | classify-sklearn置信度｜classification confidence (default: 0.7) |
| `--min-length` | `50` | int | extract-reads最小长度｜extract-reads min length (default: 50) |
| `--max-length` | `0` | int | extract-reads最大长度(0=不限)｜extract-reads max length (0=none) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--validate-level` | `min` | min/max | tools import校验级别｜import validate level (default: min) |
| `--classifier` | — |  | 预训练分类器(.qza),省略则自动训练｜Pre-trained classifier (.qza), auto-train if omitted |
| `--database-dir` | — |  | 原始参考库目录(SILVA/UNITE)｜Raw reference DB directory (default: ~/database/qiime2) |
| `--qiime-path` | — |  | qiime可执行文件路径｜qiime executable path (default: ~/miniforge3/envs/qiime_v.2024.10.1/bin/qiime) |
| `--classifier-cache-dir` | — |  | 分类器缓存目录｜Classifier cache directory (default: <db>/classifier_cache) |
| `--r1-suffix` | `_1.clean.fq.gz` |  | R1文件后缀｜R1 file suffix (default: _1.clean.fq.gz) |
| `--r2-suffix` | `_2.clean.fq.gz` |  | R2文件后缀｜R2 file suffix (default: _2.clean.fq.gz) |
| `--skip-cutadapt` | — | store_true | 跳过引物切除(数据已去引物)｜Skip primer trimming |
| `--skip-phylogeny` | — | store_true | 跳过系统发育建树(ITS自动跳过)｜Skip phylogeny (auto for ITS) |
| `--metadata-file` | — |  | 样品元数据TSV(可选)｜Sample metadata TSV (optional) |
| `-f, --force` | — | store_true | 覆盖已有输出｜Overwrite existing outputs |
| `-v, --verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- QIIME2（conda 环境 `qiime_v.2024.10.1`，默认 `~/miniforge3/envs/qiime_v.2024.10.1/bin/qiime`）
- 参考数据库：16S 用 SILVA（默认 `~/database/qiime2/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz`）；ITS 用 UNITE（默认 `~/database/qiime2/` 下的 UNITE tgz）
- biom（QIIME2 环境内，用于 biom↔tsv 转换）

## 常见问题 | FAQ { #faq }

**Q1：支持断点续传吗？**
支持，逐步骤。每步以关键输出（如 `03_denoise/table.qza`、`05_taxonomy/taxonomy.qza`）为标记，存在即跳过。用 `-f/--force` 清空各步骤目录重跑。

**Q2：第一次跑为什么特别慢？**
因为要自动训练分类器（16S 约 30–60 分钟）。训练结果会缓存到分类器缓存目录，之后自动复用；也可 `--classifier xx.qza` 直接给现成的。

**Q3：ITS 能跑系统发育树吗？**
程序对 ITS 自动跳过系统发育（--skip-phylogeny 置真），多样性用非系统发育版核心指标。16S 默认会建树。

**Q4：抽样深度怎么定？**
默认自动取每样本 reads 总数的第 10 百分位，兼顾深度与样本保留。想自己定就 `--sampling-depth N`（N 不能超过最浅样本的总 reads）。

**Q5：切引物为什么用反向互补？**
cutadapt 的 --p-front-r 需要反向引物的反向互补序列，程序内部已用 reverse_complement 自动转换，用户只需按 5'→3' 给正常引物序列即可。

**Q6：没给 --metadata-file 会怎样？**
程序自动生成一个最小元数据（所有样本归为 all 组），满足 core-metrics 的元数据要求；要做分组比较时再传自己的 TSV。