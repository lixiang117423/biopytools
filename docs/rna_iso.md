# rna-iso | 全长转录本分析(IsoSeq3+IsoQuant)|Full-length transcript analysis (IsoSeq3+IsoQuant)

把三代测序的长读长 RNA 数据（PacBio 或 nanopore）变成一份「全长转录本清单」——每条转录本的完整序列、它属于哪个基因、有哪些可变剪接形式、表达量多高，一条命令拿到手。

## 功能概述 | Overview

- 三种输入形态自动识别，不需要记命令差异：PacBio 原始 `subreads.bam`、已生成好的 `ccs.bam`（HiFi）、ONT 的 `fastq/fasta`
- 双引擎可选：
  - **IsoQuant 4.0**（默认）：把 reads 比对到参考基因组后重建转录本模型，输出 GTF/GFF3 注释文件 + 基因/转录本表达量表 + 可变剪接分析
  - **IsoSeq3 26.2**（仅 PacBio）：不依赖参考基因组，从 reads 直接聚类出全长转录本序列（polished 高质量 fasta）
- PacBio 数据自动完成 ccs（生成 HiFi reads）→ refine（去引物、去 polyA、去嵌合）前处理
- 多个测序 run 文件可一起传入（同一样品），自动合并分析
- 断点续传：中断重跑自动跳过已完成的步骤

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools rna-iso -i sample.fq.gz --data-type ont -g genome.fa -o output/
```

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗解释<br>Plain words |
|---|---|
| 全长转录本<br>Full-length transcript | 一条 mRNA 被完整地从头读到尾——像拿到整封信，而不是撕碎的纸片。三代长读长能做到，二代短读长只能拼碎片 |
| isoform（异构体） | 同一个基因的「不同版本信」——剪接方式不同，蛋白质可能就不同。找全 isoform 是本模块的核心价值 |
| FLNC | Full-Length Non-Chimeric：既有头（5'引物）又有尾（3'引物+polyA）且中间不断裂的读长——真正可信的全长 reads |
| subreads / CCS(HiFi) | PacBio 原始下机的环化共识读子叫 subreads；把同一分子的多圈读子算平均得到的高精度读长叫 CCS/HiFi |
| ccs | subreads → HiFi 的加工步骤，非常吃 CPU（一个 SMRT cell 约几十分钟到数小时） |
| refine | 把 HiFi reads 砍头去尾（去引物/polyA）、剔除嵌合读长，得到 FLNC |
| 聚类/cluster2 | 把来自同一转录本的众多 FLNC 归成一堆，每堆生成一条代表序列=转录本（IsoSeq3 引擎的步骤） |
| 参考引导 vs de novo | IsoQuant 需要「地图」（参考基因组）把 reads 对上去重建模型；IsoSeq3 不用地图，直接按序列相似度聚类。有参考时优先 IsoQuant |

## 输入 | Input

| 形态<br>Form | 文件示例<br>Example | 模块行为<br>Behavior |
|---|---|---|
| PacBio subreads | `m54086_170204.subreads.bam`（同名 `.pbi` 建议同放，缺失仅警告） | ccs → refine → 引擎 |
| PacBio CCS | `sample.ccs.bam` / `sample.hifi_reads.bam` | refine → 引擎 |
| reads 文件 | `sample.fq.gz`、`sample.fasta`（ONT 或已有 HiFi） | 直接进引擎；**必须加 `--data-type pacbio` 或 `ont`** |

多个文件 = 同一样品的多个测序 run（如多个 SMRT cell），重复 `-i` 传入即可；**`-i` 也可以直接传目录**——自动按文件名排序收集目录内的 fasta/fastq（±gz），忽略子目录。**不同形态不要混传**（会报错）。

引物 fasta 格式（仅 PacBio 链用到；不传则用内置 Clontech SMARTer）：

```
>primer_5p
AAGCAGTGGTATCAACGCAGAGTACATGGG
>primer_3p
GTACTCTGCGTTGATACCACTGCTT
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --reads` | 必填 | Path | 输入文件或目录(可多个;目录自动收集fasta/fastq±gz):subreads.bam/ccs.bam/fasta/fastq(±gz)｜Input file(s) or dir (dir auto-collects fasta/fastq±gz) |
| `--data-type` | — | pacbio/ont | reads文件(fasta/fastq)时必填;BAM输入自动嗅探｜Required for fasta/fastq inputs; auto-sniffed for BAM |
| `-g, --reference` | — | Path | 参考基因组FASTA(isoquant引擎必填)｜Reference genome FASTA (required for isoquant engine) |
| `--genedb` | — | Path | 参考注释GTF/GFF(可选,提高精确度)｜Reference annotation GTF/GFF (optional) |
| `--engine` | `isoquant` | isoquant/isoseq3/both | 转录本重建引擎｜Transcript engine |
| `--primers` | — | Path | 引物fasta(默认内置Clontech SMARTer)｜Primer fasta (default: built-in Clontech SMARTer) |
| `--min-passes` | `1` | int | ccs最小pass数(Iso-Seq官方推荐1)｜ccs min passes |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-p, --prefix` | `rna_sample` |  | 样本前缀｜Sample prefix |
| `-o, --output-dir` | `./rna_iso_output` | Path | 输出目录｜Output directory |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --reads` | 必填 |  | 输入文件或目录(可多个;目录自动收集fasta/fastq±gz):subreads.bam/ccs.bam/fasta/fastq(±gz)｜Input file(s) or dir (dir auto-collects fasta/fastq±gz) |
| `--data-type` | — | pacbio/ont | reads文件(fasta/fastq)时必填;BAM输入自动嗅探｜Required for fasta/fastq inputs; auto-sniffed for BAM |
| `-g, --reference` | — |  | 参考基因组FASTA(isoquant引擎必填)｜Reference genome FASTA (required for isoquant engine) |
| `--genedb` | — |  | 参考注释GTF/GFF(可选,提高精确度)｜Reference annotation GTF/GFF (optional, improves precision) |
| `--engine` | `isoquant` | isoquant/isoseq3/both | 转录本重建引擎(默认isoquant)｜Transcript engine (default: isoquant) |
| `--primers` | — |  | 引物fasta(默认内置Clontech SMARTer)｜Primer fasta (default: built-in Clontech SMARTer) |
| `--min-passes` | `1` | int | ccs最小pass数(默认1,Iso-Seq官方推荐)｜ccs min passes (default: 1) |
| `-t, --threads` | `12` | int | 线程数(默认12)｜Number of threads (default: 12) |
| `-p, --prefix` | `rna_sample` |  | 样本前缀(默认rna_sample)｜Sample prefix |
| `-o, --output-dir` | `./rna_iso_output` |  | 输出目录(默认./rna_iso_output)｜Output directory |

<!-- END PARAMS:auto -->

## 分析流程 | Pipeline

```
PacBio subreads.bam ──> 01 ccs(HiFi reads) ──> 02 refine(FLNC) ──┬──> 03 IsoQuant(GTF+表达量)
PacBio ccs.bam ──────────────────────────────────────────────────┤
ONT/HiFi fasta|fastq ────────────────────────────────────────────┘
                                                                 └──> 04 cluster2(de novo转录本)
                                                                      (仅PacBio, --engine isoseq3|both)
```

## 输出 | Output

```
output_dir/
├── 00_pipeline_info/software_versions.yml   # 工具版本与运行参数
├── 01_ccs/<movie>.ccs.bam                   # PacBio subreads输入时:HiFi reads
├── 02_refine/<movie>.flnc.bam + .flnc.fa    # PacBio输入时:全长非嵌合reads
├── 03_isoquant/<prefix>/                    # IsoQuant标准输出树:
│   ├── <prefix>.transcript_models.gtf       #   最终转录本模型(喂下游注释/IGV)
│   ├── <prefix>.extended_annotation.gtf     #   原注释+新发现转录本合并(给了genedb时)
│   ├── <prefix>.transcript_counts.tsv       #   转录本表达计数
│   ├── <prefix>.gene_counts.tsv             #   基因表达计数
│   ├── <prefix>.*_tpm.tsv                   #   TPM表达量(基因/转录本/新发现)
│   └── <prefix>.polyA_prediction.tsv        #   polyA位点预测
├── 04_isoseq3/<prefix>.transcripts.bam/.fasta  # isoseq3引擎:de novo转录本
└── 99_logs/rna_iso.log                      # 模块总日志
```

（条件性步骤目录按输入形态创建；上面全列出仅为示意。）

## 结果解读 | Interpreting Results

**IsoQuant（默认引擎）**
- `<prefix>.transcript_models.gtf`：最终产物（e2e 实测核心产物）。每个转录本一条记录，含外显子坐标。好结果判据：转录本数在物种预期量级（如植物几万条）、绝大多数落在已知基因位点（若给了 genedb，另有 `<prefix>.extended_annotation.gtf` 合并新发现转录本）
- `gene_counts.tsv` / `transcript_counts.tsv`：表达量表，TPM 与 raw counts 双列；下游差异分析用 raw counts
- 日志里 `Reads assigned to genes` 比例高（>60%）说明数据与参考匹配良好；偏低要检查物种是否对应、reads 是否被污染

**isoseq3 引擎**
- `<prefix>.transcripts.fasta`：de novo 转录本序列。条数与 IsoQuant 同量级为正常；数量级偏少通常意味着引物不匹配（检查 `--primers`）
- 适合作为 SQANTI 等下游工具的输入，或无参考物种的主结果

## 参数选择建议 | Parameter Guidance

- **`--engine`**：有参考基因组用默认 `isoquant`（GTF 直接喂注释流程）；想同时拿 de novo 序列对比用 `both`；无参考只能 `isoseq3`（仅 PacBio）
- **`--data-type`**：只在输入是 fasta/fastq 时需要——模块无法从序列内容可靠区分 PacBio 和 ONT；BAM 输入时自动嗅探，传了也会被忽略
- **`--genedb`**：**给了更好**——把已有注释（哪怕是不完整的）交给 IsoQuant，新转录本发现和表达定量都更准
- **`--primers`**：**一般不用动**。内置的是 Clontech SMARTer 通用引物（绝大多数 Iso-Seq 文库用它）。只有自建库用别的引物时才需要换
- **`--min-passes`**：**不用动**（Iso-Seq 官方推荐 1）。调高更保守但丢数据多
- **`--threads`**：ccs 步骤极吃 CPU，建议能给多少给多少（如 64）；IsoQuant/cluster2 受益递减

## 依赖 | Dependencies

- conda env `isoseq_v.4.0.0`：ccs (pbccs) 6.4.0、isoseq3 26.2.0、IsoQuant 4.0.0、samtools 1.24
- 路径可用环境变量覆盖：`CCS_PATH` / `ISOSEQ3_PATH` / `ISOQUANT_PATH` / `SAMTOOLS_PATH`

## 常见问题 | FAQ

**Q1: isoseq3 引擎什么时候不可用？**
输入是 fasta/fastq（无 BAM）或 ONT 数据时。isoseq3 的 refine/cluster2 只吃 PacBio BAM。此时 `--engine isoseq3` 会直接报错，`--engine both` 会报错（不是静默跳过——避免你以为跑了其实没跑）。

**Q2: 旧文档里的 lima 和 polish 去哪了？**
isoseq3 26.2 移除了 polish（HiFi reads 无需 subreads 再打磨），去引物并入 `refine` 子命令（不再需要单独的 lima 步骤）。本模块直接使用新版 CLI。

**Q3: ccs 要跑多久？**
非常吃 CPU：一个 SMRT cell 官方实测约 50 分钟 × 70 核。模块的 ccs 步骤把 `--threads` 全部用上，建议提交大线程作业。

**Q4: 输入混了 subreads.bam 和 fastq 为什么报错？**
同一条命令只接受同一种形态（它们属于不同加工阶段）。分开两次跑。

**Q5: `.pbi` 文件缺失会怎样？**
仅 WARNING，ccs 可无索引运行（速度下降）。建议从下机目录一起拷贝。

**Q6: IsoQuant 报 reads 没有 polyA？**
IsoQuant 期望 reads 保留 polyA 尾（用于定位 3' 端）。如果上游做过去 polyA 处理，重建质量会下降；raw reads 直接输入即可。

**Q7: 想要多个样品各自出结果？**
一条命令一个样品（`-p` 前缀区分），多样品分别运行后汇总。`-i` 传多个文件是「同一样品的多个 run」，会合并成一个结果。
