# EviAnn 基因组注释 | EviAnn Genome Annotation

一句话理解：**把「基因组里哪些地方是基因、每个基因长什么样」这件事，用 RNA-seq 和蛋白质证据自动推断出来**——输入基因组序列，再喂给它转录组（二代/三代）和/或同源蛋白，输出基因结构注释（GFF）和蛋白/转录本序列。

## 功能概述 | Overview { #overview }

- 基于证据的真核生物基因组注释，证据可来自二代 RNA-seq、三代全长转录组、转录本 FASTA、蛋白 FASTA 任意组合
- 调用 EviAnn 主脚本（eviann.sh），一条命令完成基因结构预测、转录本与蛋白序列输出
- 支持 lncRNA 注释（`--lncrna-tpm` 阈值）、部分 CDS、现有 CDS GFF 校正、功能注释（`--functional`）
- 输出直接落在基因组文件所在目录，用基因组文件名做前缀，无需指定输出目录
- 纯 Python 包装器，外部计算全部交给 EviAnn 完成

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools eviann -g genome.fa --long-reads longreads.fq.gz -p proteins.fa -t 12
```

最小要求：一个基因组 FASTA（`-g`）+ 至少一种证据（RNA-seq 或转录本 `-e`）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| 基因组注释 | 在一条长序列上标出「哪里是基因、外显子和内含子怎么切」，像给一长串密码做批注 |
| 基因结构 | 一个基因由多段外显子（编码段）被内含子（非编码段）隔开组成，注释就是还原这个拼接关系 |
| 二代测序 RNA-seq | 把 RNA 打成小碎片测序，量大价廉，适合拼出外显子边界 |
| 三代测序（Iso-Seq/ONT） | 直接读出整条全长转录本，不用拼，能看清完整拼接结构 |
| 同源蛋白 | 别的物种里已知的蛋白序列，用来「按图索骥」帮新物种找基因 |
| CDS | 真正编码蛋白的那段序列，外显子里去掉两端非翻译区后的部分 |
| lncRNA | 不编码蛋白的 RNA，但有调控功能；TPM 是衡量它表达量的单位 |
| 倍性（ploidy） | 一个细胞里有多少套染色体，人类是 2（二倍体），注释时告诉软件以免误判 |

## 输入 | Input { #input }

### 基因组 FASTA（必需，`-g`）

标准 FASTA 格式。输出文件会写到这个文件所在的目录，文件名以基因组文件名为前缀。

### 证据数据（至少一种）

- **二代转录组** `--short-reads`：FASTQ 文件或目录（识别 `.fq`/`.fastq`/`.fq.gz`/`.fastq.gz`），自动标记为 fastq 类型
- **三代转录组** `--long-reads`：FASTQ 或 FASTA 文件/目录（额外识别 `.fa`/`.fa.gz`/`.fasta`/`.fasta.gz`），自动标记为 isoseq 类型
- **转录本 FASTA** `-e/--transcripts`：已有的转录本序列，作为直接证据
- **蛋白质 FASTA** `-p/--proteins`：同源蛋白，辅助预测基因结构（可选）

二代/三代转录组输入会自动生成 EviAnn 描述文件（`rnaseq_inputs.txt`，写到基因组目录），把每个文件的绝对路径和类型登记进去。

## 参数说明 | Parameters { #parameters }

### 证据输入 | Evidence inputs

**通俗理解|In plain words:** 告诉 EviAnn「用哪些证据来猜基因」。给得越全、质量越高，注释越准；但至少要有一种 RNA-seq 或转录本证据，否则程序拒绝运行。蛋白（`-p`）和 UniProt（`-s`）是辅助，可给可不给。

相关参数：`-g/--genome`（必需）、`--short-reads`、`--long-reads`、`-e/--transcripts`、`-p/--proteins`、`-s/--uniprot`、`-c/--cds-gff`、`--extra-gff`、`--mito-contigs`。

### 运行与注释选项 | Run & annotation options

**通俗理解|In plain words:** `-t` 线程数越大越快（默认 12 一般够用）；`-d` 倍性默认 2 绝大多数真核不用动；`--max-intron` 是允许的最大内含子长度，不填就交给 EviAnn 自动判断，只在物种内含子特别长/特别短时才需要手调；`--lncrna-tpm` 越低保留的 lncRNA 越多、越高越严格，默认 1.0。

相关参数：`-t/--threads`（默认 12）、`-d/--ploidy`（默认 2）、`-m/--max-intron`（默认自动）、`--lncrna-tpm`（默认 1.0）。

### 开关选项 | Flags

**通俗理解|In plain words:** `--partial` 打开后会把不完整（缺头少尾）的 CDS 也算进去，注释更全但可能掺入假阳性；`--functional` 打开后额外跑功能注释，更慢；`--debug`/`--verbose` 是排错用的，出问题时再开。

相关参数：`--partial`、`--functional`、`--debug`、`--verbose`。

## 分析流程 | Pipeline { #pipeline }

```text
基因组 FASTA + 证据数据
    │
    ▼
生成 RNA-seq 描述文件(二代/三代自动登记)
    │
    ▼
构建 EviAnn 命令(conda run 调用 eviann.sh)
    │
    ▼
EviAnn 执行:基因结构预测 + (可选)功能注释
    │
    ▼
输出 .gff / .proteins.fasta / .transcripts.fasta
```

## 输出 | Output { #output }

输出写在基因组文件所在目录（无独立输出目录参数），以基因组文件名为前缀：

```text
<基因组所在目录>/
├── genome.pseudo_label.gff      # 基因结构注释(主结果)
├── genome.proteins.fasta        # 预测蛋白序列
├── genome.transcripts.fasta     # 预测转录本序列
├── genome.eviann.log            # 运行日志
└── rnaseq_inputs.txt            # 自动生成的 RNA-seq 描述文件(给了二代/三代时)
```

`--functional` 打开时，EviAnn 会额外输出功能注释相关文件（随 EviAnn 版本而定）。

## 结果解读 | Interpreting Results { #interpreting-results }

- `*.pseudo_label.gff`：核心产物，含每个预测基因的坐标、外显子/内含子结构，可加载到 IGV 或 JBrowse 可视化，也可喂给下游工具
- `*.proteins.fasta`：预测的蛋白序列，条数即预测的蛋白编码基因数，可用 BUSCO 等评估完整性
- `*.transcripts.fasta`：预测的转录本序列
- 日志里会打印每个输出文件是否存在，缺失会以 warning 提示，据此判断哪一步没产出

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 只有二代 RNA-seq：`--short-reads`；只有三代全长：`--long-reads`；两者都有就都加上，互补
- 有现成转录本/蛋白：优先 `-e` 或 `-p` 直接喂，比从头组装更省事
- 默认参数即可覆盖多数物种；内含子极长的物种（如哺乳动物）可调大 `--max-intron`
- 需要后续做功能富集：加 `--functional`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 | Path | 基因组FASTA文件｜Genome FASTA file (required) |
| `--short-reads` | — | Path | 二代转录组数据（文件或目录）｜Short-read RNA-seq data (file or directory) |
| `--long-reads` | — | Path | 三代转录组数据（文件或目录）｜Long-read RNA-seq data (file or directory) |
| `-e, --transcripts` | — | Path | 转录本FASTA文件｜Transcripts FASTA file |
| `-p, --proteins` | — | Path | 蛋白质FASTA文件｜Proteins FASTA file |
| `-s, --uniprot` | — | Path | UniProt-SwissProt FASTA｜UniProt-SwissProt FASTA |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-m, --max-intron` | — | int | 最大内含子长度｜Maximum intron length (default: auto) |
| `-d, --ploidy` | `2` | int | 基因组倍性｜Genome ploidy |
| `-c, --cds-gff` | — | Path | 现有CDS的GFF文件｜GFF file with existing CDS |
| `--lncrna-tpm` | `1.0` | float | lncRNA最小TPM｜Minimum TPM for lncRNA |
| `--partial` | `False` |  | 包含部分CDS｜Include partial CDS |
| `--functional` | `False` |  | 执行功能注释｜Perform functional annotation |
| `--mito-contigs` | — | Path | 线粒体contig列表｜File with mitochondrial contigs |
| `--extra-gff` | — | Path | 额外的GFF特征｜Extra features from GFF |
| `--debug` | `False` |  | 调试模式｜Debug mode |
| `--verbose` | `False` |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--short-reads` | — |  | 二代转录组数据（文件或目录）｜Short-read RNA-seq data (file or directory) |
| `--long-reads` | — |  | 三代转录组数据（文件或目录）｜Long-read RNA-seq data (file or directory) |
| `-e, --transcripts` | — |  | 转录本FASTA文件｜Transcripts FASTA file |
| `-p, --proteins` | — |  | 蛋白质FASTA文件｜Proteins FASTA file |
| `-s, --uniprot` | — |  | UniProt-SwissProt FASTA｜UniProt-SwissProt FASTA |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-m, --max-intron` | — | int | 最大内含子长度｜Maximum intron length |
| `-d, --ploidy` | `2` | int | 基因组倍性｜Genome ploidy |
| `--lncrna-tpm` | `1.0` | float | lncRNA最小TPM｜Minimum TPM for lncRNA |
| `--partial` | — | store_true | 包含部分CDS｜Include partial CDS |
| `--functional` | — | store_true | 执行功能注释｜Perform functional annotation |
| `--debug` | — | store_true | 调试模式｜Debug mode |
| `--verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- EviAnn 软件：`~/miniforge3/envs/eviann_v.2.0.5`（conda 环境 `eviann_v.2.0.5`，通过 `bin/eviann.sh` 调用）
- Python 3（仅用于包装，无额外 Python 依赖）

## 常见问题 | FAQ { #faq }

**Q1：为什么没有 `-o` 输出目录参数？**
输出固定写在基因组文件所在目录，用基因组文件名做前缀，这是为了跟 EviAnn 的原生行为保持一致。想换个地方，就把基因组文件复制或软链到目标目录再跑。

**Q2：报「必须提供 RNA-seq 数据或转录本数据」？**
`-g` 之外必须至少给 `--short-reads`、`--long-reads` 或 `-e/--transcripts` 之一，只给 `-p` 蛋白不够。

**Q3：会断点续传吗？**
包装器本身不管理断点续传，每次重跑都重新调用 eviann.sh。若中途失败，重新运行同一条命令即可（EviAnn 内部对已完成的子步骤有自己的缓存）。

**Q4：输出文件没生成怎么办？**
先看 `<genome>.eviann.log` 日志里的报错；常见原因是 EviAnn 环境路径不对或证据文件格式不被识别。
