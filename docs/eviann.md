# EviAnn 基因组注释 | EviAnn Genome Annotation

一句话理解：**把「基因组里哪些地方是基因、每个基因长什么样」这件事，用 RNA-seq 和蛋白质证据自动推断出来**——输入基因组序列，再喂给它转录组（二代/三代）和/或同源蛋白，输出基因结构注释（GFF）和蛋白/转录本序列。本模块负责把杂乱的原始测序文件自动整理成 EviAnn 需要的格式并跑完全程。

## 功能概述 | Overview { #overview }

- 基于证据的真核生物基因组注释，证据可来自二代 RNA-seq、三代全长转录组、转录本 FASTA、蛋白 FASTA 任意组合
- **自动整理输入**：给一批文件或一个目录，自动识别二代/三代、自动配对双端、自动合并多样本多次测序、自动把同一样本的二代+三代配成混合（mix）模式
- 三种输入模式：`--rnaseq-data` 全自动 / `--sample-sheet` 精确清单 / `-r` 原生描述文件透传
- 独立输出目录（by-step 结构），不污染基因组所在目录；最终结果齐全时重跑自动跳过（断点续传友好）
- 支持 lncRNA 注释（`--lncrna-tpm` 阈值）、部分 CDS、现有 CDS GFF 校正、功能注释（`-f`）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools eviann -g genome.fa --rnaseq-data ./rna_data/ -p proteins.fa -t 12 -o out/
```

最小要求：一个基因组 FASTA（`-g`）+ 输出目录（`-o`）+ 至少一种证据（`--rnaseq-data` 或 `-e`）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| 基因组注释 | 在一条长序列上标出「哪里是基因、外显子和内含子怎么切」，像给一长串密码做批注 |
| 基因结构 | 一个基因由多段外显子（编码段）被内含子（非编码段）隔开组成，注释就是还原这个拼接关系 |
| 二代测序 RNA-seq | 把 RNA 打成小碎片测序，量大价廉，适合拼出外显子边界 |
| 三代测序（Iso-Seq/ONT） | 直接读出整条全长转录本，不用拼，能看清完整拼接结构 |
| mix 模式 | 同一个样本既有二代又有三代时，两者合并喂给 EviAnn 联合组装，效果比分开用好 |
| 同源蛋白 | 别的物种里已知的蛋白序列，用来「按图索骥」帮新物种找基因 |
| CDS | 真正编码蛋白的那段序列，外显子里去掉两端非翻译区后的部分 |
| lncRNA | 不编码蛋白的 RNA，但有调控功能；TPM 是衡量它表达量的单位 |
| 倍性（ploidy） | 一个细胞里有多少套染色体，人类是 2（二倍体），注释时告诉软件以免误判 |

## 输入 | Input { #input }

### 基因组 FASTA（必需，`-g`）

标准 FASTA 格式。模块会在输出目录 `work/` 下为它建软链（保持原文件名），EviAnn 的全部中间文件和最终结果都以基因组文件名为前缀落在 `work/`。

### 转录组数据（三种模式任选其一，或用 `-e` 转录本）

**模式A：`--rnaseq-data` 全自动（推荐）**

接受单个文件、多个文件或目录（逗号分隔可混合传多个目录/文件）。模块自动完成：

1. **类型识别**：文件名含 `isoseq/ccs/hifi/pacbio/nanopore/ont/dorado` 等关键词或为 FASTA → 三代；BAM 按关键词区分二代/三代
2. **命名惯例识别**：本实验室默认命名——三代 `*.clean.fq.gz`、二代 `*_1.clean.fq.gz` / `*_2.clean.fq.gz`；无关键词时按此规则分类
3. **双端自动配对**：支持 `_R1/_R2`、`_1/_2`、`.R1.`、`read1/read2`、`_R1_001` 等常见后缀
4. **多样本合并**：同一样本多次测序的文件（如多个三代 cell）自动二进制拼接为一个文件（gzip 可直接拼接）
5. **mix 自动配对**：二代样本名与三代样本名做前缀匹配，配上的组成混合行（效果最好），配不上的各自独立
6. 自动生成的样本清单写到 `01_inputs/sample_sheet.tsv`，**可人工查看修改后用模式B重跑**

**模式B：`--sample-sheet` 精确清单**

TSV 格式，每行一个样本，`#` 开头为注释行，列间用 Tab 分隔：

```text
# sample_id	r1	r2	long_reads	tag(可省)
S1	/x/S1_1.clean.fq.gz	/x/S1_2.clean.fq.gz	/x/S1.clean.fq.gz
S2	/x/S2_1.clean.fq.gz	/x/S2_2.clean.fq.gz
S3		/x/S3.clean.fq.gz
```

- 列可留空（纯三代样本 r1/r2 留空）；每格支持逗号分隔多个文件（自动合并）
- 省略 tag 列时按文件后缀自动推断（fastq/isoseq/mix/bam/bam_isoseq/bam_mix）
- 模式A 自动生成的清单就是本格式，改完直接喂回来即可

**模式C：`-r` 原生描述文件透传**

高级用户直接提供 EviAnn 原生 `-r` 描述文件（一行一个实验：文件路径们 + 类型标签），模块不做任何加工。

### 其他证据（可选）

- **转录本 FASTA** `-e/--transcripts`：近缘物种的转录本序列，作为直接证据（与 `-r` 至少给一个）
- **蛋白质 FASTA** `-p/--proteins`：近缘物种蛋白（推荐 10+ 物种、总量约为预期基因数 5-10 倍）
- **UniProt** `-s/--uniprot`：功能注释用 SwissProt 库

## 参数说明 | Parameters { #parameters }

### 证据输入 | Evidence inputs

**通俗理解|In plain words:** 告诉 EviAnn「用哪些证据来猜基因」。给得越全、质量越高，注释越准；`--rnaseq-data` / `--sample-sheet` / `-r` 三选一（互斥），都不给时必须给 `-e` 转录本，否则程序拒绝运行。蛋白（`-p`）和 UniProt（`-s`）是辅助，可给可不给。

相关参数：`-g/--genome`（必需）、`-o/--output-dir`（必需）、`--rnaseq-data`、`--sample-sheet`、`-r/--rnaseq`、`-e/--transcripts`、`-p/--proteins`、`-s/--uniprot`、`-c/--cds-gff`、`--extra-gff`、`--mito-contigs`。

### 运行与注释选项 | Run & annotation options

**通俗理解|In plain words:** `-t` 线程数越大越快（默认 12 一般够用）；`-d` 倍性默认 2 绝大多数真核不用动；`-m` 是允许的最大内含子长度，不填就交给 EviAnn 自动判断，只在物种内含子特别长/特别短时才需要手调；`--lncrna-tpm` 越低保留的 lncRNA 越多、越高越严格，默认 1.0；`--min-prot` 是无同源证据时开放式找 ORF 的最小蛋白长度，默认 75 一般不用动。

相关参数：`-t/--threads`（默认 12）、`-d/--ploidy`（默认 2）、`-m/--max-intron`（默认自动）、`--lncrna-tpm`（默认 1.0）、`--min-prot`（默认 75）。

### 开关选项 | Flags

**通俗理解|In plain words:** `--partial` 打开后会把不完整（缺头少尾）的 CDS 也算进去，注释更全但可能掺入假阳性；`-f/--functional` 打开后额外跑功能注释（比对 SwissProt），更慢；`--debug`/`--verbose` 是排错用的，出问题时再开。

相关参数：`--partial`、`-f/--functional`、`--debug`、`--verbose`。

## 分析流程 | Pipeline { #pipeline }

```text
基因组 FASTA + 转录组数据(--rnaseq-data / --sample-sheet / -r)
    │
    ▼
自动识别与整理(类型识别/双端配对/多样本合并/mix匹配)
    │
    ▼
生成 EviAnn -r 描述文件 + 样本清单(01_inputs/, 可人工修改重跑)
    │
    ▼
conda run 调用 eviann.sh(work/ 目录内运行, 中间文件不外泄)
    │
    ▼
最终结果复制到 results/(work/ 保留供断点续传)
```

## 输出 | Output { #output }

```text
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml      # 模块与 EviAnn 版本
├── 01_inputs/
│   ├── sample_sheet.tsv           # 自动生成的样本清单(可修改后重跑)
│   └── rnaseq_list.txt            # EviAnn -r 描述文件
├── work/                          # EviAnn 运行目录(中间文件, 断点续传状态)
│   ├── genome.fa -> /原路径       # 基因组软链
│   ├── genome.fa.*                # EviAnn 各阶段中间产物
│   └── merged/                    # 多样本合并的 fastq
├── results/                       # 最终结果(复制自 work/)
│   ├── genome.fa.pseudo_label.gff # 基因结构注释(主结果)
│   ├── genome.fa.proteins.fasta   # 预测蛋白序列
│   └── genome.fa.transcripts.fasta# 预测转录本序列
└── 99_logs/
    └── eviann.log                 # 运行日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

- `*.pseudo_label.gff`：核心产物，含每个预测基因的坐标、外显子/内含子结构。`Evidence` 属性标注证据类型（`complete`=转录本+蛋白双证据最可靠 / `protein_only` / `transcript_only`）；`pseudo=true` 为假基因（无 CDS 输出）；`type=lncRNA` 为长非编码 RNA
- `*.proteins.fasta`：预测的蛋白序列，条数即蛋白编码转录本数，可用 BUSCO 等评估完整性
- `*.transcripts.fasta`：预测的转录本序列
- 日志 `99_logs/eviann.log` 记录了完整命令行、每类文件识别数量、mix 配对结果、每个实验行内容，排查问题先看它

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 文件名规范（`*.clean.fq.gz` 三代 / `*_1.clean.fq.gz`+`*_2.clean.fq.gz` 二代）：直接 `--rnaseq-data 目录/` 全自动
- 文件名不规范或 mix 配对结果不对：看 `01_inputs/sample_sheet.tsv`，改对后用 `--sample-sheet` 重跑
- 只有二代 RNA-seq：`--rnaseq-data` 指向目录即可，自动配对成双端行；只有三代全长：同样全自动识别
- 有近缘物种注释：用 `gffread` 提取转录本+蛋白，`-e` + `-p` 直接喂（lift-over 模式），无需 RNA-seq
- 默认参数即可覆盖多数物种；内含子极长的物种（如哺乳动物）可调大 `-m`
- 需要后续做功能富集：加 `-f`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 | Path | 基因组FASTA文件｜Genome FASTA file (required) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory (required) |
| `--rnaseq-data` | — |  | 转录组数据文件或目录(逗号分隔多个),自动识别二代/三代｜RNA-seq file(s) or dir(s), comma-separated |
| `--sample-sheet` | — | Path | 样本清单TSV｜Sample sheet TSV |
| `-r, --rnaseq` | — | Path | EviAnn原生-r描述文件(透传)｜EviAnn native -r file |
| `-e, --transcripts` | — | Path | 近缘物种转录本FASTA｜Transcripts FASTA |
| `-p, --proteins` | — | Path | 近缘物种蛋白质FASTA｜Proteins FASTA |
| `-s, --uniprot` | — | Path | UniProt-SwissProt FASTA｜UniProt-SwissProt FASTA |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-m, --max-intron` | — | int | 最大内含子长度｜Maximum intron length (default: auto) |
| `-d, --ploidy` | `2` | int | 基因组倍性｜Genome ploidy |
| `-c, --cds-gff` | — | Path | 含现有CDS的GFF｜GFF with existing CDS |
| `--lncrna-tpm` | `1.0` | float | lncRNA最小TPM｜Minimum TPM for lncRNA |
| `--min-prot` | — | int | 无同源证据时ab initio ORF最小蛋白长度(aa)｜Min protein length for ab initio ORF |
| `--partial` | `False` |  | 包含部分CDS｜Include partial CDS |
| `-f, --functional` | `False` |  | 执行功能注释｜Perform functional annotation |
| `--mito-contigs` | — | Path | 线粒体contig列表文件｜File with mitochondrial contigs |
| `--extra-gff` | — | Path | 额外GFF特征｜Extra features from external GFF |
| `--debug` | `False` |  | 保留中间文件｜Keep intermediate files |
| `--verbose` | `False` |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `--rnaseq-data` | — |  | 转录组数据文件或目录(逗号分隔多个),自动识别二代/三代｜RNA-seq file(s) or dir(s), comma-separated; auto-classified |
| `--sample-sheet` | — |  | 样本清单TSV(sample_id/r1/r2/long_reads,逗号分隔多文件)｜Sample sheet TSV |
| `-r, --rnaseq` | — |  | EviAnn原生-r描述文件(透传)｜EviAnn native -r description file |
| `-e, --transcripts` | — |  | 近缘物种转录本FASTA｜Transcripts FASTA from related species |
| `-p, --proteins` | — |  | 近缘物种蛋白质FASTA｜Proteins FASTA from related species |
| `-s, --uniprot` | — |  | UniProt-SwissProt FASTA｜UniProt-SwissProt FASTA |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `-m, --max-intron` | — | int | 最大内含子长度｜Maximum intron length (default: auto) |
| `-d, --ploidy` | `2` | int | 基因组倍性｜Genome ploidy (default: 2) |
| `-c, --cds-gff` | — |  | 含现有CDS的GFF｜GFF with existing CDS |
| `--lncrna-tpm` | `1.0` | float | lncRNA最小TPM｜Minimum TPM for lncRNA (default: 1.0) |
| `--min-prot` | — | int | 无同源证据时ab initio ORF最小蛋白长度(aa)｜Min protein length for ab initio ORF (default: 75) |
| `--partial` | — | store_true | 包含部分CDS｜Include partial CDS |
| `-f, --functional` | — | store_true | 执行功能注释｜Perform functional annotation |
| `--mito-contigs` | — |  | 线粒体contig列表文件｜File with mitochondrial contigs |
| `--extra-gff` | — |  | 额外GFF特征｜Extra features from external GFF |
| `--debug` | — | store_true | 保留中间文件｜Keep intermediate files |
| `--verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- EviAnn 软件：conda 环境 `eviann_v.2.0.5`（实际版本 2.0.4，纯证据驱动、无 de novo 基因预测器），通过 `bin/eviann.sh` 调用，环境内自带 minimap2/HISAT2/StringTie/gffread/BLAST/TransDecoder 等全套依赖
- Python 3（仅用于包装，无额外 Python 依赖）

## 常见问题 | FAQ { #faq }

**Q1：三种输入模式能给多个吗？**
不能，`--rnaseq-data` / `--sample-sheet` / `-r` 互斥，程序会报错。全自动结果不满意时，改 `01_inputs/sample_sheet.tsv` 换 `--sample-sheet` 模式重跑。

**Q2：我的三代文件叫 `S1.clean.fq.gz`，没带 isoseq 关键词，能识别吗？**
能。无关键词时按本实验室命名惯例识别：`*.clean.fq.gz`（无 `_1`/`_2` 配对痕迹）→ 三代；`*_1.clean.fq.gz`/`*_2.clean.fq.gz` → 二代双端。

**Q3：同一样本多个三代文件会怎样？**
样本名一致时（如 `S3_iso1.clean.fq.gz`/`S3_iso2.clean.fq.gz`，剥掉关键词后都叫 S3）自动二进制拼接成一个文件（gzip 流可直接拼接，安全），合并文件在 `work/merged/`。注意：如果批次后缀会残留在样本名里（如 `S3_iso_a` → 样本名 `S3_a`、`S3_iso_b` → `S3_b`），自动聚类会拆成两个实验——这不会报错但会拆散证据，此时在 `sample_sheet.tsv` 里把两个文件写进同一行逗号分隔，用 `--sample-sheet` 重跑。

**Q4：会断点续传吗？**
两层：模块层面 `results/` 三件套齐全时重跑直接跳过；EviAnn 层面 `work/` 中间文件保留，中断后重跑同一条命令会从最后完成的阶段继续。

**Q5：输出文件没生成怎么办？**
先看 `99_logs/eviann.log`：常见原因是基因组过大超内存、蛋白文件太小（建议 10+ 物种）、或 `-s` UniProt 库缺失（`-f` 时会自动联网下载，超算断网时需手动提供）。
