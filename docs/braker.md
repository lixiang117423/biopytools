# braker - BRAKER3 基因组注释 | BRAKER3 genome annotation

一句话理解：**把「基因组里哪些位置是基因」这件事交给 BRAKER3 自动学出来——给它基因组和至少一种证据(近缘蛋白 / 二代 RNA-seq / 三代转录本)，它就输出一套基因结构注释(GTF/GFF3 + 蛋白序列)**。
它解决的是「新测了一个基因组，怎么把基因位置标出来」的问题，全程自动完成重复序列屏蔽、转录组比对、基因结构预测三大步。

## 功能概述 | Overview { #overview }

- 基于 BRAKER3(整合 GeneMark-ETP 与 AUGUSTUS)的真核生物基因结构预测流程
- 自动完成重复序列屏蔽：RepeatModeler 建库 + RepeatMasker 软屏蔽
- 支持三类证据：近缘蛋白、三代全长转录本(Iso-Seq)、二代 RNA-seq，可单独或组合使用
- 内置 repeat 库过滤(剔除被误判成重复序列的基因家族)与证据还原(方案 1/2)
- 用 Singularity 镜像运行 BRAKER3，自动处理中文路径等非 ASCII 路径问题
- 断点续传：四个步骤各自按产物存在性跳过

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools braker -g genome.fa -s my_species -p proteins.fa
```

最小输入：基因组 FASTA + 物种名 + 至少一种证据(这里用近缘蛋白)。物种名用于输出文件命名；没有蛋白时可用 `--rnaseq_dirs`(二代 RNA-seq 目录)或 `--isoseq`(三代转录本)代替。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 基因结构注释 | 在基因组上标出基因的起止、外显子边界，并给出能翻译出的蛋白序列 |
| BRAKER3 | 主流基因预测软件，用「蛋白相似 + 转录本比对」自动训练模型再预测基因 |
| 证据 | 帮软件判断「这里是不是基因」的参考信息：近缘蛋白、转录本、RNA-seq |
| 重复序列屏蔽 | 把转座子等高重复区域盖住，避免干扰基因预测 |
| 软屏蔽 | 不是删掉重复区，而是把重复碱基写成小写字母，让软件知道「这里别太当真」 |
| 真菌模式 | 针对真菌/卵菌这类基因密集、内含子短的物种的特殊设置 |
| GTF / GFF3 | 两种基因注释标准格式，GTF 更常用，GFF3 兼容性更好 |
| Singularity | 一种容器技术，把软件和依赖打包成 .sif 文件，免去手动装依赖 |
| 断点续传 | 中断后重跑，已完成的步骤自动跳过 |

## 输入 | Input { #input }

### 基因组(-g，必填)

FASTA 格式的基因组序列。程序内部会先做重复屏蔽再用于预测。

### 证据(至少一种)

- **近缘蛋白(`--prot_seq`)**：一个 FASTA 文件或目录。给目录自动识别合并，给文件自动清理非标准氨基酸字符(如 `.`)
- **二代 RNA-seq(`--rnaseq_dirs`)**：逗号分隔的目录列表，目录里是成对 FASTQ，默认按 `_1.clean.fq.gz` / `_2.clean.fq.gz` 识别(可用 `--read1_pattern` / `--read2_pattern` 改)
- **三代转录本(`--isoseq`)**：三代全长转录本文件或目录

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 基因组 + 物种名是必填。物种名只是给输出文件起名用，取什么都行，但建议用拉丁名或可读的短名。

### 输入证据 | Input evidence

**通俗理解|In plain words:** 决定「靠什么来预测基因」，至少给一种。蛋白证据最快最省事，但会漏掉物种特有的新基因；RNA-seq 证据更全但要有测序数据。三者都给的组合效果最好。`--read1_pattern` / `--read2_pattern` 是二代数据双端文件的命名规律，**默认值适配 fastp 清洗后的标准命名，一般不用动**。

### 输出与运行参数 | Output & runtime

**通俗理解|In plain words:** `-o` 指定输出目录(默认 `./braker_output`)；`-t` 线程数只影响速度。`--fungus` 是真菌/卵菌模式开关——疫霉、霜霉这类卵菌物种**建议加 `--fungus`**，它们基因密集、内含子短，普通模式会预测不准。

### Singularity 容器 | Singularity

**通俗理解|In plain words:** BRAKER3 依赖复杂(含 GeneMark-ETP、AUGUSTUS、TSEBRA)，默认用 Singularity 镜像运行。镜像路径**一般不用动**；只有在自己已经装好 BRAKER3 时，才用 `--no_singularity` 走本机安装版。

### 步骤控制 | Step control

**通俗理解|In plain words:** 三个 skip 开关分别跳过「重复屏蔽」「三代转录本处理」「二代 RNA-seq 处理」。用于断点续传时跳过某一步，或想用外部准备好的屏蔽基因组/BAM 时。**正常完整运行不用动**。

### BRAKER3 高级参数 | BRAKER3 advanced

**通俗理解|In plain words:** 透传给 BRAKER3 的可选参数。`--busco_lineage` 指定 BUSCO 谱系(用于评估预测质量)、`--utr` 预测 UTR、`--training_genes` 提供训练基因集、`--use_existing` 复用已有参数。**新手一般不用动**。

### repeat 库过滤与证据还原 | repeat library filter & rescue

**通俗理解|In plain words:** 重复序列库里常混入被误判的真基因家族，`--skip_repeat_filter` 跳过这一步过滤(默认过滤)。`--skip_rescue` 控制「证据还原」——把被屏蔽步骤误盖掉的真基因区找回来，默认关闭(还原)。`--pfam_db` 是过滤时用的 Pfam 结构域库，**默认路径一般不用动**。另有 `--te_domain_evalue`、`--rescue_min_*`、`--prot_homology_*` 等细调参数只在直接调用 `python -m biopytools.braker` 时暴露，click 包装器用默认值。

## 分析流程 | Pipeline { #pipeline }

```text
输入: 基因组 + 证据(蛋白/二代RNA-seq/三代转录本)
    |
    v
步骤1: 重复序列屏蔽
  - RepeatModeler 构建重复库 -> RepeatMasker 软屏蔽 -> 屏蔽后基因组
  - (可选) repeat 库过滤: 剔除误判的基因家族
    |
    v
步骤2: 三代转录本处理(给了 isoseq 才跑)
  - minimap2 比对 -> isoseq.sorted.bam
    |
    v
步骤3: 二代 RNA-seq 处理(给了 rnaseq_dirs 才跑)
  - HISAT2 建索引 + 比对 -> rnaseq.sorted.bam
    |
    v
步骤3.5: 证据还原(可选, --no-skip_rescue 开启)
  - 用蛋白/转录证据找回被误屏蔽的真基因区
    |
    v
步骤4: BRAKER3 注释(Singularity 容器内)
  - 训练基因模型 -> 预测 -> TSEBRA 合并
  - 输出 braker.gtf / braker.gff3 / braker.aa / braker.codingseq
```

## 输出 | Output { #output }

```text
output_dir/
├── 01_repeat_masking/               # 步骤1: 重复屏蔽产物
│   ├── genome.fa.masked             # 屏蔽后的基因组(软屏蔽, 小写=重复区)
│   └── genome.fa.out                # RepeatMasker 报告
├── 02_long_reads/                   # 步骤2: 三代转录本比对(给了 isoseq 才有)
│   └── isoseq.sorted.bam
├── 03_short_reads/                  # 步骤3: 二代 RNA-seq 比对(给了 rnaseq_dirs 才有)
│   └── rnaseq.sorted.bam            # 排序后的 BAM + 索引
├── 04_braker_annotation/            # 步骤4: 最终注释结果(核心产物)
│   ├── braker.gtf                   # 基因结构 GTF 格式
│   ├── braker.gff3                  # 基因结构 GFF3 格式
│   ├── braker.aa                    # 预测蛋白序列
│   ├── braker.codingseq             # 预测 CDS 序列
│   ├── braker.log                   # BRAKER3 运行日志
│   ├── hintsfile.gff                # 证据提示文件
│   ├── Augustus/                    # AUGUSTUS 训练产物
│   └── GeneMark-ETP/                # GeneMark-ETP 训练产物
└── logs/                            # 流程日志
    ├── braker_main.log
    └── braker_pipeline.log
```

## 结果解读 | Interpreting Results { #results }

### 1. 基因注释(`04_braker_annotation/braker.gtf`)

**通俗理解|In plain words:** 这是最终答案——每条记录标出一个基因/转录本/外显子的坐标和结构。下游功能注释、比较基因组都从这里开始。

- `braker.gtf`：标准基因结构，交给下游工具(如 eggnog-mapper、interproscan)做功能注释
- `braker.aa`：每个预测基因的蛋白序列，做功能注释时直接用它
- `braker.codingseq`：每个预测基因的 CDS(编码区)核酸序列
- 基因数量是否合理：可用 `--busco_lineage` 让 BRAKER3 顺便评估完整性，或单独跑 `biopytools busco` 用 `braker.aa` 验证(完整度 >80% 通常算合格)

### 2. 日志与质量

**通俗理解|In plain words:** 日志记录每一步做了什么、用了什么命令。BRAKER 阶段最关键的是看它有没有正常收敛(训练+预测完成)，`braker.log` 末尾有预测统计。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规真核物种**：`-g` + `-s` + 蛋白/RNA-seq 任一证据即可，全默认
- **真菌/卵菌**：加 `--fungus`(疫霉、霜霉等卵菌必加)
- **只有近缘蛋白、没测序数据**：只给 `--prot_seq`，最快，但会漏物种特有基因
- **有 RNA-seq 且想注释更全**：`--prot_seq` + `--rnaseq_dirs` 一起给
- **重复区基因重要(如效应子)**：用 `biopytools braker4phyto` 或加 `--skip_repeat`，避免屏蔽掉真基因
- **换参数重跑**：先删对应步骤的产物，否则断点续传会复用旧结果(见 FAQ)

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--species, -s` | 必填 |  | 物种名称｜Species name (for BRAKER output naming) |
| `--prot_seq, -p` | — |  | 近缘物种蛋白质序列文件或文件夹｜Protein sequences file or directory |
| `--isoseq, -l` | — |  | 三代全长转录本文件夹｜Long-read transcript directory |
| `--rnaseq_dirs` | — |  | 二代RNA-seq目录列表(逗号分隔)｜Comma-separated RNA-seq directories |
| `--read1_pattern` | `_1.clean.fq.gz` |  | R1文件模式｜R1 file pattern |
| `--read2_pattern` | `_2.clean.fq.gz` |  | R2文件模式｜R2 file pattern |
| `--output_dir, -o` | `./braker_output` |  | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--fungus` | — |  | 使用真菌模式｜Use fungus mode (suitable for oomycetes) |
| `--singularity_image` | `~/software/singularity/braker3_devel.sif` |  | Singularity镜像路径｜Singularity image path |
| `--no_singularity` | — |  | 不使用Singularity镜像｜Do not use Singularity image |
| `--skip_repeat` | — |  | 跳过重复序列屏蔽｜Skip repeat masking |
| `--skip_long_reads` | — |  | 跳过三代转录本处理｜Skip long-read processing |
| `--skip_short_reads` | — |  | 跳过二代RNA-seq处理｜Skip short-read processing |
| `--busco_lineage` | — |  | BUSCO谱系｜BUSCO lineage |
| `--utr` | — |  | 预测UTR｜Predict UTR |
| `--training_genes` | — |  | 训练基因集｜Training gene set file |
| `--use_existing` | — |  | 使用已有参数｜Use existing parameters |
| `--skip_repeat_filter` | — |  | 跳过repeat库过滤(方案1,默认开启)｜Skip repeat library filtering |
| `--skip_rescue/--no-skip_rescue` | `True` |  | 跳过证据还原(默认关闭,--no-skip_rescue开启)｜Skip rescue (default off) |
| `--pfam_db` | `~/database/eggnog/pfam/Pfam-A.hmm` |  | Pfam-A HMM 库路径｜Pfam-A HMM DB path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--species, -s` | 必填 |  | 物种名称｜Species name (用于BRAKER输出命名)｜used for BRAKER output naming |
| `--prot_seq, -p` | — |  | 近缘物种蛋白质序列文件或文件夹｜Protein sequences file or directory |
| `--isoseq, -l` | — |  | 三代全长转录本文件夹｜Long-read transcript directory |
| `--rnaseq_dirs` | — |  | 二代RNA-seq数据目录列表，逗号分隔｜Comma-separated list of RNA-seq directories |
| `--read1_pattern` | `_1.clean.fq.gz` |  | R1文件模式｜R1 file pattern (default: _1.clean.fq.gz) |
| `--read2_pattern` | `_2.clean.fq.gz` |  | R2文件模式｜R2 file pattern (default: _2.clean.fq.gz) |
| `--output_dir, -o` | `./braker_output` |  | 输出目录｜Output directory (default: ./braker_output) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--fungus` | — | store_true | 使用真菌模式｜Use fungus mode (suitable for oomycetes) |
| `--singularity_image` | `~/software/singularity/braker3_devel.sif` |  | Singularity镜像路径｜Singularity image path |
| `--no_singularity` | — | store_true | 不使用Singularity镜像｜Do not use Singularity image |
| `--skip_repeat` | — | store_true | 跳过重复序列屏蔽｜Skip repeat masking step |
| `--skip_long_reads` | — | store_true | 跳过三代转录本处理｜Skip long-read processing step |
| `--skip_short_reads` | — | store_true | 跳过二代RNA-seq处理｜Skip short-read processing step |
| `--busco_lineage` | — |  | BUSCO谱系｜BUSCO lineage |
| `--utr` | — | store_true | 预测UTR｜Predict UTR regions |
| `--training_genes` | — |  | 训练基因集文件｜Training gene set file |
| `--use_existing` | — | store_true | 使用已有参数｜Use existing parameters |
| `--skip_repeat_filter` | — | store_true | 跳过repeat库过滤(方案1,默认开启)｜Skip repeat library filtering |
| `--skip_rescue` | `True` |  | 跳过证据还原(默认关闭,filter库级已处理;--no-skip_rescue开启)｜Skip rescue (default off) |
| `--pfam_db` | `~/database/eggnog/pfam/Pfam-A.hmm` |  | Pfam-A HMM 库路径｜Pfam-A HMM DB path |
| `--te_domain_evalue` | `1e-05` | float | TE domain hmmscan E-value 阈值｜TE domain E-value cutoff |
| `--filter_min_orf_len` | `30` | int | 过滤用最小ORF长度(aa)｜Min ORF length (aa) for filter |
| `--rescue_min_cds_len` | `100` | int | rescue蛋白证据最小覆盖长度(bp)｜Min CDS overlap (bp) for rescue |
| `--rescue_min_identity` | `70` | float | rescue蛋白最小identity(%%)｜Min protein identity (%%) for rescue |
| `--rescue_min_depth` | `5` | int | rescue RNA-seq最小覆盖度｜Min RNA-seq depth for rescue |
| `--prot_homology_evalue` | `1e-05` | float | prot_seq 同源 E-value 阈值｜Protein homology E-value |
| `--prot_homology_pident` | `50.0` | float | prot_seq 同源 identity 阈值(%%)｜Protein homology identity (%%) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- **BRAKER3**(含 GeneMark-ETP、AUGUSTUS、TSEBRA)：Singularity 镜像 `~/software/singularity/braker3_devel.sif`
- **Singularity**：`~/miniforge3/envs/singularity_v.3.8.7/bin/singularity`
- **RepeatModeler / RepeatMasker / BuildDatabase**：conda 环境 `repeat`
- **minimap2**：conda 环境 `align`
- **HISAT2 / hisat2-build**：conda 环境 `rna`
- **hmmscan**：conda 环境 `braker_v.3.0.8`(repeat 库过滤用)
- **mmseqs / miniprot**：conda 环境 `annot`
- **samtools**：系统自动检测
- **Pfam-A.hmm 数据库**：默认 `~/database/eggnog/pfam/Pfam-A.hmm`(repeat 库过滤用)

## 常见问题 | FAQ { #faq }

**Q1：中断后重跑会从头再来吗？**
不会。重复屏蔽、三代比对、二代比对、BRAKER 注释四步各自按产物存在性判断，已完成的自动跳过。

**Q2：换参数重跑，为什么结果没变？**
断点续传按产物存在性判断。换参数(如加/去 `--fungus`、换证据)重跑前，先删除对应步骤的产物(如 `04_braker_annotation/braker.gtf`)，否则会复用旧结果。

**Q3：报「至少提供一种证据」是什么问题？**
`-g` + `-s` 之外必须再给一种证据(蛋白/三代/二代任选其一)，否则程序不知道靠什么预测基因。

**Q4：`--te_domain_evalue`、`--rescue_min_cds_len` 这些参数命令行里为什么没有？**
这些 repeat 库过滤/证据还原的细调参数只在直接调用 `python -m biopytools.braker` 时暴露，click 包装器用默认值。需要细调时走模块直调。

**Q5：为什么默认用 Singularity？可以不用吗？**
BRAKER3 依赖 GeneMark-ETP、AUGUSTUS、TSEBRA 等一堆软件，版本要求严格，容器化最省事。已装好本机版时用 `--no_singularity` 走本机。

**Q6：运行很久、BRAKER 阶段失败怎么办？**
看 `logs/braker_pipeline.log` 和 `braker.log`。常见原因：中文路径(程序已自动用 ~/tmp 软链接规避)、species 名残留(重跑需清 AUGUSTUS species 目录，程序已自动处理)、镜像缺失(检查 `--singularity_image` 路径)。
