
# ANNOVAR - 变异功能注释 | Variant Functional Annotation

一句话理解：**把 VCF 里的每个变异（SNP/INDEL）翻译成「它落在基因的哪个部位、会不会改变蛋白」**，解决「测序找到了一堆变异，却不知道哪个真正影响基因功能」的问题。

## 功能概述 | Overview { #overview }

- 输入 VCF + GFF3 注释 + 参考基因组，输出每个变异的功能注释（外显子/内含子/剪接位点/基因间区等）
- 四步流水线：GFF3 转 GenePred → 提取转录本/蛋白/CDS 序列 → VCF 处理与格式转换 → ANNOVAR 注释
- 内置 GFF3 自动清理与修复（坐标、CDS phase），全程在工作副本上操作，不修改用户的输入文件
- 支持 `--step` 单独重跑某一步；断点续传按输出文件存在性自动跳过已完成步骤
- 注释结果自动后处理成 TSV 表格，方便直接查看与下游分析

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools annovar -i variants.vcf -g annotation.gff3 -G genome.fa -b OV
```

最小输入：一个 VCF（变异）+ 一个 GFF3（基因注释）+ 一个参考基因组 FASTA + 一个基因组版本标识（如 OV）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| VCF | 变异清单文件，每行记一个「基因组哪里、从什么变成了什么」 |
| GFF3 | 基因注释文件，像基因组的「户型图」，标出每个基因/外显子在哪 |
| 外显子(exon) | 真正编码蛋白的片段，变异落在这里影响最大 |
| 内含子(intron) | 外显子之间的「间隔区」，转录后会被剪掉 |
| 剪接位点 | 外显子与内含子交界处的「铆点」，变异可能破坏正常剪接 |
| GenePred/refGene | ANNOVAR 认识的一种基因注释格式，由 GFF3 转换而来 |
| buildver | 基因组版本标识，给这组注释文件命名用（如 OV、KY131） |
| ANNOVAR | 本工具的核心注释引擎，把变异和基因注释对上号 |

## 输入 | Input { #input }

三个必需文件 + 一个版本标识：

- **VCF**（`-i`）：标准 VCF，支持 `.vcf` / `.vcf.gz`
- **GFF3**（`-g`）：基因注释文件，需含 gene 等特征
- **参考基因组**（`-G`）：FASTA 格式，序列名须与 GFF3、VCF 的染色体名一致
- **buildver**（`-b`）：基因组版本标识，仅作为输出文件名前缀，不能含路径分隔符（如 `OV`、`KY131`）

格式示例（VCF）：

```text
#CHROM  POS     ID  REF  ALT  QUAL  FILTER  INFO
Chr01   732093  .   G    A    99    PASS    .
```

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 这是四条缺一不可的输入：变异在哪（VCF）、基因怎么画（GFF3）、基因组的碱基序列（FASTA）、以及给这批注释起个版本名（buildver）。buildver 只是文件名前缀，取个有意义的名字即可。

### 软件与输出路径 | Paths

**通俗理解|In plain words:** `--annovar-path` 是 ANNOVAR 软件的安装目录（里面是一堆 perl 脚本），一般用默认值；`--output-dir` 既放最终结果、又充当 ANNOVAR 的「数据库目录」（refGene 基因注释文件会生成在这里），所以别把它和别的项目混用。

### VCF 处理 | VCF processing

**通俗理解|In plain words:** `--qual-threshold` 决定「质量分低于多少的变异先丢掉」，只在启用 VCF 过滤时才生效。**默认不启用 VCF 过滤**，因为多数输入的 VCF 已经过滤过了；只有输入是「原始未过滤」的 VCF 时才用 `--enable-vcf-filter` 打开。

### GFF3 处理 | GFF3 processing

**通俗理解|In plain words:** 程序会先对 GFF3 做一次「清理 + 修复」（修坐标、补 CDS phase），再拿去转换，且全程在工作副本上做、绝不改你的原始文件。只有当你确信自己的 GFF3 完全规范、或修复反而出错时，才用 `--skip-gff-cleaning` / `--skip-gff-fix` 跳过。**一般不用动。**

### 步骤控制 | Step control

**通俗理解|In plain words:** `--step 1~4` 让你只重跑某一步（1=转 GenePred，2=提序列，3=VCF 处理，4=注释），适合中间某步失败后单独补跑，不用从头再来。

## 分析流程 | Pipeline { #pipeline }

```text
输入 VCF + GFF3 + 参考基因组 + buildver
    │
    ▼
步骤1: GFF3 清理修复 → 转 GenePred ({buildver}_refGene.txt)
    │
    ▼
步骤2: 提取转录本/蛋白/CDS 序列 (retrieve_seq_from_fasta.pl + gffread)
    │
    ▼
步骤3: VCF 处理(可选过滤 QUAL) → 转 ANNOVAR 格式 ({vcf}.annovar.vcf)
    │
    ▼
步骤4: annotate_variation.pl --geneanno 基因注释
    │
    ▼
结果后处理 → 注释 TSV 表格 + annotation_summary.txt
```

## 输出 | Output { #output }

```text
output_dir/
├── {buildver}_refGene.txt            # 转换后的基因注释(ANNOVAR 数据库)
├── {buildver}.cleaned.gff3           # GFF3 清理后的工作副本
├── {buildver}_refGeneMrna.fa         # 转录本序列
├── {buildver}_refGenePep.fa          # 蛋白序列
├── {buildver}_refGeneCds.fa          # CDS 序列
├── {vcf}.annovar.vcf                 # 转成 ANNOVAR 格式的变异
├── {vcf}.exonic_variant_function     # 外显子区变异注释(核心结果)
├── {vcf}.variant_function            # 全变异功能注释
├── {vcf}.log                         # ANNOVAR 运行日志
├── {vcf}_processed_exonic.tsv        # 后处理的外显子注释表
├── {vcf}_processed_all.tsv           # 后处理的全量注释表
└── annotation_summary.txt            # 注释总结报告
```

其中 {vcf} 是输入 VCF 去掉 .vcf[.gz] 后的名字，{buildver} 是 -b 指定的版本标识。

## 结果解读 | Interpreting Results { #interpreting }

### 1. {vcf}.exonic_variant_function（外显子注释）

**通俗理解|In plain words:** 这是最核心的结果——落在「蛋白编码区」的变异清单，每个变异标注它会不会改变蛋白（错义/无义/移码等）。做「致病/重要变异」筛选时先看这里。

### 2. {vcf}.variant_function（全量注释）

**通俗理解|In plain words:** 所有变异的「住址」分类表：外显子、内含子、剪接位点、UTR、基因间区等，一网打尽。适合看整体分布。

### 3. {vcf}_processed_*.tsv（后处理表）

**通俗理解|In plain words:** 把 ANNOVAR 的原始输出整理成规整的表格，方便用 Excel/R 打开、筛选、排序。

### 4. annotation_summary.txt（总结报告）

**通俗理解|In plain words:** 记录了本次运行的输入、参数、输出了哪些文件，相当于一张回执单。

## 参数选择建议 | Parameter Guidance { #guidance }

- **-b buildver**：取个能代表「这个物种/这个版本注释」的名字即可，如物种缩写；不要带路径分隔符。
- **--enable-vcf-filter**：只有当输入 VCF 是「未过滤原始结果」时才打开，配合 --qual-threshold 设置质量线（默认 20）。
- **--skip-gff-cleaning / --skip-gff-fix**：GFF3 转换失败报错、或你确定注释很规范时才考虑，平时保持默认（不跳过）。
- **--step**：某步失败时用 --step N 单独补跑；注意断点续传是「看输出文件在不在」，换参数重跑要先删对应旧产物。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | VCF变异文件｜VCF variant file path |
| `--gff3, -g` | 必填 |  | GFF3注释文件｜GFF3 annotation file path |
| `--genome, -G` | 必填 |  | 基因组序列文件｜Genome sequence file path |
| `--build-ver, -b` | 必填 |  | 基因组版本(例如: OV, KY131)｜Genome build version identifier (e.g., OV, KY131) |
| `--annovar-path, -a` | `~/software/annovar/annovar` |  | ANNOVAR软件路径｜ANNOVAR software installation path |
| `--output-dir, -o` | `./annovar_output` | Path | 输出目录(同时作为ANNOVAR数据库目录)｜Output directory (also the ANNOVAR database dir) |
| `--qual-threshold, -q` | `20` | int | VCF质量阈值｜VCF quality filtering threshold |
| `--step, -s` | — | IntRange | 运行指定步骤｜Run only specified step: 1: GFF3转换｜GFF3 conversion 2: 提取序列｜Extract sequences 3: VCF处理｜VCF processing 4: 变异注释｜Variant annotation |
| `--skip-gff-cleaning` | — |  | 跳过GFF3清理｜Skip GFF3 file format cleaning |
| `--skip-gff-fix` | — |  | 跳过GFF3修复｜Skip automatic GFF3 file fixes |
| `--enable-vcf-filter` | — |  | 启用VCF过滤(默认跳过)｜Enable VCF filtering step (skipped by default) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | VCF变异文件路径｜VCF variant file path |
| `-g, --gff3` | 必填 |  | GFF3注释文件路径｜GFF3 annotation file path |
| `-G, --genome` | 必填 |  | 基因组序列文件路径｜Genome sequence file path |
| `-b, --build-ver` | 必填 |  | 基因组构建版本标识符(如: OV, KY131) - 不应包含路径分隔符｜Genome build version identifier (e.g., OV, KY131) - should not contain path separators |
| `-a, --annovar-path` | `~/software/annovar/annovar` |  | ANNOVAR软件安装路径｜ANNOVAR software installation path |
| `-o, --output-dir` | `./annovar_output` |  | 输出目录(同时作为ANNOVAR数据库目录,refGene在此生成)｜Output directory (also the ANNOVAR database dir; refGene is built here) |
| `-q, --qual-threshold` | `20` | int | VCF质量过滤阈值(仅在启用VCF过滤时生效)｜VCF quality filtering threshold (only effective when VCF filtering is enabled) |
| `-s, --step` | — | 1/2/3/4 | 只运行指定步骤｜Run only specified step (1: gff3转换｜gff3 conversion, 2: 提取序列｜extract sequences, 3: VCF处理｜VCF processing, 4: 注释｜annotation) |
| `--skip-gff-cleaning` | — | store_true | 跳过GFF3文件的格式清理(attributes清理和坐标修复)｜Skip GFF3 file format cleaning (attributes cleaning and coordinate fixing) |
| `--skip-gff-fix` | — | store_true | 跳过GFF3文件的自动修复(CDS phase等问题)｜Skip automatic GFF3 file fixes (CDS phase and other issues) |
| `--skip-vcf-filter` | `True` | store_true | 跳过VCF过滤步骤，直接使用输入的VCF文件(默认启用)｜Skip VCF filtering step, use input VCF file directly (enabled by default) |
| `--enable-vcf-filter` | — | store_true | 启用VCF过滤步骤(使用bcftools)｜Enable VCF filtering step (using bcftools) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- ANNOVAR（默认 ~/software/annovar/annovar，内含 annotate_variation.pl、retrieve_seq_from_fasta.pl、convert2annovar.pl 等 perl 脚本；可用 --annovar-path 或 ANNOVAR_PATH 环境变量覆盖）
- gff3ToGenePred（UCSC kent 工具，经 conda 环境调用）
- gffread（提取蛋白/CDS 序列，经 conda 环境调用；可用 GFFREAD_PATH 覆盖）
- seqkit（结果后处理，可用 SEQKIT_PATH 覆盖）
- bcftools（仅启用 VCF 过滤时需要）
- Python 3

## 常见问题 | FAQ { #faq }

**Q1：换参数重跑，结果没变？**
断点续传按「输出文件是否存在」判断。换过滤参数（如 --qual-threshold）重跑前，先删除旧的 {vcf}.annovar.vcf 等对应产物，否则会复用旧结果。

**Q2：报错「GenePred 缺失某些染色体的基因模型」？**
GFF3 转 GenePred 不完整，通常是 GFF3 格式不规范（坐标、feature 类型问题）。先别急着 --skip-gff-fix，检查输入 GFF3；程序默认会自动清理修复。

**Q3：buildver 能写路径吗？**
不能。它只是文件名前缀，含路径分隔符会被自动剥成 basename。

**Q4：为什么默认跳过 VCF 过滤？**
因为多数人输入的 VCF 已经是过滤后的；--qual-threshold 只在 --enable-vcf-filter 打开时才生效。

**Q5：会改我的 GFF3 吗？**
不会。所有清理/修复都在 output_dir 内的 {buildver}.cleaned.gff3 工作副本上进行，输入文件保持原样。
