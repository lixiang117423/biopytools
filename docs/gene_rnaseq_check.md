# 候选基因 RNA-seq 转录验证 | Candidate Gene RNA-seq Transcriptional Validation

一句话理解：**给一批「候选基因」做 RNA 层面的体检**——用 RNA-seq 数据逐个检查这些基因到底有没有表达、结构对不对、边界准不准，最后给每个基因打一个结论标签。

## 功能概述 | Overview { #overview }

- 针对目标基因列表做定向转录验证：比对(HISAT2)、覆盖度(bedtools)、剪接 junction、StringTie 组装 + gffcompare
- 自动检测建库链特异性（infer_experiment.py），链特异时自动重比对
- 每个基因打 5 类标签：完整表达 / 部分表达 / 边界问题 / 未表达 / 无法判定
- 输出主报告 TSV + 需复核基因的 BED + 边界修改建议
- 分步骤执行（可 `--steps` 指定子集），断点续传（`--force` 强制重跑）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools gene-rnaseq-check -g genome.fa -a anno.gff -e genes.txt -r ./reads/ -o out
```

最小输入：基因组 FASTA + GFF3 注释 + 基因 ID 列表 + RNA-seq reads 目录（配对 FASTQ）。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| 候选基因 | 你重点怀疑或关注的基因列表，工具只检查这些，不全基因组铺开 |
| 链特异性建库 | RNA-seq 建库时记录每条 read 来自哪条链；比对时搞反链会出错 |
| 覆盖度(breadth/depth) | 基因被 reads 盖住多少（广度）、盖得有多厚（深度），衡量有没有表达 |
| junction read | 横跨两个外显子接缝的 read，是「这个剪接方式真实存在」的直接证据 |
| StringTie | 从比对结果重新组装转录本的软件，用来发现注释里没有的新转录本 |
| gffcompare class_code | 组装结果与参考注释的对应关系代码（`=` 完全一致，`c` 部分一致等） |
| 5 类结论标签 | 工具给每个基因的最终体检结论（见结果解读） |

## 输入 | Input { #input }

- **基因组** `-g`：FASTA
- **注释** `-a`：GFF3（含 gene / mRNA / exon / CDS 特征）
- **基因列表** `-e`：纯文本，每行一个基因 ID（`#` 开头忽略），ID 须与 GFF3 里 gene 的 `ID` 一致（兼容 `gene:xxx` 前缀格式）
- **reads 目录** `-r`：配对 FASTQ，自动识别常见命名（`_1/_2`、`R1/R2`、`.1/.2` 等）；命名特殊时用 `--reads-pattern`

## 参数说明 | Parameters { #parameters }

### 运行资源 | Resources

**通俗理解|In plain words:** `-t` 线程数默认 12；`--sample-timeout` 是单样本超时（秒，默认 6 小时），样本大或机器慢时调大，避免被误杀。**一般不用动。**

相关参数：`-t/--threads`（默认 12）、`--sample-timeout`（默认 21600）。

### 分析窗口与容差 | Window & tolerance

**通俗理解|In plain words:** `--flanking-window` 是基因上下游各取多长来判断「边界外有没有不该有的表达」（默认 500 bp），调大能发现更远的边界渗漏但更保守；`--junction-tolerance` 是判断剪接位点是否匹配注释的容差（bp），越大越宽松；`--strandness-confidence` 是链特异性判定的置信度门槛（%），低于它就用不区分链模式。

相关参数：`--flanking-window`（默认 500）、`--junction-tolerance`（默认 5）、`--strandness-confidence`（默认 70.0）。

### 步骤控制与日志 | Steps & logging

**通俗理解|In plain words:** `--steps` 默认 `all` 跑全流程，也可只跑子集（如 `parse_gff,coverage,classify`）快速出结论；`--force` 强制重跑所有已完成的步骤；`--reads-pattern` 用于非常规 FASTQ 命名；`-v`/`--quiet` 控制日志详略。

相关参数：`--steps`（默认 all）、`--force`、`--reads-pattern`、`-v/--verbose`、`--quiet`。

## 分析流程 | Pipeline { #pipeline }

```text
基因组 + GFF3 + 基因列表 + reads
    │
    ▼
解析 GFF3 目标基因 → HISAT2 建索引
    │
    ▼
HISAT2 比对 → 链特异性检测(必要时重比对)
    │
    ▼
bedtools 覆盖度 + junction reads 分析
    │
    ▼
StringTie 组装 + gffcompare 分类
    │
    ▼
综合分类 → 主报告 + 复核 BED + 边界建议
```

## 输出 | Output { #output }

```text
out/
├── 00_pipeline_info/                    # 软件版本与参数
├── 01_align/                            # 索引、各样本 {sample}_sorted.bam/.bai/.flagstat
├── 02_coverage/
│   ├── effector_*.bed                   # exon/gene/上下游 BED
│   ├── coverage_summary.tsv             # 每基因覆盖度汇总
│   └── coverage_per_exon.tsv            # 每外显子覆盖度
├── 03_junctions/
│   └── junction_analysis.tsv            # 剪接支持/新剪接
├── 04_stringtie/
│   ├── {sample}.gtf / merged.gtf        # 组装转录本
│   └── gffcompare/                      # 与参考注释比较
├── 05_results/
│   ├── gene_validation_report.tsv       # 主报告(核心)
│   ├── all_target_genes.bed             # 全部目标基因 BED(IGV 查看)
│   ├── boundary_suggestions.tsv         # 边界修改建议
│   └── flagged_for_review/              # 需复核基因,按类别分 BED
└── 99_logs/                             # 运行日志
```

主报告列：`gene_id`、`chrom`、`start/end`、`strand`、`gene_length`、`num_exons`、`overall_coverage_pct`、`mean_exon_coverage_depth`、`min_exon_coverage_pct`、`upstream/downstream_500bp_mean_depth`、`junction_support`、`novel_junctions`、`stringtie_class_code`、`classification`、`needs_review`、`notes`。

## 结果解读 | Interpreting Results { #interpreting-results }

- **EXPRESSED_COMPLETE（完整表达）**：覆盖度达标且 gffcompare 完全/部分一致，注释可信，无需处理
- **EXPRESSED_PARTIAL（部分表达）**：有些外显子覆盖低或整体覆盖不完整，可能是注释多了/少了外显子，建议复核
- **BOUNDARY_ISSUE（边界问题）**：基因边界外有表达或出现新剪接，暗示起止位置可能标错
- **NOT_EXPRESSED（未表达）**：基因本体与两侧都几乎没 reads，可能是假基因、该条件不表达，或注释错了
- **AMBIGUOUS（无法判定）**：证据不足以归类，需人工判断
- `needs_review=True` 的基因（部分表达/边界问题/无法判定）会自动整理进 `flagged_for_review/` 的 BED，可直接拖进 IGV 逐个看

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 快速粗筛：`--steps parse_gff,align,coverage,classify`（跳过 StringTie，最快拿到覆盖度结论）
- 完整验证：默认 `all`，含 StringTie 组装与边界建议
- 换参数重跑：加 `--force` 强制刷新（否则断点续传会复用旧结果）
- 命名非标准：用 `--reads-pattern`（如 `*_1.fq.gz`，`*` 前为样本名前缀）

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-a, --annotation` | 必填 |  | GFF3注释文件｜GFF3 annotation file |
| `-e, --gene-list` | 必填 |  | 目标基因ID列表文件｜Target gene ID list file |
| `-r, --reads-dir` | 必填 |  | RNA-seq reads目录｜RNA-seq reads directory |
| `-o, --output-dir` | `./gene_rnaseq_check_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` |  | 线程数｜Threads |
| `--reads-pattern` | `` |  | FASTQ命名模式｜FASTQ naming pattern |
| `--steps` | `all` |  | 执行步骤(逗号分隔或all)｜Steps (comma-separated or all) |
| `--sample-timeout` | `21600` |  | 单样本超时时间(秒)｜Sample timeout in seconds |
| `--flanking-window` | `500` |  | 上下游分析窗口(bp)｜Flanking analysis window (bp) |
| `--junction-tolerance` | `5` |  | Junction容差(bp)｜Junction tolerance (bp) |
| `--strandness-confidence` | `70.0` |  | 链特异性判定置信度(%%)｜Strandness confidence threshold (%%) |
| `-v, --verbose` | — |  | 详细模式｜Verbose mode |
| `--quiet` | — |  | 静默模式｜Quiet mode |
| `--force` | — |  | 强制重新运行｜Force re-run |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-a, --annotation` | 必填 |  | GFF3注释文件｜GFF3 annotation file |
| `-e, --gene-list` | 必填 |  | 目标基因ID列表文件｜Target gene ID list file |
| `-r, --reads-dir` | 必填 |  | RNA-seq reads目录｜RNA-seq reads directory |
| `-o, --output-dir` | `./gene_rnaseq_check_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--reads-pattern` | `` |  | FASTQ命名模式｜FASTQ naming pattern (e.g. *_1.fq.gz) |
| `--steps` | `all` |  | 执行步骤(逗号分隔或all)｜Steps (comma-separated or all) |
| `--sample-timeout` | `21600` | int | 单样本超时时间(秒)｜Sample timeout in seconds |
| `--flanking-window` | `500` | int | 上下游分析窗口(bp)｜Flanking analysis window (bp) |
| `--junction-tolerance` | `5` | int | Junction容差(bp)｜Junction tolerance (bp) |
| `--strandness-confidence` | `70.0` | float | 链特异性判定置信度(%%)｜Strandness confidence threshold (%%) |
| `-v, --verbose` | — | store_true | 详细模式｜Verbose mode |
| `--quiet` | — | store_true | 静默模式｜Quiet mode |
| `--force` | — | store_true | 强制重新运行｜Force re-run |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- HISAT2、hisat2-build、extract_splice_sites.py、extract_exons.py、StringTie、infer_experiment.py（conda 环境 `rna`，默认 `~/miniforge3/envs/rna/bin/`）
- samtools、bedtools（conda 环境 `align`）
- gffcompare（conda 环境 `annot`）
- Python 3 + pysam（junction 分析读 BAM）

## 常见问题 | FAQ { #faq }

**Q1：会断点续传吗？**
会。HISAT2 索引、BAM、StringTie GTF 等按输出文件存在性跳过；`--force` 强制全部重跑。换参数（如 `--flanking-window`）重跑前建议 `--force`，否则会复用旧覆盖度结果。

**Q2：报「未检测到 RNA-seq 样本」？**
reads 目录里的命名不在默认配对模式里，且没给 `--reads-pattern`。确认文件是成对 FASTQ，命名类似 `sample_1.fq.gz`/`sample_2.fq.gz`，或显式指定 `--reads-pattern`。

**Q3：基因列表里的 ID 匹配不上？**
列表 ID 必须等于 GFF3 中 gene 行的 `ID`（工具兼容 `gene:xxx` 这种带前缀的写法）。日志会打印「未在 GFF3 中找到」的缺失 ID，据此核对。

**Q4：为什么有的基因结论是「无法判定」？**
覆盖度既不够「完整」也不够「未表达」、又没有明确边界/剪接问题时，归为 AMBIGUOUS。属于正常兜底分类，人工在 IGV 里看覆盖图判断即可。
