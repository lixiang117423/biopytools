# gffcompare - GFF/GTF 两两比较分析 | GFF/GTF Pairwise Comparison

一句话理解：**把多个（2 到 5 个）GFF/GTF 注释文件两两互相比较，用 gffcompare 算出每个注释相对另一个被归入哪一类（完全一致、部分重叠、新转录本……），再汇总成一张表**，用来评估不同注释之间的差异和一致性。

## 功能概述 | Overview

- 对 2 到 5 个 GFF/GTF 文件做两两「双向」比较（A 比 B、B 比 A 各一次，共 n×(n-1) 对）
- 底层调用 gffcompare，输出每个比较对的分类统计（stats）、匹配追踪（tracking）、合并注释（annotated.gtf）等
- 自动识别文件或目录输入，目录下自动扫描 .gff/.gtf/.gff3 及 .gz 压缩版本
- 把全部比较对的 .stats 合并成一张 all_stats.tsv，方便横向对比
- 支持断点续传：已完成的比较对自动跳过，--force 强制重跑
- 生成 software_versions.yml 记录软件版本与运行参数

## 快速开始 | Quick Start

```bash
biopytools gffcompare -i sampleA.gff -i sampleB.gtf -o ./output
```

也可传一个目录，自动识别其中的所有 GFF/GTF 文件：

```bash
biopytools gffcompare -i ./gff_files/ -o ./output
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| GFF / GTF | 两种基因组注释格式，GTF 可看作 GFF 的一个子集，gffcompare 两者都能读 |
| query / reference | 每次比较的两个角色：query 是「被测的」，reference 是「当标准的」；A 比 B 和 B 比 A 是两次不同比较 |
| 分类 class code | gffcompare 给每个 query 转录本打的一个字母标签（如 = 完全一致、j 新亚型、u 未知间区等），表示它相对 reference 属于什么情况 |
| stats 文件 | gffcompare 输出的分类统计表，每种 class code 各有多少条 |
| tracking 文件 | 追踪表，记录每个位点上 query 与 reference 转录本的对应关系 |
| tmap / refmap | 转录本级别和位点级别的匹配清单 |
| 端部外显子范围 exon range | 判断「两个转录本端部是否对齐」允许的容差 |

## 输入 | Input

### 输入文件

2 到 5 个 GFF/GTF 文件，或包含这些文件的目录（可混用文件与目录）。支持扩展名 .gff、.gtf、.gff3 及其 .gz 压缩版本。

要求：每个文件名（去掉扩展名）互不相同，因为文件名会作为样本名出现在结果里，重复会导致校验失败。

```bash
# 三个文件的比较：会生成 3x2 = 6 个比较对
biopytools gffcompare -i a.gff -i b.gff -i c.gtf -o ./output
```

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** -i 可重复多次、每个值可以是一个文件或目录；-o 指定输出目录（默认 ./gffcompare_output）。文件总数 2 到 5 个，太少没得比，太多组合爆炸。

### 比较行为 | Comparison behavior

**通俗理解|In plain words:** 这组参数直接透传给 gffcompare，控制「怎么算一致、要不要更严格」。-e 端部外显子范围、-d 转录本起始位点分组距离，默认由 gffcompare 自定；-M/-N 丢弃单外显子转录本；-R/-Q 只比较有重叠的；--strict-match 更严格匹配；--cds-match 校验 CDS 链。一般不用动，除非你明确知道这些 gffcompare 选项的含义。

### 输出控制 | Output control

**通俗理解|In plain words:** -T 不生成 .tmap/.refmap（省空间时用）；-s/--genome-seq 提供参考 FASTA（需要 gffcompare 计算序列时才用）；-p/--cprefix 给合并 GTF 里的转录本加前缀；-V 详细模式。一般不用动。

### 运行控制 | Runtime

**通俗理解|In plain words:** --force 强制重新运行所有比较对（覆盖断点续传）；--gffcompare-path 指定 gffcompare 软件路径（默认自动找 annot 环境）。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先把所有输入两两配对（双向），再逐对跑 gffcompare，最后把所有统计合并成一张总表。

```text
输入 2-5 个 GFF/GTF 文件/目录
    |
    v
解析输入 -> 得到样本名列表 (检查唯一性, 2-5 个)
    |
    v
生成两两双向比较对 (n x (n-1) 对, 如 a_vs_b, b_vs_a)
    |
    v
逐对运行 gffcompare (断点续传: 已完成的跳过, --force 重跑)
    |  每对产出 gffcmp.stats / .tracking / .annotated.gtf / .tmap / .refmap / .loci
    v
合并所有 .stats -> 02_summary/all_stats.tsv
    |
    v
生成 software_versions.yml + 运行日志
```

## 输出 | Output

```text
output/
├── 00_pipeline_info/
│   └── software_versions.yml          # gffcompare 版本与运行参数
├── 01_gffcompare/
│   ├── sampleA_vs_sampleB/
│   │   ├── gffcmp.stats              # 分类统计 (核心)
│   │   ├── gffcmp.tracking           # 位点匹配追踪表
│   │   ├── gffcmp.annotated.gtf      # 合并后的注释 GTF
│   │   ├── gffcmp.tmap               # 转录本匹配清单 (默认生成)
│   │   ├── gffcmp.refmap             # 参考位点清单 (默认生成)
│   │   └── gffcmp.loci               # 位点汇总
│   └── ...                           # 每个比较对一个子目录
├── 02_summary/
│   └── all_stats.tsv                 # 所有比较对的统计汇总表
└── 99_logs/
    └── gffcompare.log                # 运行日志
```

比较对目录命名规则：query样本名_vs_reference样本名，例如 sampleA_vs_sampleB 表示「以 A 为 query、B 为 reference」。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 每个 .stats 文件回答「这个 query 注释相对 reference 注释，各个类别分别有多少条」；all_stats.tsv 把所有这些答案拼成一张大表，方便横向对比。

- .stats 文件：gffcompare 标准分类统计，含 class code 与条数（如 = 完全一致、j 新亚型、u 未知间区、x 反义重叠等），末尾有总数
- all_stats.tsv：每行一个「比较对 × 一个 class code」，前三列是 pair_name、query、reference，后面是该行的统计内容，可直接导入 Excel 做透视
- 注意「双向」：A_vs_B 与 B_vs_A 通常不完全对称，因为 reference 不同；要看「谁覆盖谁」得两个方向都看
- 一致性判断：完全一致（= 类）占比越高，两个注释越接近；u 类（query 独有、参考无对应）越多，说明 query 注释里多出的内容越多

## 参数选择建议 | Parameter Guidance

- 常规评估：直接传 2-3 个文件，其余全默认，看 all_stats.tsv 即可
- 只想看分类统计、不关心匹配清单：加 -T 跳过 .tmap/.refmap，省磁盘
- 明确要更严格的「端部对齐才叫一致」：用 --strict-match 并按需设 -e 范围
- 需要 gffcompare 计算序列信息（如估算序列一致性）：加 -s/--genome-seq 传参考 FASTA
- 换参数重跑：先加 --force，否则断点续传会复用旧结果

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入GFF/GTF文件或文件夹(自动识别)｜Input GFF/GTF file(s) or directory (auto-detect) |
| `--output-dir, -o` | `./gffcompare_output` |  | 输出目录｜Output directory |
| `--exon-range, -e` | — | int | 端部外显子最大允许变异范围｜Max terminal exon range |
| `--tss-distance, -d` | — | int | 转录本起始位点分组距离｜TSS grouping distance |
| `-M, --discard-single-exon-query` | — |  | 丢弃单外显子query转录本｜Discard single-exon query transcripts |
| `-N, --discard-single-exon-ref` | — |  | 丢弃单外显子reference转录本｜Discard single-exon reference transcripts |
| `-R, --ref-overlap-only` | — |  | 仅考虑与query重叠的reference｜Only consider reference overlapping query |
| `-Q, --query-overlap-only` | — |  | 仅考虑与reference重叠的query｜Only consider query overlapping reference |
| `-T, --no-tmap-refmap` | — |  | 不生成.tmap和.refmap文件｜Skip .tmap and .refmap files |
| `--strict-match` | — |  | 严格匹配模式｜Strict match mode |
| `--cds-match` | — |  | 启用CDS链匹配验证｜Enable CDS chain matching validation |
| `--genome-seq, -s` | — |  | 基因组序列路径(FASTA)｜Genome sequence path (FASTA) |
| `--cprefix, -p` | — |  | 合并GTF中转录本前缀｜Transcript prefix in combined GTF |
| `-V, --verbose-mode` | — |  | 详细处理模式｜Verbose processing mode |
| `--force, -f` | — |  | 强制重新运行｜Force re-run |
| `--gffcompare-path` | — |  | gffcompare软件路径｜gffcompare software path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入GFF/GTF文件或文件夹，支持混用(自动识别.gff/.gtf/.gff3)｜Input GFF/GTF file(s) or directory, supports mixing (auto-detect .gff/.gtf/.gff3) |
| `-o, --output-dir` | `./gffcompare_output` |  | 输出目录｜Output directory |
| `-e, --exon-range` | — | int | 端部外显子最大允许变异范围(默认100)｜Max terminal exon range (default 100) |
| `-d, --tss-distance` | — | int | 转录本起始位点分组距离(默认100)｜TSS grouping distance (default 100) |
| `-M, --discard-single-exon-query` | — | store_true | 丢弃单外显子query转录本｜Discard single-exon query transcripts |
| `-N, --discard-single-exon-ref` | — | store_true | 丢弃单外显子reference转录本｜Discard single-exon reference transcripts |
| `-R, --ref-overlap-only` | — | store_true | 仅考虑与query重叠的reference｜Only consider reference overlapping query |
| `-Q, --query-overlap-only` | — | store_true | 仅考虑与reference重叠的query｜Only consider query overlapping reference |
| `-T, --no-tmap-refmap` | — | store_true | 不生成.tmap和.refmap文件｜Skip .tmap and .refmap files |
| `--strict-match` | — | store_true | 严格匹配模式(考虑端部外显子范围)｜Strict match mode (consider -e range) |
| `--cds-match` | — | store_true | 启用CDS链匹配验证｜Enable CDS chain matching validation |
| `-s, --genome-seq` | — |  | 基因组序列路径(FASTA)｜Genome sequence path (FASTA) |
| `-p, --cprefix` | — |  | 合并GTF中转录本前缀(默认TCONS)｜Transcript prefix in combined GTF (default TCONS) |
| `-V, --verbose-mode` | — | store_true | 详细处理模式｜Verbose processing mode |
| `--force, -f` | — | store_true | 强制重新运行(覆盖断点续传)｜Force re-run (override checkpoint resume) |
| `--gffcompare-path` | — |  | gffcompare软件路径｜gffcompare software path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3（需 pyyaml）
- gffcompare，默认路径 ~/miniforge3/envs/annot/bin/gffcompare，conda 环境 annot

## 常见问题 | FAQ

**Q1：为什么我的比较对数量是 n×(n-1) 而不是 n×(n-1)/2？**
因为做的是「双向」比较：A 比 B 和 B 比 A 是两次独立比较（reference 不同）。3 个文件就是 6 对。

**Q2：能断点续传吗？**
能。每个比较对按 gffcmp.annotated.gtf 是否存在且非空判断是否已完成，已完成的自动跳过；用 --force 强制全部重跑。

**Q3：最多能传几个文件？**
最多 5 个（超过会校验失败）。因为两两双向的组合数随文件数平方增长，文件太多比较对爆炸。

**Q4：报「重复的样本名称」是什么原因？**
两个输入文件名去掉扩展名后相同。结果里靠样本名区分，必须唯一，请改文件名。

**Q5：目录输入会识别哪些文件？**
扩展名为 .gff、.gtf、.gff3 及其 .gz 压缩版本，按文件名排序去重后全部纳入比较。