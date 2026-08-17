# gffcompare | GFF/GTF 两两双向比较分析

**使用 gffcompare 对多个 GFF/GTF 注释文件执行两两双向比较并汇总分类统计 | Pairwise bidirectional comparison of GFF/GTF annotations using gffcompare**

## 功能概述 | Overview

本模块包装 [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)（StringTie2 套件成员），对一组 GFF/GTF 文件两两双向运行比较：每一对 `(A, B)` 都跑两次（A vs B、B vs A），输出每个转录本相对另一注释集的分类（`=`, `c`, `k`, `j`, `e`, `o` 等）、坐标重叠、外显子匹配等信息，并自动汇总所有 `.stats` 文件为单一 TSV。

常用于评估不同基因预测/注释流程（如 BRAKER、EVidenceModeler、StringTie 等）的一致性，或与参考注释（如 Ensembl/NCBI）做比对，定位新增、缺失或差异转录本。

## 快速开始 | Quick Start

```bash
# 比较两个注释文件（双向各跑一次）
biopytools gffcompare -i sampleA.gff -i sampleB.gtf -o ./output

# 比较整个目录下的所有 GFF/GTF
biopytools gffcompare -i ./gff_files/ -o ./output
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入 GFF/GTF 文件或目录（可重复；目录会自动识别） |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./gffcompare_output` | 输出目录 |
| `-e, --exon-range` | - | 端部外显子最大允许变异范围 |
| `-d, --tss-distance` | - | 转录本起始位点（TSS）分组距离 |
| `-s, --genome-seq` | - | 基因组序列路径（FASTA，用于 CDS 验证等） |
| `-p, --cprefix` | - | 合并 GTF 中转录本前缀（默认 TCONS） |
| `--gffcompare-path` | - | gffcompare 软件路径（可置于 conda 环境） |
| `-f, --force` | 关闭 | 强制重新运行 |
| `-V, --verbose-mode` | 关闭 | 详细处理模式 |

### 比较模式选项 | Comparison Mode Flags

| 参数 | 描述 |
|------|------|
| `-M, --discard-single-exon-query` | 丢弃单外显子 query 转录本 |
| `-N, --discard-single-exon-ref` | 丢弃单外显子 reference 转录本 |
| `-R, --ref-overlap-only` | 仅考虑与 query 重叠的 reference |
| `-Q, --query-overlap-only` | 仅考虑与 reference 重叠的 query |
| `-T, --no-tmap-refmap` | 不生成 `.tmap` 和 `.refmap` 文件 |
| `--strict-match` | 严格匹配模式 |
| `--cds-match` | 启用 CDS 链匹配验证 |

（运行 `biopytools gffcompare -h` 查看完整参数列表）

## 输出 | Output

```
gffcompare_output/
├── 01_gffcompare/
│   ├── A_vs_B/
│   │   ├── gffcmp.stats         # 分类计数统计
│   │   ├── gffcmp.tracking      # 转录本对应关系
│   │   ├── gffcmp.annotated.gtf # 带分类标签的 query GTF
│   │   ├── gffcmp.tmap          # 转录本映射
│   │   └── gffcmp.refmap        # 参考 vs query 映射
│   └── B_vs_A/                  # 反向比较结果
└── 02_summary/
    └── all_stats.tsv            # 所有 .stats 汇总（含 pair/query/ref 列）
```

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

- [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)（必需，StringTie2 套件成员）
- 可选：conda 环境（模块支持通过 `--gffcompare-path` 指定 conda 环境中的 gffcompare）

## 引用 | Citation

- Pertea M, Kim D, Pertea GM, Leek JT, Salzberg SL. Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. *Nature Protocols*, 2016, 11(9): 1650-1667. doi:10.1038/nprot.2016.095
- Kovaka S, Zimin AV, Pertea GM, et al. Transcriptome assembly from long-read RNA-seq alignments with StringTie2. *Genome Biology*, 2019, 20: 278. doi:10.1186/s13059-019-1910-1

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [gffcompare 官方页面](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
- [gffcompare 分类代码说明](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#class_codes)
