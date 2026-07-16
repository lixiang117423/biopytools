# BRAKER后效应子查漏补缺 | Post-BRAKER Effector Gap-filling

**BRAKER注释后, 用miniprot比对近缘蛋白, 补回效应子多拷贝与漏注位点 | Fill gaps after BRAKER using miniprot vs. close-reference proteins: recover effector copies & missed loci**

## 功能概述 | Overview

ps_gene_anno 模块针对 BRAKER 注释后的基因组, 用 miniprot 将近缘蛋白(如已鉴定的效应子、同源物种蛋白)比对到**未 mask 原始基因组**, 找出 BRAKER 漏注的位点与被错误合并的多拷贝基因, 并补建基因模型, 输出与 BRAKER 合并后的 GFF3。

适用场景: 疫霉菌等病原中, 效应子(RxLR/CRN 等)多拷贝基因家族常落在 TE 区、被 BRAKER 漏注或错误合并, 本模块用蛋白证据补回。

## 快速开始 | Quick Start

```bash
# 基础用法(基因组 + BRAKER GFF3 + 近缘蛋白)
biopytools ps-gene-anno -g genome.fa -b braker.gff3 -p effectors.faa -o out/

# 提供 RNA-seq BAM 与 RepeatMasker .out(用于 gap 验证报告 + TE 重叠)
biopytools ps-gene-anno -g genome.fa -b braker.gff3 -p effectors.faa -o out/ \
  --rnaseq-bam rnaseq.bam --repeat-out rm.out
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-g, --genome` | 未 mask 原始基因组 |
| `-b, --braker-gff3` | BRAKER 输出 GFF3 |
| `-p, --prot-seq` | 近缘蛋白(证据) |
| `-o, --output-dir` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `--prefix` | genome 文件名 | 输出前缀 |
| `--rnaseq-bam` | — | RNA-seq BAM(逗号分隔; gap 报告用) |
| `--repeat-out` | — | RepeatMasker .out(TE 重叠分析) |
| `--gap-min-identity` | `70` | miniprot 命中最低 identity % |
| `--gap-min-coverage` | `80` | 最低覆盖蛋白比例 % |
| `--gap-min-cds-len` | `100` | 最小 CDS 长度(bp) |
| `--overlap-cutoff` | `0` | 漏检判定: 与 BRAKER CDS 重叠 < 此 % 才算漏检(0=零重叠) |
| `--exclude-te-gap` | `False` | 质控排除 TE 区 gap(默认不排, 疫霉效应子常在 TE 区) |
| `--no-split` | `False` | 关闭错误合并基因拆分 |
| `--split-min-hits` | `2` | 合并判定: gene 内完整独立拷贝数下限 |
| `--split-min-copy-coverage` | `80` | 拷贝完整性覆盖 % |
| `--skip-merge` | `False` | 跳过与 BRAKER 合并(只输出 gap_filled) |

(运行 `biopytools ps-gene-anno -h` 查看完整参数列表)

## 输出 | Output

- `{prefix}.gap_filled.gff3`: 补建的 gap 基因模型(gene/mRNA/exon/CDS)
- `{prefix}.merged.gff3`: 与 BRAKER 合并后的完整 GFF3(主输出)
- `{prefix}.gap_report.tsv`: gap 验证报告(蛋白证据 + RNA-seq depth + StringTie FPKM/TPM + TE 重叠/family)
- miniprot 证据 GFF3、99_logs/

## 依赖 | Dependencies

- **miniprot**: 蛋白→基因组比对
- **samtools**: RNA-seq depth 统计(gap 报告)
- **StringTie**: FPKM/TPM 定量(gap 报告)
- **gffread**: 序列提取(按需)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
