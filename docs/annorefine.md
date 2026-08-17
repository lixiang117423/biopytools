# BRAKER后效应子查漏补缺 | Post-BRAKER Effector Gap-filling

**BRAKER注释后, 用miniprot比对近缘蛋白, 补回效应子多拷贝与漏注位点 | Fill gaps after BRAKER using miniprot vs. close-reference proteins: recover effector copies & missed loci**

## 功能概述 | Overview

annorefine 模块针对 BRAKER 注释后的基因组, 用 miniprot 将近缘蛋白(如已鉴定的效应子、同源物种蛋白)比对到**未 mask 原始基因组**, 找出 BRAKER 漏注的位点与被错误合并的多拷贝基因, 并补建基因模型, 输出与 BRAKER 合并后的 GFF3。

适用场景: 疫霉菌等病原中, 效应子(RxLR/CRN 等)多拷贝基因家族常落在 TE 区、被 BRAKER 漏注或错误合并, 本模块用蛋白证据补回。

## 快速开始 | Quick Start

```bash
# 基础用法(基因组 + BRAKER GFF3 + 近缘蛋白)
biopytools annorefine -g genome.fa -b braker.gff3 -p effectors.faa -o out/

# 提供 RNA-seq BAM 与 RepeatMasker .out(用于 gap 验证报告 + TE 重叠)
biopytools annorefine -g genome.fa -b braker.gff3 -p effectors.faa -o out/ \
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

(运行 `biopytools annorefine -h` 查看完整参数列表)

## 输出 | Output

- `{prefix}.gap_filled.gff3`: 补建的 gap 基因模型(gene/mRNA/exon/CDS)
- `{prefix}.merged.gff3`: 与 BRAKER 合并后的完整 GFF3(主输出)
- `{prefix}.gap_report.tsv`: gap 验证报告(蛋白证据 + RNA-seq depth + StringTie FPKM/TPM + TE 重叠/family)
- miniprot 证据 GFF3、99_logs/

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 未mask原始基因组｜Unmasked genome |
| `-s, --species` | 必填 |  | 物种名｜Species name |
| `-p, --prot-seq` | 必填 |  | 近缘蛋白(文件或目录)｜Protein file/dir |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output dir |
| `--rnaseq-dirs` | — |  | 二代RNA-seq目录(逗号分隔)｜RNA-seq dirs |
| `--isoseq` | — |  | 三代转录本(文件或目录)｜Iso-seq file/dir |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--fungus/--no-fungus` | `True` |  | 真菌模式(疫霉适用)｜Fungus mode |
| `--no-singularity` | — |  | 不用Singularity｜No singularity |
| `--skip-repeat` | — |  | 跳过repeat屏蔽｜Skip repeat masking |
| `--skip-repeat-filter` | — |  | 跳过repeat库过滤｜Skip repeat filter |
| `--skip-rescue/--no-skip-rescue` | `True` |  | 证据还原(默认关)｜Rescue (default off) |
| `--split-min-copy-coverage` | `80` | float | 保守合并判据:完整拷贝覆盖率%｜Split copy coverage |
| `--no-split` | — |  | 关闭合并拆分｜Disable split |
| `--repeat-out` | — |  | RepeatMasker .out(filling真TE排除)｜RepeatMasker out |
| `--exclude-te-gap` | — |  | 质控排除TE区gap(默认不排)｜exclude TE-overlap gaps |
| `--no-real-orf` | — |  | 关闭真实完整ORF检查(默认开)｜disable real-ORF check (default on) |
| `--no-coord-zero-overlap` | — |  | 关闭gap坐标零重叠(默认开)｜disable coord-zero-overlap (default on) |
| `--no-unique-reads` | — |  | 关闭唯一比对过滤(默认开)｜disable unique-read filter (default on) |
| `--min-unique-mapq` | `20` | int | 唯一比对MAPQ兜底阈值｜unique MAPQ fallback |
| `--min-expression-depth` | `1.0` | float | 唯一reads平均深度下限｜min unique-read depth |
| `--min-coverage-breadth` | `50.0` | float | CDS覆盖广度%下限｜min coverage breadth |
| `--no-gap-fill` | — |  | 关闭纯漏检填补(只保留合并拆分)｜disable pure gap-fill (split only) |
| `--recover-small-proteins` | — |  | 开启小蛋白回收通道(默认关, 通用)｜enable small-protein lane (default off) |
| `--small-max-cds-len` | `450` | int | 小蛋白CDS上限bp｜small max CDS len |
| `--small-min-identity` | `50.0` | float | 小蛋白放宽identity%(有表达时)｜small min identity |
| `--small-min-coverage` | `50.0` | float | 小蛋白放宽coverage%(有表达时)｜small min coverage |
| `--small-min-expression-depth` | `1.0` | float | 小蛋白表达深度下限(effector低表达可调低如0.1)｜small min expression depth |
| `--small-min-coverage-breadth` | `60.0` | float | 小蛋白CDS覆盖广度%下限｜small min coverage breadth |
| `--no-small-exclude-te` | — |  | 关闭小蛋白TE区排除(effector常在TE区可关)｜disable small-protein TE exclusion |
| `--small-strong-homology-identity` | `95.0` | float | 强同源直通identity%阈值(≥此值绕过TE/表达过滤)｜strong-homology bypass identity |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 未mask原始基因组(braker 内部 mask, filling 用未mask)｜Unmasked genome |
| `-s, --species` | 必填 |  | 物种名(braker 输出命名)｜Species name |
| `-p, --prot-seq` | 必填 |  | 近缘蛋白(文件或目录, braker+filling 共用)｜Protein file/dir |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output dir |
| `--rnaseq-dirs` | — |  | 二代RNA-seq目录(逗号分隔)｜RNA-seq dirs |
| `--isoseq` | — |  | 三代转录本(文件或目录)｜Iso-seq file/dir |
| `-t, --threads` | `12` | int | 线程数｜Threads (default 12) |
| `--fungus` | `True` |  | 真菌模式(默认开, --no-fungus 关)｜Fungus mode (default on) |
| `--singularity-image` | `~/software/singularity/braker3_devel.sif` |  | Singularity镜像｜Singularity image |
| `--no-singularity` | — | store_true | 不用Singularity｜No singularity |
| `--skip-repeat` | — |  |  |
| `--skip-repeat-filter` | — | store_true | 跳过repeat库过滤(默认开)｜Skip repeat filter |
| `--skip-rescue` | `True` |  | 跳过证据还原(默认关, --no-skip-rescue 开)｜Skip rescue (default on) |
| `--split-min-copy-coverage` | `80` | float | 保守合并判据:完整拷贝覆盖率%%｜Split copy coverage (default 80) |
| `--no-split` | — | store_true | 关闭合并拆分｜Disable merged-gene split |
| `--repeat-out` | — |  | RepeatMasker .out(默认自动找braker产物)｜RepeatMasker out |
| `--exclude-te-gap` | — | store_true | 质控排除TE区gap(默认不排)｜exclude TE-overlap gaps |
| `--gap-min-identity` | `70` | float | filling identity%%(default 70) |
| `--gap-min-coverage` | `80` | float | filling coverage%%(default 80) |
| `--no-real-orf` | — | store_true | 关闭真实完整ORF检查(ATG+stop+3倍数,默认开)｜disable real-ORF check (default on) |
| `--no-coord-zero-overlap` | — | store_true | 关闭gap坐标零重叠(默认开:与BRAKER基因坐标相交不算新基因)｜disable coord-zero-overlap (default on) |
| `--no-unique-reads` | — | store_true | 关闭唯一比对过滤(默认开:多比对reads不算表达)｜disable unique-read filter (default on) |
| `--min-unique-mapq` | `20` | int | 唯一比对MAPQ兜底阈值(samtools无-e时)｜unique MAPQ fallback (default 20) |
| `--min-expression-depth` | `1.0` | float | 唯一reads平均深度下限(>0)｜min unique-read depth (default 1.0) |
| `--min-coverage-breadth` | `50.0` | float | CDS被唯一reads覆盖广度%%下限｜min coverage breadth (default 50) |
| `--no-gap-fill` | — | store_true | 关闭纯漏检填补(只保留合并拆分)｜disable pure gap-fill (split only) |
| `--recover-small-proteins` | — | store_true | 开启小蛋白回收通道(默认关, 放宽长度找回短蛋白)｜enable small-protein lane (default off) |
| `--small-max-cds-len` | `450` | int | 小蛋白CDS上限bp(默认450=150aa)｜small max CDS len (default 450) |
| `--small-min-identity` | `50.0` | float | 小蛋白放宽identity%%(默认50, 有表达时)｜small min identity (default 50, with expr) |
| `--small-min-coverage` | `50.0` | float | 小蛋白放宽coverage%%(默认50, 有表达时)｜small min coverage (default 50, with expr) |
| `--small-min-expression-depth` | `1.0` | float | 小蛋白表达深度下限(默认1.0; effector低表达可调低如0.1)｜small min expression depth (default 1.0) |
| `--small-min-coverage-breadth` | `60.0` | float | 小蛋白CDS覆盖广度%%下限(默认60)｜small min coverage breadth (default 60) |
| `--no-small-exclude-te` | — | store_true | 关闭小蛋白TE区排除(默认排; effector常在TE区可关)｜disable small-protein TE exclusion (default on) |
| `--small-strong-homology-identity` | `95.0` | float | 强同源直通identity%%阈值(默认95, ≥此值绕过TE/表达过滤)｜strong-homology bypass identity (default 95) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **miniprot**: 蛋白→基因组比对
- **samtools**: RNA-seq depth 统计(gap 报告)
- **StringTie**: FPKM/TPM 定量(gap 报告)
- **gffread**: 序列提取(按需)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
