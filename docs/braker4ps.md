# braker + 查漏补缺端到端 | braker + ps-gene-anno End-to-End

**端到端编排器: 阶段1 BRAKER注释 + 阶段2 ps-gene-anno效应子查漏补缺 | End-to-end orchestrator: BRAKER annotation (Phase 1) + ps-gene-anno gap-filling (Phase 2)**

## 功能概述 | Overview

braker4ps 是一个端到端编排器, 仅通过 import 调用 braker 与 ps_gene_anno 两个模块(不改其源码), 串成完整流程:

- **阶段 1**: 用 BRAKR3 做基因组注释(RNA-seq + 蛋白证据)
- **阶段 2**: 调 ps_gene_anno 对 BRAKER 结果做效应子查漏补缺——自动探测阶段 1 产出的 RepeatMasker `.out` 与 RNA-seq BAM, 供 gap 验证报告使用

统一日志, 两阶段输出分层存放。

## 快速开始 | Quick Start

```bash
# 基础用法
biopytools braker4ps -g genome.fa -s psojae -p prot.fa --rnaseq-dirs r1,r2 -o out/

# 关闭效应子救援相关的合并拆分
biopytools braker4ps -g genome.fa -s psojae -p prot.fa --rnaseq-dirs r1,r2 -o out/ --no-split
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-g, --genome` | 基因组 FASTA |
| `-s, --species` | 物种名 |
| `-p, --prot-seq` | 蛋白证据 |
| `-o, --output-dir` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--rnaseq-dirs` | — | 二代 RNA-seq 目录(逗号分隔) |
| `--isoseq` | — | 三代转录本(文件或目录) |
| `-t, --threads` | `12` | 线程数 |
| `--fungus` / `--no-fungus` | `True` | BRAKER fungus 模式 |
| `--singularity-image` | 自动 | singularity 镜像 |
| `--no-singularity` | `False` | 不用 singularity |
| `--skip-repeat` | `False` | 跳过 repeat 屏蔽 |
| `--skip-repeat-filter` | `False` | 跳过 repeat 库过滤 |
| `--skip-rescue` | `True` | 跳过 braker 效应子救援 |
| `--repeat-out` | 自动探测 | RepeatMasker .out(阶段 2 用) |
| `--exclude-te-gap` | `False` | 阶段 2 质控排除 TE 区 gap |
| `--no-split` | `False` | 关闭阶段 2 合并基因拆分 |
| `--split-min-copy-coverage` | `80` | 合并判定拷贝完整覆盖 % |
| `--gap-min-identity` | `70` | 阶段 2 filling identity % |
| `--gap-min-coverage` | `80` | 阶段 2 filling coverage % |

(运行 `biopytools braker4ps -h` 查看完整参数列表)

## 输出 | Output

- **阶段 1**(BRAKER): 注释 GFF3、蛋白、RepeatMasker `.out`、RNA-seq BAM 等
- **阶段 2**(`05_gap_filling/`): `gap_filled.gff3`、`merged.gff3`、`gap_report.tsv`
- 流程元数据与日志

## 依赖 | Dependencies

- **BRAKER3**: 阶段 1 基因组注释
- **ps_gene_anno 的依赖**: miniprot / samtools / StringTie / gffread(阶段 2)
- (可选) **singularity**: BRAKER 运行容器

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
