# Minimap2 比对与未比对区间提取 | Minimap2 Alignment and Unmapped Region Extraction

**使用 Minimap2 进行全基因组比对，并自动提取未比对到的区间和序列 | Perform whole-genome alignment with Minimap2 and automatically extract unmapped regions and sequences.**

## 功能概述 | Overview

- 调用 Minimap2 完成查询基因组到目标基因组的全基因组比对，输出 PAF
- 解析 PAF，按匹配长度与 tp 类型 (primary/secondary) 过滤有效比对
- 扫描每条 query 序列，找出未被有效比对覆盖的区间，生成 BED
- 调用 seqkit 根据 BED 提取未比对序列 FASTA
- 生成详细的总结报告，按 query 统计未比对区间数量与总长度

## 快速开始 | Quick Start

```bash
# 默认asm5预设，适合近缘基因组
biopytools minimap2 -t target_genome.fasta -q query_genome.fasta -o results/

# 使用map-ont预设把ONT reads比对到基因组
biopytools minimap2 -t genome.fasta -q reads.fasta -o ont_align -x map-ont
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-t, --target` | 目标基因组文件路径 |
| `-q, --query` | 查询基因组文件路径 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./minimap2_output` | 输出目录 |
| `-x, --preset` | `asm5` | Minimap2 预设：`asm5`/`asm10`/`asm20`/`map-ont`/`map-pb` |
| `-p, --threads` | `12` | 线程数 |
| `-m, --min-match` | `1000` | 最小匹配长度阈值 (PAF 中 number_match) |
| `-u, --min-unmapped` | `1000` | 最小未比对区间长度阈值 |
| `--tp-type` | `P` | 保留的比对类型：`S`(secondary)/`P`(primary)/`SP`(both) |
| `-M, --minimap2-path` | `minimap2` | minimap2 可执行文件路径 |
| `-S, --seqkit-path` | `seqkit` | seqkit 可执行文件路径 |

（运行 `biopytools minimap2 -h` 查看完整参数列表）

## 输出 | Output

输出目录下生成（`{query_base}` 为查询文件主名）：

- `{query_base}_alignment.paf`：Minimap2 原始 PAF 比对结果
- `{query_base}_unmapped.bed`：未比对区间 BED 文件
- `{query_base}_unmapped.fa`：未比对序列 FASTA（由 seqkit 提取）
- `minimap2_summary.txt`：包含参数、输出文件与各 query 统计的总结报告
- 运行日志文件

## 依赖 | Dependencies

- minimap2
- seqkit
- Python 包：pandas（用于 PAF 解析）

## 引用 | Citation

- Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 2018.
- Shen W. et al. seqkit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLoS ONE*, 2016.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
