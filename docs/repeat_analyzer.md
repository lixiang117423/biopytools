# 重复序列分析 | Repeat Sequence Analyzer

**整合 RepeatModeler + LTR 流水线 + RepeatMasker + TEsorter 的一站式重复序列注释 | All-in-one repeat annotation pipeline integrating RepeatModeler, LTR pipeline, RepeatMasker, and TEsorter.**

## 功能概述 | Overview

`repeat_analyzer` 模块对一个基因组的重复序列进行完整注释：先用 `RepeatModeler` 构建从头重复库；并行运行 `LTR_FINDER` + `LTRharvest` 并由 `LTR_retriever` 整合得到 LTR 库；将两个库合并后，用 `RepeatMasker` 在基因组上注释并定量，最后用 `TEsorter` 对库中的完整转座子进行结构分类。流程支持跳过任意阶段，可自定义所有工具路径，自动检测 PATH 中的可执行文件，并输出基因组统计与完整的分析摘要报告。

## 快速开始 | Quick Start

```bash
# 完整流程
biopytools repeat-analyzer -i genome.fa -o repeat_results -t 64

# 仅基于已有的 RepeatModeler 库 + RepeatMasker + TEsorter（跳过 LTR）
biopytools repeat-analyzer -i genome.fa -o results --skip-ltr

# 仅运行 RepeatMasker + TEsorter（使用外部库）
biopytools repeat-analyzer -i genome.fa -o results --skip-modeler --skip-ltr

# 指定工具路径
biopytools repeat-analyzer -i genome.fa -o results \
    --repeatmasker-path ~/.local/bin/RepeatMasker \
    --tesorter-path ~/.local/bin/TEsorter -t 32
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入基因组 FASTA 文件路径 |
| `-o, --output` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `88` | 线程数 |
| `--skip-modeler` | off | 跳过 RepeatModeler 步骤 |
| `--skip-ltr` | off | 跳过 LTR 分析步骤（LTR_FINDER/LTRharvest/LTR_retriever） |
| `--repeatmodeler-path` | `RepeatModeler` | RepeatModeler 路径 |
| `--ltr-finder-path` | `ltr_finder` | LTR_FINDER 路径 |
| `--ltrharvest-path` | `gt ltrharvest` | LTRharvest 路径 |
| `--ltr-retriever-path` | `LTR_retriever` | LTR_retriever 路径 |
| `--repeatmasker-path` | `RepeatMasker` | RepeatMasker 路径 |
| `--tesorter-path` | `TEsorter` | TEsorter 路径 |

> 注：CLI 默认线程数为 `88`，如果通过 `-t` 指定非 88 值才会传给底层。工具路径只要在 PATH 中即可使用默认值。

（运行 `biopytools repeat-analyzer -h` 查看完整参数列表）

## 流程步骤 | Pipeline Steps

1. **基因组统计** | Genome statistics — 序列数、总长、N50、GC 含量
2. **依赖检查** | Dependency check — 自动检测所有外部工具
3. **RepeatModeler** | 从头发现重复序列（可 `--skip-modeler`）
4. **LTR 流水线** | LTR_FINDER + LTRharvest → LTR_retriever（可 `--skip-ltr`）
5. **合并重复库** | Combine libraries — RepeatModeler 库 + LTR 库
6. **RepeatMasker** | 基因组注释与定量
7. **TEsorter** | 完整转座子结构分类（超家族/order 注释）
8. **生成摘要报告** | `repeat_analysis_summary.txt`

## 输出 | Output

```
repeat_results/
├── repeat_analysis_summary.txt   # 汇总报告（基因组统计、库文件、输出文件清单）
├── *-families.fa / *.lib         # RepeatModeler/LTR_retriever 库
├── *_combined.fa                 # 合并后的重复库
├── *.out / *.gff / *.tbl         # RepeatMasker 注释结果
├── *.cls / *.cls.tsv             # TEsorter 分类结果
└── 99_logs/                      # 运行日志
```

## 依赖 | Dependencies

- `RepeatModeler` + `RECON`/`RepeatScout`（RepeatModeler 内部依赖）+ `Tandem Repeats Finder`
- `LTR_FINDER`、`LTRharvest` (GenomeTools)、`LTR_retriever`（含 `CD-HIT`）
- `RepeatMasker`（含 `rmblast`/`cross_match` 之一 + RepBase 或 Dfam）
- `TEsorter`（含 `HMMER`）
- `BEDTools`、`seqkit`（部分子步骤使用）

## 引用 | Citation

- Flynn JM et al. RepeatModeler2 for automated genomic discovery of transposable element families. PNAS. 2020.
- Ou S. & Jiang N. LTR_retriever: A highly accurate and sensitive program for identification of long terminal repeat retrotransposons. Plant Physiology. 2018.
- Smit AFA & Green P. RepeatMasker. (RepeatMasker)
- Zhang RG et al. TEsorter: an accurate and fast method to classify LTR-retrotransposons in plant genomes. Horticulture Research. 2022.

## 相关链接 | References

- [RepeatMasker](https://www.repeatmasker.org/)
- [TEsorter GitHub](https://github.com/zhangrengang/TEsorter)
- [项目主页](https://github.com/lixiang117423/biopytools)
