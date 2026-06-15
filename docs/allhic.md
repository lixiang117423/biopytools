# ALLHiC 染色体挂载 | ALLHiC Chromosome Scaffolding

**基于 Hi-C 数据的染色体级别基因组挂载流程 | Hi-C-assisted chromosome-level genome scaffolding**

## 功能概述 | Overview

ALLHiC 是专为多倍体/复杂基因组设计的 Hi-C 辅助挂载工具，能够在已知染色体数目的前提下，通过修剪、分区、优化、救援等步骤将 contigs 组装成染色体级别的 pseudomolecule。本模块封装了 ALLHiC v5.4 的完整流程（含 Asmkit/JBAT 后处理）。

流程分 8 个步骤：1) Hi-C 比对过滤；2) 等位基因检测；3) contig 修剪；4) 染色体分区；5) 提取矩阵；6) Rescue 救援；7) 构建染色体 FASTA；8) 热图绘制与 JBAT 报告生成。每个步骤均可单独跳过，便于断点续跑和调参。适用于异源多倍体作物（如小麦、棉花、油菜）和大型二倍体基因组。

## 快速开始 | Quick Start

```bash
# 基本用法：12 条染色体的挂载
biopytools allhic -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -k 12

# 指定酶切位点和更多线程
biopytools allhic -r asm.fa -1 R1.fq.gz -2 R2.fq.gz -k 21 \
    -e HindIII -t 64 -o ./allhic_results

# 已有比对结果，跳过 mapping 步骤
biopytools allhic -r asm.fa -1 R1.fq.gz -2 R2.fq.gz -k 12 --skip-mapping
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-r, --reference` | 参考基因组（contig 级别 FASTA）|
| `-1, --read1` | Hi-C 读段 R1 文件 |
| `-2, --read2` | Hi-C 读段 R2 文件 |
| `-k, --chr-num` | 目标染色体数目 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-e, --enzyme` | `GATC` | 酶切位点 motif |
| `-t, --threads` | `12` | 线程数 |
| `-w, --workdir` | `./allhic_output` | 工作目录 |
| `--mapq-step1` | `1` | 步骤 1 比对质量阈值 |
| `--bin-size` | `500k` | Bin 大小 |
| `--min-bin-size` | `50k` | 最小 Bin 大小 |
| `--skip-mapping` | 关 | 跳过步骤 1（Hi-C 比对）|
| `--skip-allele` | 关 | 跳过步骤 1.5（等位基因检测）|
| `--skip-prune` | 关 | 跳过步骤 2（修剪）|
| `--skip-partition` | 关 | 跳过步骤 3（分区）|
| `--skip-extract` | 关 | 跳过步骤 3.5（提取矩阵）|
| `--skip-rescue` | 关 | 跳过步骤 4（救援）|
| `--skip-optimize` | 关 | 跳过步骤 5（优化）|
| `--skip-build` | 关 | 跳过步骤 6（构建 FASTA）|
| `--skip-plot` | 关 | 跳过步骤 7（绘制热图）|
| `--skip-asmkit` | 关 | 跳过步骤 8（JBAT 生成）|
| `--diagnose` | 关 | 诊断模式 |
| `-v, --verbose` | 关 | 详细输出 |

（运行 `biopytools allhic -h` 查看完整参数列表）

## 输出 | Output

```
allhic_output/
├── 01_mapping/          # Hi-C 比对与过滤
├── 02_allele/           # 等位基因检测结果
├── 03_prune/            # contig 修剪
├── 04_partition/        # 染色体分区（clusters）
├── 05_extract/          # 提取的 contact 矩阵
├── 06_rescue/           # 救援结果
├── 07_optimize/         # 优化后的 ordering
├── 08_build/            # 染色体级别 FASTA
├── 09_plot/             # Hi-C 热图
├── 10_asmkit/           # JBAT/Asmkit 报告
└── pipeline.log         # 流程日志
```

## 依赖 | Dependencies

- ALLHiC 工具包（ALLHiC_pipeline、ALLHiC_rescue、ALLHiC_optimize 等）
- SAMtools（BAM 处理）
- BWA（Hi-C 比对，ALLHiC 内部调用）
- Juicer Tools / Asmkit（JBAT 报告生成，步骤 8）
- Python + matplotlib（热图绘制）

## 引用 | Citation

- Zhang, X. et al. Allele-aware genome assembly of the allopolyploid common wheat. *Nature Genetics* 50, 952-961 (2018).
- Zhang, L. et al. ALLHiC: scaffolding large ploidy genomes using Hi-C data. *Nature Methods* 16, 1325-1326 (2019).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [ALLHiC 官方仓库](https://github.com/tangerzhang/ALLHiC)
