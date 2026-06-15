# 叶绿体基因组组装 | Plastome Assembly

**基于 GetOrganelle 从二代全基因组测序数据中自动组装叶绿体基因组 | Assemble plastome from WGS data using GetOrganelle**

## 功能概述 | Overview

本模块包装了 GetOrganelle，可从常规 Illumina 双端 WGS 数据中通过 k-mer 扩展算法特异性地"钓取"并组装细胞器基因组（叶绿体、线粒体、核糖体 DNA 等）。支持批量处理多个样品目录，自动检测 R1/R2 文件并调度组装。

GetOrganelle 通过种子 k-mer 在 reads pool 中迭代延伸、去冗余、环化，对参考基因组依赖低，特别适合植物叶绿体基因组（plastome）的高质量参考构建。默认目标为 `embplant_pt`（植物叶绿体），也可切换到线粒体等其它细胞器类型。

## 快速开始 | Quick Start

```bash
# 基本用法：批量处理目录中的样品
biopytools get-plastome -i fastq_folder -o plastome_output

# 单样品模式（隐式参数，默认批量）
biopytools get-plastome -i ./sample_reads -o ./plastome -p my_sample

# 切换为植物线粒体
biopytools get-plastome -i ./reads -o ./mt_output --organelle-type embplant_mt
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入目录（包含 reads 文件，批量模式时每个子目录为一个样品）|

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./plastome_output` | 输出目录 |
| `-p, --prefix` | - | 输出前缀（单样品模式）|
| `--organelle-type` | `embplant_pt` | 细胞器类型（embplant_pt/embplant_mt/embplant_nr/other_pt/animal_mt/fungus_mt/fungus_nr）|
| `-R, --max-rounds` | `15` | 最大扩展轮数 |
| `-k, --kmer-list` | `21,45,65,85,105` | K-mer 列表（逗号分隔）|
| `-t, --threads` | `12` | 线程数 |
| `--getorganelle-path` | `~/miniforge3/envs/getorganelle_v.1.7.71/bin/get_organelle_from_reads.py` | GetOrganelle 脚本路径 |
| `-v, --verbose` | 关 | 详细输出 |
| `--log-file` | - | 日志文件路径 |

（运行 `biopytools get-plastome -h` 查看完整参数列表）

## 输出 | Output

```
plastome_output/
├── <sample_name>/
│   ├── GetOrganelle_<type>/   # GetOrganelle 原始工作目录
│   │   ├── assembled_graph/   # 组装图（GFA）
│   │   └── filtered_fastas/   # 过滤后的 FASTA
│   └── <sample>.fa            # 最终环化叶绿体序列（如有）
└── plastome_assembly.log      # 流程日志
```

## 依赖 | Dependencies

- GetOrganelle（默认通过 conda 环境调用，`--getorganelle-path` 可指定）
- Bowtie2（GetOrganelle 内部依赖）
- SPAdes（GetOrganelle 内部依赖）

## 引用 | Citation

- Jin, J.-J. et al. GetOrganelle: a fast and versatile toolkit for accurate de novo assembly of organelle genomes. *Genome Biology* 21, 241 (2020).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [GetOrganelle 官方仓库](https://github.com/Kinggerm/GetOrganelle)
