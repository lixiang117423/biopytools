# BUSCO 基因组质量评估工具

**Benchmarking Universal Single-Copy Orthologs 质量评估分析 | BUSCO Quality Assessment Analysis**

## 📖 功能概述 | Overview

BUSCO (Benchmarking Universal Single-Copy Orthologs) 是评估基因组组装和注释质量的权威工具。通过搜索物种特异性单拷贝直系同源基因，BUSCO可以评估基因组完整性、注释质量以及与其他物种的比较分析。本工具基于BUSCO软件封装，提供批处理、结果汇总和可视化功能。

## ✨ 主要特性 | Key Features

- **🎯 多种分析模式**: 支持基因组、转录组、蛋白质序列三种模式
- **📊 批量处理**: 支持目录批量分析，自动汇总结果
- **🧬 多谱系数据库**: 支持所有BUSCO官方谱系数据库
- **⚙️ 灵活配置**: 支持多种基因预测工具（Augustus/Metaeuk/Miniprot）
- **🔄 自动谱系选择**: 可根据序列自动选择最佳谱系数据库
- **📈 结果汇总**: 自动生成汇总表格和统计报告
- **💾 多格式输出**: 支持TXT、CSV、XLSX格式结果输出
- **🚀 高性能**: 多线程并行分析，充分利用计算资源

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基因组质量评估（推荐）
biopytools busco \
    -i genome.fa \
    -l eukaryota_odb12 \
    -t 16 \
    -o ./busco_results
```

### 批量分析 | Batch Analysis

```bash
# 批量评估多个基因组
biopytools busco \
    -i ./genomes/ \
    -l brassicales_odb12 \
    -m genome \
    -t 24 \
    -o ./batch_results
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input` | 输入文件或目录路径 | `-i genome.fa` |
| `-l, --lineage` | BUSCO数据库谱系名称 | `-l eukaryota_odb12` |

### 基本参数 | Basic Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./busco_output` | 📁 输出目录路径 |
| `-m, --mode` | `genome` | 🔬 分析模式 (genome/transcriptome/proteins) |
| `-t, --threads` | `12` | ⚙️ CPU线程数 |
| `--sample-suffix` | `*.fa` | 🏷️ 样本名提取后缀模式 |
| `--output-format` | `txt` | 📄 输出文件格式 (txt/csv/xlsx) |

### 分析模式 | Analysis Modes

| 模式 | 简写 | 描述 | 适用场景 |
|------|------|------|----------|
| `genome` | `geno` | 基因组模式 | 基因组组装质量评估 |
| `transcriptome` | `tran` | 转录组模式 | 转录组组装完整性评估 |
| `proteins` | `prot` | 蛋白质模式 | 蛋白质序列完整性评估 |

### 数据库参数 | Database Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--datasets-version` | `odb12` | 📚 数据集版本 (odb10/odb12) |
| `--download-path` | `None` | 💾 数据集下载路径 |
| `--offline` | `False` | 🔌 离线模式（使用本地数据库） |

### 自动谱系选择 | Auto Lineage Selection

| 参数 | 描述 |
|------|------|
| `--auto-lineage` | 🔍 自动选择最佳谱系数据库 |
| `--auto-lineage-euk` | 🧬 自动选择真核生物谱系 |
| `--auto-lineage-prok` | 🦠 自动选择原核生物谱系 |

### 基因预测工具 | Gene Prediction Tools

| 参数 | 描述 |
|------|------|
| `--augustus` | 🧬 使用Augustus基因预测器 |
| `--augustus-parameters` | ⚙️ Augustus额外参数 |
| `--augustus-species` | 🏷️ Augustus物种名称 |
| `--metaeuk` | 🔬 使用Metaeuk基因预测器 |
| `--metaeuk-parameters` | ⚙️ Metaeuk额外参数 |
| `--metaeuk-rerun-parameters` | ⚙️ Metaeuk重运行参数 |
| `--miniprot` | 🧪 使用Miniprot基因预测器 |

### 性能参数 | Performance Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-e, --evalue` | `1e-3` | 📊 BLAST E值阈值 |
| `--limit` | `3` | 🔢 候选区域数量限制 |
| `--contig-break` | `10` | 🔗 Contig打断的N数量 |
| `--long` | `False` | ⏱️ 启用Augustus长模式优化 |

### 其他选项 | Other Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-f, --force` | `False` | 🔴 强制重写现有文件 |
| `-r, --restart` | `False` | 🔄 重启未完成的分析 |
| `--skip-bbtools` | `False` | ⏭️ 跳过BBTools统计 |
| `--scaffold-composition` | `False` | 📊 生成scaffold组成文件 |
| `--tar` | `False` | 📦 压缩子目录 |
| `-q, --quiet` | `False` | 🔇 静默模式 |
| `--busco-path` | `busco` | 🛠️ BUSCO软件路径 |

## 📚 BUSCO谱系数据库 | Lineage Datasets

### 常用谱系数据库 | Common Lineage Datasets

| 谱系名称 | 物种范围 | 基因数量 |
|----------|----------|----------|
| **eukaryota_odb12** | 真核生物通用 | 255 |
| **bacteria_odb12** | 细菌通用 | 148 |
| **archaea_odb12** | 古菌通用 | 162 |
| **metazoa_odb12** | 动物界 | 954 |
| **vertebrata_odb12** | 脊椎动物 | 5,306 |
| **arthropoda_odb12** | 节肢动物 | 1,066 |
| **insecta_odb12** | 昆虫纲 | 1,367 |
| **plantae_odb12** | 植物界 | 1,614 |
| **fungi_odb12** | 真菌界 | 758 |

### 物种特异性谱系 | Species-Specific Lineages

```bash
# 十字花科
brassicales_odb12

# 禾本科
poales_odb12

# 茄科
solanales_odb12

# 查看所有可用谱系
busco --list-datasets
```

## 💡 使用示例 | Usage Examples

### 示例1：植物基因组评估 | Example 1: Plant Genome Assessment

```bash
# 使用植物特异性谱系评估
biopytools busco \
    -i plant_genome.fa \
    -l plantae_odb12 \
    -m genome \
    -t 24 \
    -o ./plant_busco
```

### 示例2：自动选择谱系 | Example 2: Auto Lineage Selection

```bash
# 让BUSCO自动选择最佳谱系
biopytools busco \
    -i unknown_genome.fa \
    --auto-lineage \
    -t 16 \
    -o ./auto_results
```

### 示例3：转录组完整性评估 | Example 3: Transcriptome Completeness

```bash
# 评估转录组组装完整性
biopytools busco \
    -i transcriptome.fa \
    -l metazoa_odb12 \
    -m transcriptome \
    -t 12 \
    -o ./rna_busco
```

### 示例4：使用Augustus基因预测 | Example 4: Use Augustus Gene Prediction

```bash
# 使用Augustus预测工具，指定物种
biopytools busco \
    -i genome.fa \
    -l eukaryota_odb12 \
    --augustus \
    --augustus-species arabidopsis \
    -t 24 \
    -o ./augustus_results
```

### 示例5：批量分析并导出Excel | Example 5: Batch Analysis with Excel Export

```bash
# 批量分析多个样本，导出Excel格式结果
biopytools busco \
    -i ./genomes/*.fa \
    -l brassicales_odb12 \
    -t 32 \
    --output-format xlsx \
    -o ./batch_excel_results
```

### 示例6：离线模式分析 | Example 6: Offline Mode Analysis

```bash
# 使用本地数据库进行离线分析
biopytools busco \
    -i genome.fa \
    -l eukaryota_odb12 \
    --offline \
    --download-path /path/to/busco_downloads \
    -t 16 \
    -o ./offline_results
```

## 📊 输出结果 | Output Results

### 输出目录结构 | Output Directory Structure

```
busco_output/
├── sample1_busco/              # 单个样本结果目录|Individual sample result directory
│   ├── run_busco.tsv          # BUSCO主结果|Main BUSCO results
│   ├── short_summary.txt      # 简要汇总|Short summary
│   ├── full_summary.txt       # 完整汇总|Full summary
│   ├── busco_sequences.fasta  # BUSCO基因序列|BUSCO gene sequences
│   └── missing_busco_list.txt # 缺失基因列表|Missing genes list
├── busco_results.txt          # 汇总结果表格|Summary results table
├── busco_summary.txt          # 统计汇总报告|Statistical summary report
└── busco_analysis.log         # 运行日志|Run log
```

### BUSCO指标说明 | BUSCO Metrics

| 指标 | 描述 | 计算方式 |
|------|------|----------|
| **Complete (C)** | 完整BUSCO基因 | 单拷贝+ duplicated |
| **Single Copy (S)** | 单拷贝完整基因 | 1:1直系同源 |
| **Duplicated (D)** | 重复完整基因 | 多拷贝直系同源 |
| **Fragmented (F)** | 片段化基因 | 部分匹配 |
| **Missing (M)** | 缺失基因 | 未检测到 |

### 质量评估标准 | Quality Assessment Standards

**基因组完整性评级**：

| 级别 | 完整性 | 说明 |
|------|--------|------|
| **A** | >95% | 完整性高，适合高质量分析 |
| **B** | 90-95% | 完整性良好，基本满足需求 |
| **C** | 80-90% | 完整性中等，需注意质量 |
| **D** | <80% | 完整性较差，建议改善组装 |

**组装质量评估**：
- Single Copy比例高：基因组重复区域少，组装质量好
- Duplicated比例高：可能存在heterozygosity或组装冗余
- Fragmented比例高：基因组contiguity较差
- Missing比例高：基因组覆盖不完整

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **BUSCO** (版本 5.0 或更新)
  - 安装方法: `conda install -c bioconda busco`
  - 或: `pip install busco`

- **核心依赖软件**:
  - [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (>=2.12) - 序列比对
  - [HMMER3](http://hmmer.org/) (>=3.3) - 结构域搜索
  - [Augustus](http://bioinf.uni-greifswald.de/augustus/) (>=3.5) - 基因预测
  - [Metaeuk](https://github.com/soedinglab/MMseqs2-app) (可选) - 宏基因组基因预测
  - [Miniprot](https://github.com/mengyao/Miniprot) (可选) - 蛋白质比对
  - [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (可选) - 统计分析

### Python依赖 | Python Dependencies

- Python >= 3.7
- pandas
- numpy

### 环境配置 | Environment Setup

```bash
# 创建Conda环境并安装BUSCO及依赖
conda create -n busco_env python=3.9
conda activate busco_env

# 安装BUSCO
conda install -c bioconda busco

# 验证安装
busco --version

# 查看可用谱系
busco --list-datasets
```

## ⚠️ 注意事项 | Important Notes

1. **谱系选择**: 选择与待评估物种进化关系最近的谱系数据库
2. **内存需求**: BUSCO分析内存消耗较大，建议至少50GB RAM
3. **磁盘空间**: 临时文件会占用较大空间，确保有足够磁盘空间
4. **线程设置**: 根据服务器CPU核心数合理设置线程数
5. **数据库下载**: 首次运行会自动下载数据库，可提前下载以节省时间
6. **离线模式**: 服务器无外网时需使用离线模式

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "HMMER3 not found" 错误**

```bash
# 安装HMMER3
conda install -c bioconda hmmer

# 验证安装
hmmsearch --version
```

**Q: "Augustus not found" 错误**

```bash
# 安装Augustus
conda install -c bioconda augustus

# 配置Augustus物种模型
# --augustus-species参数需要已有对应物种模型
```

**Q: 内存不足错误**

```bash
# 减少线程数以降低内存使用
biopytools busco -i genome.fa -l eukaryota_odb12 -t 8 -o ./results

# 或增加系统交换空间
```

**Q: 数据库下载失败**

```bash
# 方法1: 手动下载数据库到指定目录
busco --download-path /path/to/datasets --download eukaryota_odb12

# 方法2: 使用离线模式
biopytools busco -i genome.fa -l eukaryota_odb12 --offline \
    --download-path /path/to/busco_downloads
```

**Q: 批量分析部分样本失败**

```bash
# 查看日志文件确定失败原因
cat busco_output/busco_analysis.log | grep -i error

# 单独重新运行失败的样本
biopytools busco -i failed_sample.fa -l lineage -o ./rerun
```

## 📚 相关资源 | Related Resources

- [BUSCO官方网站](https://busco.ezlab.org/)
- [BUSCO GitHub仓库](https://github.com/ezlab/BUCKy)
- [BUSCO官方文档](https://busco.ezlab.org/busco_userguide.html)
- [BUSCO谱系数据库](https://busco.ezlab.org/portal2/)
- [BUSCO论文 (Nature Methods, 2015)](https://doi.org/10.1038/nmeth.3311)
- [BUSCO v5论文 (Molecular Biology and Evolution, 2021)](https://doi.org/10.1093/molbev/msab053)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

BUSCO软件本身遵循其原始许可证。

---

## 🔬 引用信息 | Citation

如果在学术研究中使用BUSCO工具，请引用原始文献：

```
Simão FA, Waterhouse RM, Ioannidis P, Kriventseva EV, Zdobnov EM.
BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs.
Bioinformatics, 2015, 31(19): 3210-3212.
doi: 10.1093/bioinformatics/btv351

Manni M, Berkeley MR, Glover NM, Waterhouse RM, Ioannidis P, Kriventseva EV, Zdobnov EM.
BUSCO: The essential tool for assessing genome assembly and annotation quality.
Molecular Biology and Evolution, 2021, 38(9): 3882-3895.
doi: 10.1093/molbev/msab053
```

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入文件或目录｜Input file or directory |
| `--lineage, -l` | 必填 | str | BUSCO数据库谱系｜BUSCO database lineage |
| `--output-dir, -o` | `./busco_output` | Path | 输出目录｜Output directory |
| `--mode, -m` | `genome` | genome/geno/transcriptome/tran/proteins/prot | BUSCO分析模式｜BUSCO analysis mode |
| `--threads, -t` | `12` | int | CPU线程数｜Number of CPU threads |
| `--sample-suffix` | `*.fa` |  | 样本名提取后缀｜Sample name extraction suffix pattern |
| `--output-format` | `txt` | txt/csv/xlsx | 输出文件格式｜Output file format |
| `--force, -f` | — |  | 强制重写现有文件｜Force rewriting existing files |
| `--augustus` | — |  | 使用Augustus基因预测器｜Use Augustus gene predictor |
| `--augustus-parameters` | — | str | Augustus额外参数｜Additional Augustus parameters |
| `--augustus-species` | — | str | Augustus物种名｜Augustus species name |
| `--auto-lineage` | — |  | 自动选择谱系｜Automatically select lineage |
| `--auto-lineage-euk` | — |  | 自动选择真核生物谱系｜Automatically select eukaryote lineage |
| `--auto-lineage-prok` | — |  | 自动选择原核生物谱系｜Automatically select prokaryote lineage |
| `--contig-break` | `10` | int | Contig打断的N数量｜Number of Ns for contig break |
| `--datasets-version` | `odb12` |  | 数据集版本｜Dataset version |
| `--download-path` | — | str | 数据集下载路径｜Dataset download path |
| `--evalue, -e` | `0.001` | float | BLAST E值阈值｜BLAST E-value threshold |
| `--limit` | `3` | int | 候选区域限制｜Candidate region limit |
| `--long` | — |  | 启用Augustus长模式优化｜Enable Augustus long mode optimization |
| `--metaeuk` | — |  | 使用Metaeuk基因预测器｜Use Metaeuk gene predictor |
| `--metaeuk-parameters` | — | str | Metaeuk额外参数｜Additional Metaeuk parameters |
| `--metaeuk-rerun-parameters` | — | str | Metaeuk重新运行参数｜Metaeuk rerun parameters |
| `--miniprot` | — |  | 使用Miniprot基因预测器｜Use Miniprot gene predictor |
| `--skip-bbtools` | — |  | 跳过BBTools统计｜Skip BBTools statistics |
| `--offline` | — |  | 离线模式｜Offline mode |
| `--restart, -r` | — |  | 重启未完成的分析｜Restart incomplete analysis |
| `--quiet, -q` | — |  | 静默模式｜Quiet mode |
| `--scaffold-composition` | — |  | 生成scaffold组成文件｜Generate scaffold composition file |
| `--tar` | — |  | 压缩子目录｜Compress subdirectories |
| `--busco-path` | — |  | BUSCO软件路径｜BUSCO software path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入文件或目录路径｜Input file or directory path |
| `-l, --lineage` | 必填 |  | BUSCO数据库/谱系名称｜BUSCO database lineage name |
| `-o, --output-dir` | `./busco_output` |  | 输出目录｜Output directory |
| `-m, --mode` | `genome` | genome/geno/transcriptome/tran/proteins/prot | BUSCO分析模式｜BUSCO analysis mode |
| `-t, --threads` | `12` | int | CPU线程数｜Number of CPU threads |
| `--sample-suffix` | `*.fa` |  | 样本名称提取后缀模式｜Sample name extraction suffix pattern |
| `--output-format` | `txt` | txt/csv/xlsx | 输出文件格式｜Output file format |
| `-f, --force` | — | store_true | 强制重写现有文件｜Force rewriting of existing files |
| `--augustus` | — | store_true | 使用Augustus基因预测器｜Use Augustus gene predictor |
| `--augustus-parameters` | — |  | Augustus附加参数｜Additional Augustus parameters |
| `--augustus-species` | — |  | Augustus物种名称｜Augustus species name |
| `--auto-lineage` | — | store_true | 自动选择谱系｜Automatically select lineage |
| `--auto-lineage-euk` | — | store_true | 自动选择真核生物谱系｜Automatically select eukaryote lineage |
| `--auto-lineage-prok` | — | store_true | 自动选择原核生物谱系｜Automatically select prokaryote lineage |
| `--contig-break` | `10` | int | Contig断点N数量｜Number of Ns for contig break |
| `--datasets-version` | `odb12` |  | 数据集版本｜Dataset version |
| `--download-path` | — |  | 数据集下载路径｜Dataset download path |
| `-e, --evalue` | `0.001` | float | BLAST E值阈值｜BLAST E-value threshold |
| `--limit` | `3` | int | 候选区域数量限制｜Candidate region limit |
| `--long` | — | store_true | 启用Augustus长模式优化｜Enable Augustus long mode optimization |
| `--metaeuk` | — | store_true | 使用Metaeuk基因预测器｜Use Metaeuk gene predictor |
| `--metaeuk-parameters` | — |  | Metaeuk附加参数｜Additional Metaeuk parameters |
| `--metaeuk-rerun-parameters` | — |  | Metaeuk重运行参数｜Metaeuk rerun parameters |
| `--miniprot` | — | store_true | 使用Miniprot基因预测器｜Use Miniprot gene predictor |
| `--skip-bbtools` | — | store_true | 跳过BBTools统计｜Skip BBTools statistics |
| `--offline` | — | store_true | 离线模式｜Offline mode |
| `-r, --restart` | — | store_true | 重启未完成的分析｜Restart incomplete analysis |
| `-q, --quiet` | — | store_true | 静默模式｜Quiet mode |
| `--scaffold-composition` | — | store_true | 生成scaffold组成文件｜Generate scaffold composition file |
| `--tar` | — | store_true | 压缩子目录｜Compress subdirectories |
| `--busco-path` | — |  | BUSCO软件路径｜BUSCO software path |

<!-- END PARAMS:auto -->
