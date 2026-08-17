# Assembly 转 AGP 格式 | Assembly to AGP Converter

**将 ALLHiC 等流程输出的 .assembly 文件转换为标准 AGP 格式，并生成染色体列表 | Convert .assembly output into standard AGP format and a chromosome list**

## 功能概述 | Overview

`assembly2agp` 是一个轻量级的组装后处理工具，用于将 ALLHiC / JBAT 等挂载流程产生的 `.assembly` 文本文件（描述 contig 在 scaffold/染色体上的排列与方向）转换为通用的 AGP（Assembly-AGP）格式。AGP 是 NCBI、Ensembl 等数据库广泛采用的组装描述标准，便于后续提交和可视化。

该模块同时会根据用户指定的染色体数量，按 scaffold 长度从大到小排序，取前 N 条作为目标染色体，输出对应的 `chr.list` 文件，可直接作为 ALLHiC 等下游工具的输入。Gap 大小默认 100 bp，可自定义。典型使用场景包括：基因组组装挂载后整理、染色体提交前格式转换、与 Juicebox / ALLHiC 流程衔接。

## 快速开始 | Quick Start

```bash
# 基本用法：12 条染色体
biopytools assembly2agp -a corrected_asm.FINAL.assembly -p output_prefix -n 12

# 指定输出目录和 gap 大小
biopytools assembly2agp -a final.assembly -p chr_asm -n 21 \
    -o ./agp_results -g 200 -f
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-a, --assembly` | 输入 `.assembly` 文件路径（ALLHiC/JBAT 风格）|
| `-p, --prefix` | 输出前缀，用于 AGP 和 chr.list 两个文件 |
| `-n, --num-chromosomes` | 目标染色体数量（必须大于 0）|

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `.` | 输出目录 |
| `-g, --gap` | `100` | Scaffold 之间 gap 大小（bp）|
| `-f, --force` | 关 | 强制覆盖已存在的输出文件 |
| `-v, --verbose` | 关 | 详细模式（`-v`: INFO，`-vv`: DEBUG）|
| `--quiet` | 关 | 静默模式（仅输出 ERROR）|
| `--log-file` | 无 | 日志文件路径 |

（运行 `biopytools assembly2agp -h` 查看完整参数列表）

## 输出 | Output

```
./
├── {prefix}.agp          # 标准 AGP 格式文件（9 列）
└── {prefix}.chr.list     # 染色体列表：scaffold名 <TAB> 长度
```

AGP 文件包含 9 列：Chromosome、Start、End、Order、Tag（W=contig / U=gap）、Contig_ID、Contig_start、Contig_end、Orientation（+/-）。chr.list 按长度从大到小给出前 N 条 scaffold 的名称与总长。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--assembly, -a` | 必填 |  | 组装文件｜Assembly file path |
| `--prefix, -p` | 必填 |  | 输出前缀｜Output prefix |
| `--output-dir, -o` | `.` | Path | 输出目录｜Output directory path |
| `--gap, -g` | `100` | int | Scaffold间隙(bp)｜Scaffold gap size in bp |
| `--num-chromosomes, -n` | 必填 | int | 染色体数量｜Number of chromosomes |
| `--force, -f` | — |  | 强制覆盖｜Force overwrite existing files |
| `--verbose, -v` | — |  | 详细模式(-v: INFO, -vv: DEBUG)｜Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — |  | 静默模式(仅ERROR)｜Quiet mode (ERROR only) |
| `--log-file` | — | Path | 日志文件｜Log file path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-a, --assembly` | 必填 |  | Input assembly file path |
| `-p, --prefix` | 必填 |  | Output prefix (will be used for both AGP and chr.list files) |
| `-o, --output-dir` | `.` |  | Output directory path (default: current directory) |
| `-g, --gap` | `100` | int | Gap size between scaffolds in bp (default: 100) |
| `-n, --num-chromosomes` | 必填 | int | Number of chromosomes to include in chr.list file |
| `-v, --verbose` | `0` | count | Verbose mode (-v: INFO, -vv: DEBUG) |
| `--quiet` | — | store_true | Quiet mode (only ERROR) |
| `--log-file` | — |  | Log file path |
| `-f, --force` | — | store_true | Force overwrite existing files |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3.7+
- pandas（DataFrame 处理）

## 引用 | Citation

- AGP 规范：NCBI Assembly-AGP Specification v2.0
- 若源自 ALLHiC 流程：Zhang, L. et al. ALLHiC: scaffolding large ploidy genomes using Hi-C data. *Nature Methods* 16, 1325-1326 (2019).

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
