# 核苷酸多样性(π)计算 | Nucleotide Diversity (π) Calculation

一句话理解：**用 vcftools 计算每个群体内部的核苷酸多样性 π，输出全基因组和滑窗两种粒度的汇总表**，用来比较不同群体的遗传多样性高低。

## 功能概述 | Overview

- 按群体计算 π（群体内核苷酸多样性），每个群体一份结果
- 全基因组模式（默认）：输出每条染色体 + 全基因组整体的 π
- 滑窗模式（`-w`）：输出逐窗口 π，看多样性沿染色体的分布
- 全基因组模式下还会同步跑一份 100kb 滑窗结果（写到 `03_windowed/`）
- 自动建 VCF 的 tabix 索引和基因组的 `.fai` 索引；建 `.fai` 失败（如目录只读）时降级为直接解析 FASTA，结果不受影响
- 断点续传：每个群体已算好的结果自动跳过

## 快速开始 | Quick Start

```bash
biopytools pi -i variants.vcf.gz -p populations.txt -g reference.fasta -o pi_output
```

最小输入：一个 bgzip 压缩的 VCF + 群体文件 + 参考基因组（`.fai` 索引缺失时自动补建）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| π（核苷酸多样性） | 群体里随机挑两个个体，平均每多少个碱基就有一个不一样；越大=多样性越高 |
| 群体 | 一组同来源的样本，程序对每个群体单独算一个 π |
| 滑窗 | 把染色体切成固定大小的「格子」分别算 π，看多样性在哪些区域偏高/偏低 |
| .fai 索引 | 记录每条染色体长度的索引文件；程序会在缺失时用 `samtools faidx` 自动补建 |
| tabix 索引 | 让程序能按位置快速随机读取压缩 VCF 的索引（`.tbi`） |

## 输入 | Input

### VCF 文件

须为 bgzip 压缩（`.vcf.gz`），样本名须与群体文件一致；缺失 `.tbi` 索引时程序会自动补建。

### 群体文件

每行两列「样本ID 群体名」，分隔符可自动识别（TAB / 逗号 / 空格均可）：

```text
sample1    wild
sample2    wild
sample3    cultivated
```

### 参考基因组

FASTA（`.fa` / `.fasta`）或直接给 `.fai`；传 FASTA 时会自动找同目录的 `.fai`，找不到则用 `samtools faidx` 自动补建；若建索引失败（如目录只读），会降级为直接解析 FASTA 拿染色体长度，结果不受影响。

## 分析流程 | Pipeline

```text
输入 VCF + 群体文件 + 参考基因组
    │
    ▼
解析群体、补建 .fai/tabix 索引、读染色体长度
    │
    ▼
对每个群体运行 vcftools：
  ├─ 全基因组模式: --site-pi → 每染色体 π + 整体 π
  └─ 滑窗模式(-w): --window-pi → 逐窗口 π
    │
    ▼
(全基因组模式) 同步跑 100kb 滑窗 → 03_windowed/
    │
    ▼
合并为 pi_merged.tsv + 写软件版本信息
```

## 输出 | Output

```text
pi_output/
├── pi_merged.tsv                    # 汇总结果(核心)
├── 00_pipeline_info/
│   └── software_versions.yml        # 软件版本与参数
├── 01_vcftools/                     # vcftools 原始输出(每群体)
│   ├── {pop}_sites.pi                # 全基因组模式: 位点级 π
│   └── {pop}_windowed.pi             # 滑窗模式: 窗口级 π
├── 03_windowed/                     # 仅全基因组模式: 100kb 滑窗结果
│   ├── 01_vcftools/{pop}_windowed.pi
│   └── pi_merged.tsv
└── 99_logs/
    └── *.log                        # 运行日志
```

## 结果解读 | Interpreting Results

- **`pi_merged.tsv`（核心表）**：
  - 全基因组模式列：`population / chromosome / pi / n_sites`，其中 `chromosome=genome_wide` 那行是该群体整体的 π
  - 滑窗模式列：`population / chromosome / window_start / window_end / pi / n_sites`
- **π 越大 = 遗传多样性越高**；群体间比较时直接比 `genome_wide` 行的 π 值即可
- **滑窗模式看趋势**：π 明显偏低的窗口可能是选择/扫荡留下的信号，明显偏高的窗口可能是平衡选择或结构变异
- **好坏判据**：`n_sites` 太小的群体，π 估计不可靠（依据不足）；比较群体时尽量保证各群体样本量相当

## 参数选择建议 | Parameter Guidance

**通俗理解|In plain words:** 默认跑全基因组模式即可；只有想看「多样性沿染色体分布」时才需要加 `-w`。

- **`-w / --window-step`（滑窗）**：不设 `-w` 就是全基因组模式（推荐先跑这个拿整体 π）；要看分布再加 `-w 100000`（100kb 窗口）
- **`--maf / --max-missing`（质控）**：默认 0.0 / 1.0（不过滤）一般不用动；数据质量差时可设 `--maf 0.05 --max-missing 0.2`
- **`--keep-intermediate`**：默认清理 `_keep.txt` 中间文件；调试时加它保留
- **`--vcftools-path`**：vcftools 不在默认 conda 环境（pop）时，用这个指定路径

## 依赖 | Dependencies

- vcftools（默认 conda 环境 pop 的 `~/miniforge3/envs/pop/bin/vcftools`，可用 `VCFTOOLS_PATH` 环境变量或 `--vcftools-path` 覆盖）
- samtools（默认 conda 环境 align 的 `~/miniforge3/envs/align/bin/samtools`，可用 `SAMTOOLS_PATH` 覆盖；仅用于自动补建 `.fai` 索引）

## 常见问题 | FAQ

**Q1：换参数（如 `--maf`）重跑，结果没变？**
断点续传按每群体的 `_sites.pi` / `_windowed.pi` 是否存在判断。换质控参数后需删除 `01_vcftools/`（和 `03_windowed/`）里的旧结果，否则会复用。

**Q2：参考基因组目录只读，建不了 `.fai` 索引怎么办？**
程序会自动尝试用 `samtools faidx` 补建索引；若失败（目录只读/无权限），会记录 WARNING 并降级为直接解析 FASTA 拿染色体长度，流程照常继续、结果不受影响。

**Q3：VCF 需要什么格式？**
必须是 bgzip 压缩的 `.vcf.gz`；普通 gzip 或未压缩 VCF 无法用 tabix 随机读取。缺失 `.tbi` 索引时程序会尝试自动补建。

**Q4：`03_windowed/` 是干什么的？**
只在全基因组模式下出现：程序顺便用 100kb 窗口（`--default-window-size/--default-window-step` 可调）跑了一份滑窗结果，方便你直接看分布，无需再单独跑一次 `-w`。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | VCF文件路径（bgzip压缩+tabix索引）｜VCF file path (bgzip-compressed + tabix-indexed) |
| `-p, --pop-file` | 必填 |  | 群体文件路径（样本ID 群体名）｜Population file path (sample_id population_name) |
| `-g, --genome` | 必填 |  | 参考基因组fasta文件路径（无.fai索引时自动创建）｜Reference genome fasta path (.fai index auto-created if missing) |
| `-o, --output-dir` | `./pi_output` | Path | 输出目录｜Output directory |
| `-w, --window-size` | — | int | 窗口大小bp（不设置则全基因组计算）｜Window size in bp (omit for genome-wide) |
| `--window-step` | — | int | 窗口步长bp（默认等于窗口大小，无重叠）｜Window step in bp (default=window_size, no overlap) |
| `--default-window-size` | `100000` | int | 自动滑窗默认窗口大小bp（全基因组模式时同步运行）｜Default windowed window size in bp |
| `--default-window-step` | `100000` | int | 自动滑窗默认步长bp｜Default windowed step in bp |
| `--maf` | `0.0` | float | 最小等位基因频率｜Minor allele frequency |
| `--max-missing` | `1.0` | float | 最大缺失率｜Maximum missing rate |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--keep-intermediate` | — |  | 保留中间文件｜Keep intermediate files |
| `--vcftools-path` | — |  | vcftools路径｜vcftools path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --vcf-file` | 必填 |  | VCF文件路径（bgzip压缩+tabix索引）｜VCF file path (bgzip-compressed + tabix-indexed) |
| `-p, --pop-file` | 必填 |  | 群体文件路径（样本ID 群体名）｜Population file path (sample_id population_name) |
| `-g, --genome` | 必填 |  | 参考基因组fasta文件路径（无.fai索引时自动创建）｜Reference genome fasta path (.fai index auto-created if missing) |
| `-o, --output-dir` | `./pi_output` |  | 输出目录｜Output directory |
| `-w, --window-size` | — | int | 窗口大小bp（不设置则全基因组计算）｜Window size in bp (omit for genome-wide) |
| `--window-step` | — | int | 窗口步长bp（默认等于窗口大小）｜Window step in bp (default=window_size) |
| `--default-window-size` | `100000` | int | 自动滑窗默认窗口大小bp｜Default windowed window size in bp (default: 100000) |
| `--default-window-step` | `100000` | int | 自动滑窗默认步长bp｜Default windowed step in bp (default: 100000) |
| `--maf` | `0.0` | float | 最小等位基因频率｜Minor allele frequency (default: 0.0) |
| `--max-missing` | `1.0` | float | 最大缺失率｜Maximum missing rate (default: 1.0) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `--vcftools-path` | — |  | vcftools路径｜vcftools path |

<!-- END PARAMS:auto -->
