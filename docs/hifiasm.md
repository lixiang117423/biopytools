# HiFiasm 基因组组装 | HiFiasm Genome Assembly

一句话理解：**用 HiFi reads 组装基因组，并自动完成质量评估和统计**。输入一份 HiFi reads，输出组装 FASTA，外加 BUSCO（完整性）、QUAST（连续性）评估报告和统计表，让「拼得好不好」一目了然。

## 功能概述 | Overview

- 基于 hifiasm 组装，支持 HiFi / ONT / Hi-C 数据混合
- 自动跑 BUSCO（评估基因完整性）和 QUAST（评估组装连续性），均可 `--skip-*` 关闭
- 自动统计 N10–N90、L50、GC、N 含量等，生成报告和对比表
- 输出格式转换：GFA → FASTA、生成 BED 和 samtools 索引
- 注意：`--resume` 标志存在但**尚未实现**，重跑会重新开始（见 FAQ）

## 快速开始 | Quick Start

```bash
biopytools hifiasm -i sample.hifi.fq.gz -o hifiasm_results -p sample
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| HiFi reads | PacBio 高保真长读长，又长又准 |
| contig | 拼出的连续序列片段 |
| N50 | 从长到短排 contig，累加到总长一半时那条 contig 的长度，越大越完整 |
| BUSCO | 用一套「几乎所有物种都该有的保守基因」做体检，看组装里能找回多少，比例越高越完整 |
| QUAST | 汇总组装的各种连续性和统计指标（N50、contig 数、GC 等） |
| 单倍型(haplotype) | 二倍体里来自父/母的两套基因组；hifiasm 可把它们分开输出 hap1/hap2 |
| GFA | 记录组装图（片段间连接关系）的文件格式，需转成 FASTA 才方便下游用 |

## 输入 | Input

### reads 文件

HiFi reads，支持 `.fq/.fastq/.fq.gz/.fastq.gz/.fa/.fasta/.fa.gz/.fasta.gz`；可选 ONT reads（`--ont-reads`）或 Hi-C（`--hi-c-1`/`--hi-c-2`，两端必须同时给）。

## 参数说明 | Parameters

### 必需参数与基本参数 | Required & basic

**通俗理解|In plain words:** `-i` 输入，`-o` 输出目录，`-p` 前缀，`-t` 线程数。

### 组装参数 | Assembly parameters

**通俗理解|In plain words:** `--hg-size` 是「预算」（基因组大小），默认 `auto` 会按输入文件大小自动估算；`-l`/`--purge-max` 控制重复序列的清理力度，`-s` 控制单倍型分离的松紧。**这些默认值（purge 3、purge-max 65、相似度 0.75）对常规二倍体已调好，一般不用动**；`--assembly-type` 用于二倍体/三倍体/多倍体切换。

### 质量评估 | Quality assessment

**通俗理解|In plain words:** 组装完自动做「体检」。BUSCO 查基因完整性，QUAST 查连续性。`--busco-lineage` 默认 `auto`（自动选谱系，植物默认用 brassicales_odb10）；`--reference-genome` 若给了参考基因组，QUAST 还能做比对评估。不想跑就用 `--skip-busco`/`--skip-quast` 关掉（提速）。

### 分析、输出与系统参数 | Analysis, output & system

**通俗理解|In plain words:** `--analyze-haplotypes` 单独分析两套单倍型；`--min-contig-length` 过滤太短的 contig；`--keep-intermediate` 保留中间文件；`--memory` 内存上限（GB），程序会先估算推荐内存并在不够时警告。

## 分析流程 | Pipeline

```text
HiFi reads
    │
    ▼
步骤1 检查依赖 → 步骤2 估算资源 → 步骤3 预处理(格式/磁盘/续传检查)
    │
    ▼
步骤4 hifiasm 组装(assembly/)
    │
    ▼
步骤5 格式转换(GFA→FASTA) → 步骤6 统计(statistics/)
    │
    ▼
步骤7 质量评估(BUSCO/QUAST, quality_assessment/)
    │
    ▼
步骤8 单倍型分析(可选) → 步骤9 生成最终结果(final_results/) → 步骤10 清理临时文件
```

## 输出 | Output

```text
hifiasm_output/
├── assembly/                          # hifiasm 原始产物
│   ├── {prefix}.bp.p_ctg.gfa          # primary contigs(主结果,GFA)
│   ├── {prefix}.bp.a_ctg.gfa          # alternate contigs
│   ├── {prefix}.log                   # hifiasm 日志
│   └── (三倍体另有 .bp.hap1/hap2.p_ctg.gfa)
├── final_results/                     # 整理后的最终结果
│   ├── {prefix}.primary.fasta         # 主组装 FASTA(最常用)
│   ├── {prefix}.alternate.fasta       # 备选组装
│   ├── assemblies/                    # 标准化命名的组装文件
│   ├── quality_assessment/            # busco/ quast/ statistics/ 汇总
│   ├── file_manifest.txt              # 文件清单
│   └── analysis_summary.txt           # 结果摘要
├── quality_assessment/                # busco/ quast/ 原始结果
├── statistics/
│   ├── assembly_statistics_report.txt # 详细统计报告
│   └── assembly_statistics_comparison.csv  # 各组装对比表
├── format_conversion/                 # BED 文件 + .fai 索引
├── logs/                              # hifiasm_时间戳.log
├── tmp/                               # 临时文件(运行结束清理)
└── hifiasm_config.json                # 本次运行的配置快照
```

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 先看 `final_results/{prefix}.primary.fasta` 的 N50 和 BUSCO 完整度——N50 大 + BUSCO 完整度高，就是好组装。

- **N50 / contig 数**（`statistics/assembly_statistics_report.txt`）：N50 越大、contig 越少越完整
- **BUSCO 完整度**（`final_results/quality_assessment/busco/` 的 CSV/TXT）：`C:`（完整）比例 >90% 算好，<80% 提示组装或谱系选择有问题
- **N 含量**：代表组装里的空洞，越低越好（接近 0 最佳）
- **GC 含量**：应与该物种预期一致，异常偏离提示污染
- `{prefix}.primary.fasta` 是「主组装」，下游分析优先用它；`alternate.fasta` 是备选/次要单倍型

## 参数选择建议 | Parameter Guidance

- `--hg-size`：默认 `auto` 自动估算；若知道确切大小，显式给出（如 `1.4g`）可优化内存
- `--busco-lineage`：默认 `auto`（植物默认 brassicales_odb10）；已知物种谱系时显式指定（如 `embryophyta_odb10`）更准
- `--assembly-type`：二倍体保持 `auto`；三倍体/多倍体显式选 `triploid`/`polyploid`
- `--skip-busco --skip-quast`：只想快速拿组装、暂不做评估时用它提速
- `--memory`：程序会按输入文件大小估算推荐内存，不足会警告；大型基因组给 100G 以上

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input-reads, -i` | 必填 |  | HiFi测序数据文件｜Input HiFi sequencing data file |
| `--output-dir, -o` | `./hifiasm_output` | Path | 输出目录｜Output directory |
| `--prefix, -p` | `sample` | str | 输出文件前缀｜Output file prefix |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--hg-size` | `auto` | str | 基因组大小估计(如1.4g, 2.1g)｜Genome size estimation (e.g., 1.4g, 2.1g) |
| `--purge-level, -l` | `3` | int | Purge级别(0-3)｜Purge level (0-3) |
| `--purge-max` | `65` | int | 最大purge覆盖度｜Maximum purge coverage |
| `--similarity-threshold, -s` | `0.75` | float | 相似性阈值｜Similarity threshold |
| `--ont-reads` | — |  | ONT长读长数据文件｜ONT long-read data file |
| `--hi-c-1` | — |  | Hi-C第一端数据文件｜Hi-C first-end data file |
| `--hi-c-2` | — |  | Hi-C第二端数据文件｜Hi-C second-end data file |
| `--extra-hifiasm-args` | `` | str | 额外的HiFiasm参数｜Additional HiFiasm arguments |
| `--skip-busco` | — |  | 跳过BUSCO质量评估｜Skip BUSCO quality assessment |
| `--busco-lineage` | `auto` | str | BUSCO谱系数据集｜BUSCO lineage dataset |
| `--busco-mode` | `genome` | genome/proteins/transcriptome | BUSCO评估模式｜BUSCO assessment mode |
| `--skip-quast` | — |  | 跳过QUAST质量评估｜Skip QUAST quality assessment |
| `--reference-genome` | — |  | 参考基因组文件(用于QUAST)｜Reference genome file (for QUAST) |
| `--analyze-haplotypes` | — |  | 分析单倍型差异｜Analyze haplotype differences |
| `--min-contig-length` | `1000` | int | 最小contig长度过滤｜Minimum contig length filter |
| `--generate-plots` | — |  | 生成可视化图表｜Generate visualization plots |
| `--assembly-type` | `auto` | auto/diploid/triploid/polyploid | 组装类型｜Assembly type |
| `--keep-intermediate` | — |  | 保留中间文件｜Keep intermediate files |
| `--compress-output` | — |  | 压缩输出文件｜Compress output files |
| `--output-formats` | `['both']` | fasta/gfa/both | 输出格式选择｜Output format selection |
| `--memory` | `100` | int | 内存大小(GB)｜Memory size (GB) |
| `--tmp-dir` | — | Path | 临时目录(默认 output_dir/tmp)｜Temporary directory (defaults to output_dir/tmp) |
| `--max-runtime` | `48` | int | 最大运行时间(小时)｜Maximum runtime (hours) |
| `--resume` | — |  | 恢复中断的分析｜Resume interrupted analysis |
| `--hifiasm-path` | `hifiasm` | str | HiFiasm软件路径｜HiFiasm software path |
| `--busco-path` | `busco` | str | BUSCO软件路径｜BUSCO software path |
| `--quast-path` | `quast` | str | QUAST软件路径｜QUAST software path |
| `--python-path` | `python3` | str | Python解释器路径｜Python interpreter path |
| `--samtools-path` | `samtools` | str | Samtools软件路径｜Samtools software path |
| `--busco-db-path` | — | Path | BUSCO数据库路径｜BUSCO database path |
| `--busco-download-path` | — | Path | BUSCO数据集下载路径｜BUSCO dataset download path |
| `--debug` | — |  | 启用调试模式｜Enable debug mode |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose output mode |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |
| `--config-file` | — |  | 配置文件路径｜Configuration file path |
| `--dry-run` | — |  | 试运行模式(不执行实际命令)｜Dry run mode (do not execute actual commands) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-reads` | 必填 |  | 输入HiFi测序数据文件｜Input HiFi sequencing data file |
| `-o, --output-dir` | `./hifiasm_output` |  | 输出目录｜Output directory |
| `-p, --prefix` | `sample` |  | 输出文件前缀｜Output file prefix |
| `-t, --threads` | `32` | int | 线程数｜Number of threads |
| `--hg-size` | `auto` |  | 基因组大小估计(如1.4g, 2.1g)｜Genome size estimation (e.g., 1.4g, 2.1g) |
| `-l, --purge-level` | `3` | int | purge级别(0-3)｜Purge level (0-3) |
| `--purge-max` | `65` | int | 最大purge覆盖度｜Maximum purge coverage |
| `-s, --similarity-threshold` | `0.75` | float | 相似性阈值｜Similarity threshold |
| `--ont-reads` | — |  | ONT长读长数据文件｜ONT long-read data file |
| `--hi-c-1` | — |  | Hi-C第一端数据文件｜Hi-C first-end data file |
| `--hi-c-2` | — |  | Hi-C第二端数据文件｜Hi-C second-end data file |
| `--extra-hifiasm-args` | `` |  | 额外的HiFiasm参数｜Additional HiFiasm arguments |
| `--skip-busco` | — | store_true | 跳过BUSCO质量评估｜Skip BUSCO quality assessment |
| `--busco-lineage` | `auto` |  | BUSCO谱系数据集(如embryophyta_odb10)｜BUSCO lineage dataset |
| `--busco-mode` | `genome` | genome/proteins/transcriptome | BUSCO评估模式｜BUSCO assessment mode |
| `--skip-quast` | — | store_true | 跳过QUAST质量评估｜Skip QUAST quality assessment |
| `--reference-genome` | — |  | 参考基因组文件(用于QUAST)｜Reference genome file (for QUAST) |
| `--analyze-haplotypes` | — | store_true | 分析单倍型差异｜Analyze haplotype differences |
| `--min-contig-length` | `1000` | int | 最小contig长度过滤｜Minimum contig length filter |
| `--generate-plots` | — | store_true | 生成可视化图表｜Generate visualization plots |
| `--assembly-type` | `auto` | auto/diploid/triploid/polyploid | 组装类型｜Assembly type |
| `--keep-intermediate` | — | store_true | 保留中间文件｜Keep intermediate files |
| `--compress-output` | — | store_true | 压缩输出文件｜Compress output files |
| `--output-formats` | `['both']` | fasta/gfa/both | 输出格式选择｜Output format selection |
| `--memory` | `64` | int | 内存大小(GB)｜Memory size (GB) |
| `--tmp-dir` | — |  | 临时目录(默认 output_dir/tmp)｜Temporary directory (defaults to output_dir/tmp) |
| `--max-runtime` | `48` | int | 最大运行时间(小时)｜Maximum runtime (hours) |
| `--resume` | — | store_true | 恢复中断的分析｜Resume interrupted analysis |
| `--hifiasm-path` | `hifiasm` |  | HiFiasm软件路径｜HiFiasm software path |
| `--busco-path` | `busco` |  | BUSCO软件路径｜BUSCO software path |
| `--quast-path` | `quast` |  | QUAST软件路径｜QUAST software path |
| `--python-path` | `python3` |  | Python解释器路径｜Python interpreter path |
| `--samtools-path` | `samtools` |  | Samtools软件路径｜Samtools software path |
| `--busco-db-path` | — |  | BUSCO数据库路径｜BUSCO database path |
| `--busco-download-path` | — |  | BUSCO数据集下载路径｜BUSCO dataset download path |
| `--debug` | — | store_true | 启用调试模式｜Enable debug mode |
| `--verbose, -v` | `0` | count | 详细输出模式(-v, -vv, -vvv)｜Verbose output mode |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |
| `--config-file` | — |  | 配置文件路径｜Configuration file path |
| `--dry-run` | — | store_true | 试运行模式(不执行实际命令)｜Dry run mode (do not execute actual commands) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

以下 conda 软件自动解析功能域环境并经 conda run 调用（默认用裸命令名，可用 `--*-path` 参数覆盖；域环境缺失时回退 PATH 直接调用）：

- hifiasm（组装，asm 域）
- BUSCO（完整性评估，busco 域，可 `--skip-busco` 跳过）
- QUAST（连续性评估，无对应功能域环境，可 `--skip-quast` 跳过）
- samtools（统计/索引，align 域）
- python3（系统自带）
- pandas、matplotlib/seaborn（统计表与绘图，缺省时跳过对应产物）

## 常见问题 | FAQ

**Q1：`--resume` 为什么不起作用？**
当前版本续传功能**尚未实现**（代码中标记「功能开发中」），即使加了 `--resume` 也会从头重新运行。中断后请直接重跑完整流程。

**Q2：BUSCO `auto` 谱系用的是什么？**
`auto` 模式下按内置推荐列表取第一个 `brassicales_odb10`（十字花科）。若你的物种不是十字花科，请用 `--busco-lineage` 显式指定（如 `embryophyta_odb10`）。

**Q3：内存/线程不够会怎样？**
程序会先估算推荐值，若你的配置低于推荐会打印警告（不阻断）。大型基因组建议给足内存（>=100G），否则 hifiasm 可能 OOM。

**Q4：Hi-C 数据两端必须都提供吗？**
必须同时给 `--hi-c-1` 和 `--hi-c-2`，只给一端会报错。

**Q5：`primary.fasta` 和 `alternate.fasta` 区别？**
hifiasm 把「主」单倍型拼进 primary，把「备选」拼进 alternate。二倍体下游通常用 primary.fasta。

