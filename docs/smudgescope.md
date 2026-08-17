# smudgescope 基因组评估 | GenomeScope2 + Smudgeplot Genome Evaluation

一句话理解：**只用测序得到的 FASTQ 数据（不需要参考基因组），就能估出基因组大小、杂合度、重复序列含量和倍性**，是组装前「先摸清家底」的一站式工具。

## 功能概述 | Overview

- 整合 GenomeScope 2.0 与 Smudgeplot，一条命令完成 k-mer 分析与基因组特征评估
- 估算基因组大小、杂合度、重复序列比例、k-mer 覆盖度、读错误率等关键指标
- Smudgeplot 倍性分析（默认开启，可跳过），自动推断基因组倍性（1-6）
- 支持单个 FASTQ 或整个目录（自动按样本分组、批量处理）
- 断点续传：已完成步骤按输出文件存在性自动跳过（见 FAQ）
- 多样本汇总：一键生成 summary 表（TSV + XLSX）与各样本结果图

## 快速开始 | Quick Start

```bash
biopytools smudgescope -i fastq_dir -o output_dir -l 150
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| k-mer | 把测序序列切成一段段固定长度（如 21 个字母）的小片段；所有小片段统计一遍出现次数，就是「k-mer 谱」 |
| k-mer 覆盖度(kcov) | 平均每个 k-mer 在数据里出现了几次；太高说明测序太深，太低说明数据不够 |
| 基因组大小 | 由 k-mer 谱反推出来的「基因组有多大」，是评估的核心输出之一 |
| 杂合度 | 两条同源染色体之间有多少位置不一样；杂合度高的基因组两条链差异大 |
| 倍性 | 一个细胞里同一套染色体有几份；人/多数动植物是二倍体(2)，很多作物是多倍体 |
| Smudgeplot | 画一张「倍性指纹图」，靠杂合 k-mer 的组合模式推断是几倍体 |

## 输入 | Input

- 单个 FASTQ 文件（`.fastq` / `.fq` / `.fastq.gz` / `.fq.gz`），或一个含 FASTQ 的目录（递归查找）。
- 目录模式下按 `--read1-suffix` 模式把文件归组为样本（默认 `*_1.clean.fq.gz`，自动配对 `_1` / `_2`）。
- 双端数据同一样本的两端文件会被合并进同一次 k-mer 计数。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 告诉工具「数据在哪、结果放哪、读长多少」。读长 `-l` 用于 GenomeScope 建模，务必按真实测序读长填写。

- `-i, --input`：输入 FASTQ 文件或目录。
- `-o, --output-dir`：输出目录。
- `-l, --read-length`：测序读长（默认 150）。

### k-mer 与计数 | K-mer & counting

**通俗理解|In plain words:** 这组决定「怎么数 k-mer」。k 越大越能区分重复序列但需要更深的测序；哈希表大小不够会中途报错；最大覆盖度截断高覆盖的重复 k-mer 防止它们淹没信号。**绝大多数项目用默认值即可。**

- `-k, --kmer-size`：k-mer 大小（默认 21）。
- `-s, --hash-size`：Jellyfish 哈希表大小（默认 `10G`），数据量大时调大（如 `100G`）。
- `-c, --max-kmer-cov`：最大 k-mer 覆盖度（默认 1000），超过的计入该档。

### 计算资源 | Compute

**通俗理解|In plain words:** 控制用多少线程并行。线程越多越快，但受机器核数限制，给太多反而可能拖慢。

- `-t, --threads`：线程数（CLI 默认 12）。

### 倍性分析 | Ploidy

**通俗理解|In plain words:** 告诉工具「这个物种大概是几倍体」。默认 2 即可，Smudgeplot 会据此上下文推断实际倍性；确定是多倍体物种时再显式指定。

- `--ploidy`：基因组倍性 1-6（默认 2）。
- `--skip-smudgeplot`：跳过 Smudgeplot 倍性分析，只做 GenomeScope（省时间）。

### FastK 与样本识别 | FastK & sample grouping

**通俗理解|In plain words:** 这组只在跑 Smudgeplot 时才涉及。FastK 表是 Smudgeplot 的输入；`--read1-suffix` 决定程序如何从文件名认出「哪个样本、哪一端」。**一般不用动。**

- `--fastk-table`：直接指定一个已生成的 FastK 表文件路径（跳过 FastK 建表，加速重跑）。
- `--fastk-memory`：FastK 内存大小（CLI 默认 `100G`）。
- `--read1-suffix`：Read1 文件后缀模式（默认 `*_1.clean.fq.gz`）。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先用 Jellyfish 数 k-mer，再用 GenomeScope 拟合基因组模型拿指标；接着用 FastK 建表、Smudgeplot 画倍性指纹图。

1. Jellyfish `count`：统计 k-mer，产出 `.jf` 计数库
2. Jellyfish `histo`：生成 k-mer 频数直方图 `.histo`
3. GenomeScope 2.0：拟合模型，产出 `model.txt` / `summary.txt` / 图，并提取 k-mer 覆盖度(kcov)
4. FastK：为 Smudgeplot 生成 k-mer 表（`fastk_table.ktab` + `.hist`）
5. Smudgeplot `hetmers`：用阈值 `int(kcov×0.5)` 提取杂合 k-mer 对，产出 `.kmerpairs.smu`
6. Smudgeplot `plot`：画倍性指纹图 `{sample}_smudgeplot.png`

## 输出 | Output

```text
output_dir/
├── {sample}/                          # 每个样本一个子目录
│   ├── 00_pipeline_info/software_versions.yml   # 软件版本信息
│   ├── 01_jellyfish/
│   │   ├── {sample}.jellyfish.jf                # k-mer 计数库
│   │   └── {sample}.jellyfish.histo             # k-mer 频数直方图
│   ├── 02_genomescope/
│   │   ├── model.txt                 # 模型参数（含 kcov、杂合度 r 等）
│   │   ├── summary.txt               # 关键指标汇总
│   │   ├── linear_plot.png           # 拟合曲线图（线性坐标）
│   │   └── log_plot.png              # 拟合曲线图（对数坐标）
│   ├── 03_smudgeplot/                # 仅未跳过 Smudgeplot 时
│   │   ├── {sample}.kmerpairs.smu    # 杂合 k-mer 对数据
│   │   └── {sample}_smudgeplot.png   # 倍性指纹图
│   ├── fastk/                        # 条件目录，仅跑 Smudgeplot 时
│   │   ├── fastk_table.ktab
│   │   └── fastk_table.hist
│   ├── 99_logs/genomescope_pipeline.log
│   └── tmp/                          # 临时文件（运行结束清理）
└── all/                              # 全局汇总
    ├── {sample}_genomescope.png      # 汇总拷贝的各样本 GenomeScope 图
    ├── {sample}_smudgeplot.png       # 汇总拷贝的各样本 Smudgeplot 图
    ├── summary.tsv                   # 多样本指标汇总表
    └── summary.xlsx                  # 同上的 Excel 版（装 openpyxl 时）
```

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看三个地方：summary 表拿数字，linear_plot 看拟合好坏，smudgeplot 看倍性。

### 1. summary.txt / model.txt（数字指标）

- `genome haploid length`：单倍体基因组大小估计值。
- `genome repeat length` / `genome unique length`：重复序列与唯一序列长度。
- `heterozygous`（杂合度）：越高说明两条链差异越大；异常高提示可能混入了污染或样品不纯。
- `read error rate`：读错误率估计。
- `model fit`：模型拟合度（越接近 100% 越好，过低说明模型和数据不符）。
- `kmercov`（k-mer 覆盖度）：用于下游判断测序深度是否足够。

### 2. linear_plot.png / log_plot.png（拟合曲线）

看实测 k-mer 分布（离散点）与模型拟合曲线（实线）是否贴合；主峰清晰、拟合线贴得好即数据质量佳。若主峰矮胖、拖尾重，可能是测序深度不足或数据有污染。

### 3. smudgeplot.png（倍性指纹图）

图上的「团块」位置对应不同倍性组合；通常看哪个团块最亮来判断倍性。纯合二倍体可能因杂合度太低而画不出图（属正常，见 FAQ）。

### 4. summary.tsv / summary.xlsx（多样本汇总）

一表对比所有样本的基因组大小、杂合度、错误率、拟合度、k-mer 覆盖度，便于批量项目横向比较。

## 参数选择建议 | Parameter Guidance

- **`-k`**：基因组 1Gb 以内用默认 21；更大基因组可适当增大（如 27、31），但要保证测序深度足够。
- **`-s`（哈希表）**：Jellyfish 报「hash 不足」相关错误时再调大（如 `100G`），平时默认 `10G`。
- **`-c`（最大覆盖度）**：默认 1000 覆盖绝大多数场景；高重复基因组或极高测序深度时可调大，避免高覆盖 k-mer 被过度截断。
- **`--ploidy`**：默认 2 即可，Smudgeplot 会推断；确定是多倍体（如四倍体小麦）时设为对应值。
- **`--skip-smudgeplot`**：只想快速拿基因组大小、不关心倍性时加上，能省掉 FastK 和 Smudgeplot 两步。
- **`--fastk-table`**：重跑多个参数（如只换 GenomeScope 参数）时，指定已生成的 FastK 表可跳过 FastK 建表。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入FASTQ文件或目录｜Input FASTQ file or directory |
| `--output-dir, -o` | 必填 | Path | 输出目录｜Output directory |
| `--read-length, -l` | `150` | int | 测序读长｜Read length |
| `--kmer-size, -k` | `21` | int | K-mer大小｜K-mer size |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--hash-size, -s` | `10G` |  | Jellyfish哈希表大小｜Jellyfish hash size |
| `--max-kmer-cov, -c` | `1000` | int | 最大k-mer覆盖度｜Max k-mer coverage |
| `--skip-smudgeplot` | `False` |  | 跳过Smudgeplot倍性分析｜Skip Smudgeplot ploidy analysis |
| `--ploidy` | `2` | int | 基因组倍性 1-6 (默认: 2，由Smudgeplot自动推断)｜Genome ploidy level 1-6 (default: 2, auto-inferred by Smudgeplot) |
| `--fastk-table` | `` |  | FastK表文件路径｜FastK table file path |
| `--fastk-memory` | `100G` |  | FastK内存大小｜FastK memory size |
| `--read1-suffix` | `*_1.clean.fq.gz` |  | Read1文件后缀模式｜Read1 file suffix pattern |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTQ文件或目录｜Input FASTQ file or directory |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-l, --read-length` | 必填 | int | 测序读长｜Read length |
| `-k, --kmer-size` | `21` | int | K-mer大小｜K-mer size |
| `-t, --threads` | `64` | int | 线程数｜Number of threads |
| `-s, --hash-size` | `10G` |  | Jellyfish哈希表大小｜Jellyfish hash size |
| `-c, --max-kmer-cov` | `1000` | int | 最大k-mer覆盖度｜Max k-mer coverage |
| `--skip-smudgeplot` | — | store_true | 跳过Smudgeplot倍性分析｜Skip Smudgeplot ploidy analysis |
| `--ploidy` | `2` | int | 基因组倍性 1-6 (默认: 2，由Smudgeplot自动推断)｜Genome ploidy level 1-6 (default: 2, auto-inferred by Smudgeplot) |
| `--genomescope-env` | `genomescope_v.2.0.1` |  | GenomeScope conda环境名称 (默认: genomescope_v.2.0.1)｜GenomeScope conda env name (default: genomescope_v.2.0.1) |
| `--fastk-table` | `` |  | FastK表文件路径｜FastK table file path |
| `--fastk-memory` | `16G` |  | FastK内存大小｜FastK memory size |
| `--read1-suffix` | `*_1.clean.fq.gz` |  | Read1文件后缀模式｜Read1 file suffix pattern |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Jellyfish（k-mer 计数）
- GenomeScope 2.0（`genomescope2`，依赖 R/Rscript）
- FastK（Smudgeplot 预处理，仅未跳过时）
- Smudgeplot（含 C 后端 Logex，仅未跳过时）
- 以上软件经 conda 环境自动检测并调用；GenomeScope 默认锁定环境名 `genomescope_v.2.0.1`（模块直调可用 `--genomescope-env` 指定）

## 常见问题 | FAQ

**Q1：支持断点续传吗？**
支持。每一步都按「输出文件是否存在」判断是否跳过（`.jf` / `.histo` / `model.txt` / `fastk_table.ktab` / `.kmerpairs.smu` / 最终图）。换参数重跑前，需先删除对应的旧产物，否则会复用旧结果。

**Q2：Smudgeplot 没生成图，报「未找到 .smu 文件」正常吗？**
正常。纯合二倍体基因组的杂合 k-mer 对太少，Smudgeplot 检测不到信号就不画图，不影响 GenomeScope 的基因组大小等结果。

**Q3：`-t` 线程数和 `--fastk-memory` 的默认值为什么和直调模块不一样？**
`biopytools smudgescope`（CLI）默认线程 12、FastK 内存 `100G`；直接 `python -m biopytools.smudgescope`（模块直调）默认线程 64、FastK 内存 `16G`。两者底层参数存在历史差异，以你实际调用方式为准。

**Q4：输出路径里有中文会报错吗？**
GenomeScope/R 对中文路径敏感，工具已尽量用相对路径规避；建议输入输出路径避免中文字符。

**Q5：kcov 提取失败会怎样？**
会告警并使用默认值 50.0 继续跑 Smudgeplot，结果可能不理想；此时优先检查 `02_genomescope/model.txt` 或 `summary.txt` 是否正常生成。

**Q6：FastK 报内存相关错误？**
FastK 的 `-M` 参数只接受纯数字（不带 G），工具已自动剥掉单位；若还报内存不足，调大 `--fastk-memory`。
