# iseq - 公共测序数据下载工具 | Public Sequencing Data Download (iSeq)

一句话理解：**输入一个数据库编号（项目/样本/实验 ID），自动到 GSA/SRA/ENA/DDBJ 等公共数据库把对应的测序数据下载下来，可加 gzip 压缩、可转 FASTQ、可走 Aspera 高速通道**。

## 功能概述 | Overview

- 基于 iSeq 软件封装，支持从 GSA/SRA/ENA/DDBJ 等公共数据库下载测序数据与元数据
- 一个 accession 即可识别项目/样本/实验/Run 各级编号，自动下载对应数据
- 支持 gzip 压缩下载、转 FASTQ、仅下载元数据、合并样本等选项
- 支持 Aspera 高速下载、多线程、多并行连接、下载限速、MD5 校验
- 自动把命令记录到日志，并生成下载汇总报告

## 快速开始 | Quick Start

```bash
biopytools iseq -i PRJNA1014406 -o ./data
```

最小输入：一个数据库编号（如 PRJNA1014406），其余参数用默认值即可。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| accession（编号） | 数据库给某条数据发的唯一编号，像「图书检索号」，输入它就能定位到数据 |
| Run ID（SRR/ERR/DRR/CRR） | 一次测序运行，是最小的下载单元，前缀代表所在数据库 |
| gzip 压缩 | 把 FASTQ 文件压成 .gz 格式，体积更小、传输更快 |
| FASTQ 格式 | 测序原始读段的标准文本格式；SRA 是压缩归档格式，转 FASTQ 便于下游分析 |
| Aspera | 一种高速传输协议，带宽利用率高，适合大文件快速下载 |
| MD5 校验 | 下载后用校验和比对文件是否完整，防止传输损坏 |

## 输入 | Input

输入是一个数据库 accession 编号，支持多种级别，前缀自动识别类型：

```text
项目 Projects：PRJEB / PRJNA / PRJDB / PRJC / GSE
研究 Studies：ERP / DRP / SRP / CRA
样本 BioSamples：SAMD / SAME / SAMN / SAMC
样本 Samples：ERS / DRS / SRS / GSM
实验 Experiments：ERX / DRX / SRX / CRX
运行 Runs：ERR / DRR / SRR / CRR
```

无法识别的前缀也不会报错，程序会尝试继续下载。

## 参数说明 | Parameters

### 必需与输出参数 | Required and output

**通俗理解|In plain words:** -i 是要下载的编号；-o 是数据放哪（默认 ./iseq_output）；--iseq-path 是 iSeq 软件本体位置；-c 是它所在的 conda 环境。后三者都有合理默认值，绝大多数情况只填 -i 就行。

### 下载格式参数 | Download format

**通俗理解|In plain words:** 这三个决定「下载成什么样」。-g 要 gzip 压缩的 FASTQ（省空间）；-q 要解压后的 FASTQ（便于直接分析）；-m 只要元数据（描述信息）不下数据。注意 -g 和 -q 互斥、不能同时开；两个都不加时下载数据库原始格式（通常是 .sra 归档）。

### 性能参数 | Performance

**通俗理解|In plain words:** -t 是处理线程数（用于转换/压缩），-p 是同时开的下载连接数（越大下载越快、但也更容易被服务器限流），-s 限制下载速度（单位 MB/s，避免占满带宽）。一般用默认值即可，网络好想加速可调大 -p。

### 数据源参数 | Source parameters

**通俗理解|In plain words:** -d 选从哪个数据库下载（ena/sra），--protocol 选传输协议（ftp/https），-e 是合并选项（ex/sa/st，特殊需求才用）。这些默认值基本覆盖常规场景，一般不用动。

### 高级选项 | Advanced options

**通俗理解|In plain words:** --use-aspera 走 Aspera 高速通道（需已安装并配置 Aspera）；--skip-md5 跳过完整性校验（只求快、不在乎校验时用）；--quiet 静默模式隐藏进度条。平时不用开。

## 分析流程 | Pipeline

```text
输入 accession
    │
    ▼
校验编号格式 → 构建 iSeq 命令（含各选项）
    │
    ▼
conda run 包装执行 iSeq（自动检测环境）
    │
    ▼
下载数据 / 元数据到输出目录
    │
    ▼
写下载汇总报告 download_summary.txt
```

## 输出 | Output

```text
iseq_output/
├── <数据文件>               # 下载的 FASTQ/SRA 文件（文件名由 iSeq 决定）
├── metadata.txt             # 元数据（启用 --metadata-only 时）
├── download_summary.txt     # 下载汇总报告（选项/性能参数记录）
└── iseq_download.log        # 运行日志（含完整命令）
```

- 数据文件：实际下载的测序数据，格式取决于 -g/-q 等选项。
- metadata.txt：仅元数据模式下输出的元数据信息。
- download_summary.txt：本次下载的基本信息、选项、性能参数汇总。
- iseq_download.log：完整运行日志，记录执行的 iSeq 命令，便于复现与排查。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看输出目录里有没有你期望的数据文件，以及日志里 iSeq 命令是否执行成功。

- 下载成功：输出目录出现数据文件，日志结尾有「iSeq下载完成」。
- 下载失败：日志里有「命令执行失败」和错误信息（如编号错误、网络不通、Aspera 未配置）。
- MD5 校验失败：提示文件损坏，建议重下或改用 --skip-md5 跳过校验继续。
- 建议核对 download_summary.txt 里的选项，确认下载格式、数据库、协议符合预期。

## 参数选择建议 | Parameter Guidance

- 常规下载：只给 -i 编号即可，其余全默认。
- 想要 gzip 压缩省空间：加 -g；想要解压后直接分析：加 -q（两者二选一）。
- 只查有哪些数据、不下数据：加 -m。
- 大文件想提速：加 --use-aspera（需先装好 Aspera）或调大 -p。
- 带宽有限：用 -s 限制速度，避免影响其他任务。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --accession` | 必填 |  | 项目/样本/实验ID｜Project/Sample/Experiment ID (e.g., PRJNA1014406) |
| `-o, --output-dir` | `./iseq_output` | Path | 输出目录｜Output directory |
| `--iseq-path` | `~/miniforge3/envs/misc/bin/iseq` |  | iSeq软件路径｜iSeq software path |
| `-c, --conda-env` | `iseq_v.1.9.8` |  | Conda环境名｜Conda environment name |
| `-m, --metadata-only` | — |  | 仅下载元数据｜Only download metadata |
| `-g, --gzip` | — |  | 下载gzip格式FASTQ｜Download FASTQ in gzip format |
| `-q, --fastq` | — |  | 转换为FASTQ格式｜Convert to FASTQ format |
| `-e, --merge` | — | ex/sa/st | 合并选项｜Merge option (ex/sa/st) |
| `-t, --threads` | `16` | int | 线程数｜Number of threads |
| `-p, --parallel` | `10` | int | 并行连接数｜Number of parallel connections |
| `-s, --speed` | — | int | 下载速度限制(MB/s)｜Download speed limit (MB/s) |
| `-d, --database` | `ena` | ena/sra | 数据库选择｜Database selection |
| `--protocol` | `ftp` | ftp/https | 协议选择｜Protocol selection |
| `-a, --use-aspera` | — |  | 使用Aspera高速下载｜Use Aspera for high-speed download |
| `--skip-md5` | — |  | 跳过MD5校验｜Skip MD5 check |
| `--quiet` | — |  | 静默模式｜Quiet mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --accession` | 必填 |  | 项目/样本/实验ID｜Project/Sample/Experiment ID (e.g., PRJNA1014406) |
| `-o, --output-dir` | `./iseq_output` |  | 输出目录｜Output directory |
| `-a, --iseq-path` | `~/miniforge3/envs/misc/bin/iseq` |  | iSeq软件路径｜iSeq software path |
| `-m, --metadata-only` | — | store_true | 仅下载元数据｜Only download metadata |
| `-g, --gzip` | — | store_true | 下载gzip格式FASTQ｜Download FASTQ in gzip format |
| `-q, --fastq` | — | store_true | 转换为FASTQ格式｜Convert to FASTQ format |
| `-e, --merge` | — | ex/sa/st | 合并选项｜Merge option (ex/sa/st) |
| `-t, --threads` | `16` | int | 线程数｜Number of threads |
| `-p, --parallel` | `10` | int | 并行连接数｜Number of parallel connections |
| `-s, --speed` | — | int | 下载速度限制(MB/s)｜Download speed limit (MB/s) |
| `-d, --database` | `ena` | ena/sra | 数据库选择｜Database selection |
| `--protocol` | `ftp` | ftp/https | 协议选择｜Protocol selection |
| `--use-aspera` | — | store_true | 使用Aspera下载｜Use Aspera for download |
| `--skip-md5` | — | store_true | 跳过MD5校验｜Skip MD5 check |
| `--quiet` | — | store_true | 静默模式｜Quiet mode (suppress progress bars) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- iSeq 软件（默认路径 ~/miniforge3/envs/misc/bin/iseq）
- conda 环境 iseq_v.1.9.8（默认，程序会自动用 conda run 包装）
- Aspera（可选，仅 --use-aspera 时需要）

## 常见问题 | FAQ

**Q1：为什么日志里有 conda run 包装？**
conda 环境里的软件直接调用常因依赖隔离失败，程序会自动检测 iSeq 所在环境并用 conda run -n <env> --no-capture-output 包装执行，属正常行为。

**Q2：-g 和 -q 能同时用吗？**
不能。二者是互斥的（gzip 压缩 vs 转 FASTQ），同时指定会报错。两个都不加则下载数据库原始格式。

**Q3：下载中断能续传吗？**
本模块本身没有断点续传逻辑，重跑会重新下载。大数据建议走 Aspera 或确保网络稳定。

**Q4：提示 iSeq 路径不存在？**
检查 --iseq-path 是否指向真实可执行文件；默认是 ~/miniforge3/envs/misc/bin/iseq，若环境名不同请用 --iseq-path 或 -c 指定。

**Q5：编号无法识别会怎样？**
无法匹配已知前缀时程序只打警告、仍会尝试下载；若 iSeq 报「编号不存在」，请核对编号拼写和所属数据库。
