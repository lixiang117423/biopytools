# ENA 数据下载 | ENA Data Download

一句话理解：**给定一个 ENA 项目编号（accession），自动抓取它的元数据、找出里面所有 FASTQ 文件的下载链接，并生成一份可直接执行的下载脚本**，省去手动去网页上一个个点下载的麻烦。

## 功能概述 | Overview

- 输入一个 ENA accession（项目或 run 编号），自动从 ENA API 拉取元数据表
- 从元数据里提取 FASTQ 下载链接，生成可执行的下载脚本（wget 或 Aspera）
- 支持两种执行方式：save（只生成脚本，自己再跑）和 run（当场直接下载）
- 支持仅下载元数据模式，不碰 FASTQ
- 元数据可存 tsv / csv / xlsx 三种格式，并附带汇总报告

## 快速开始 | Quick Start

```bash
biopytools ena-downloader -a PRJNA661210
```

默认把元数据和下载脚本写到当前目录，元数据格式 tsv，协议 ftp，生成脚本后自行执行。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>Plain meaning |
|------|----------|
| accession | ENA/NCBI 给项目或数据分配的唯一编号，像「项目工号」，用它就能查到整个项目的数据清单 |
| 元数据 metadata | 描述数据本身的信息表（样本名、物种、测序平台、文件大小、下载地址等），像一份带下载链接的目录 |
| FASTQ | 测序仪产出读段的原始格式，是绝大多数分析的第一手输入 |
| wget | 一个命令行下载工具，支持断点续传，适合普通网络 |
| Aspera / ascp | 一个高速传输工具，速度远快于 wget，但需要私钥和专门配置 |
| FTP | 一种文件传输协议，本工具默认用它（配合 wget） |

## 输入 | Input

### accession 编号

ENA accession，例如项目级 `PRJNA661210`、`PRJEB...`，或 run 级 `SRR12345678`。工具把它传给 ENA 的 filereport API 查询。

无需本地输入文件，唯一「输入」就是这一串编号。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 只要一个 accession 编号就能跑。它决定「去 ENA 查哪个项目的数据」。

### 输出设置 | Output settings

**通俗理解|In plain words:** `-o` 指定结果写到哪个目录（不写就用当前目录）；`-d`（`--create-dir`）是「给我建个专用文件夹」开关——开了它且不写 `-o`，会自动建一个叫 `<accession>.ena.download` 的目录装所有结果，适合一次下多个项目时互不混淆。

### 元数据格式 | Metadata format

**通俗理解|In plain words:** 决定元数据表存成什么格式。tsv 最省事（纯文本、好 grep）；csv 适合 Excel 直接打开；xlsx 需要 pandas + openpyxl，没装会自动退回 tsv。一般用默认 tsv 即可。

### 下载协议 | Download protocol

**通俗理解|In plain words:** 选「用普通网络还是高速通道」。ftp（默认）用 wget，稳、通用、支持断点续传；aspera 用 ascp 快得多，但必须配 `--aspera-key` 私钥（且要求私钥权限为 600），配置麻烦，大项目或跨国传输才值得上。

### 执行方式 | Execution method

**通俗理解|In plain words:** save（默认）只把下载命令写进一个 .sh 脚本，你自己 `bash 脚本名` 去跑——好处是可先检查脚本、可后台挂机；run 是当场一条条下载。推荐默认 save。

### 特殊模式 | Special modes

**通俗理解|In plain words:** `--metadata-only` 只看元数据、不下载 FASTQ，适合先摸清项目里有哪些样本再决定下不下。`--fields` 可自定义要哪些元数据字段（默认已带全套常用字段），`--max-retries` 控制 API 请求失败重试次数，一般不用动。

## 分析流程 | Pipeline

```text
输入 accession
  |
  v
1. 请求 ENA filereport API 下载元数据 -> <accession>.meta.<format>
  |
  v
2. 从元数据的 fastq_ftp / fastq_aspera 列提取下载链接
  |
  v
3. 按方法处理:
   save -> 生成 download_<accession>_fastq_by_wget.sh / _by_aspera.sh
   run  -> 当场逐条下载
  |
  v
4. 生成汇总报告 -> download_summary.txt
```

## 输出 | Output

```text
输出目录/
├── PRJNA661210.meta.tsv                          # 元数据表(格式由 -f 决定)
├── download_PRJNA661210_fastq_by_wget.sh         # 下载脚本(ftp 协议)
├── download_PRJNA661210_fastq_by_aspera.sh       # 下载脚本(aspera 协议,二选一)
├── download_summary.txt                          # 汇总报告
└── ena_download.log                              # 运行日志
```

- `*.meta.tsv / .csv / .xlsx`：从 ENA API 拉回的元数据表，含 run、样本、物种、fastq_ftp 下载链接等列
- `download_*_fastq_by_wget.sh`：save 模式生成的下载脚本，内含逐条 `wget -c` 命令（`-c` 支持断点续传），执行 `bash download_xxx.sh` 即可开始下载
- `download_*_fastq_by_aspera.sh`：aspera 协议对应的脚本（含 ascp 命令）
- `download_summary.txt`：项目编号、协议、方式、发现的 FASTQ 文件数、下一步操作提示的总览

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 先看 `download_summary.txt` 里「发现的 FASTQ 文件数量」，确认抓到的东西符合预期；再跑下载脚本，下载完核对文件数是否与元数据行数一致。

- **FASTQ 文件数**：汇总报告中的「发现的FASTQ文件数量」。为 0 通常说明该 accession 没有公开的 FASTQ（可能是项目元数据问题或选错了编号）
- **下载脚本**：save 模式下报告里会提示「执行以下命令开始下载」。脚本带 `set -e`，某条下载失败会中止，可单独重跑
- **好坏判据**：元数据表有行、链接列非空、下载脚本里的命令数等于预期文件数，即为正常

## 参数选择建议 | Parameter Guidance

- **只想看看项目有什么数据**：加 `--metadata-only`，先下元数据再决定
- **日常下载**：默认参数（ftp + save），生成脚本后后台跑
- **要下几百 GB 大项目**：换 `--protocol aspera --aspera-key 私钥路径`（先把私钥 `chmod 600`），速度快很多
- **一次处理多个项目**：加 `--create-dir`，每个项目一个独立文件夹，避免文件互相覆盖

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--accession, -a` | 必填 |  | ENA项目编号｜ENA accession number |
| `--output-dir, -o` | — | Path | 输出目录｜Output directory |
| `--create-dir, -d` | — |  | 创建专门输出目录｜Create dedicated output directory |
| `--metadata-format, -f` | `tsv` | tsv/csv/xlsx | 元数据文件格式｜Metadata file format |
| `--protocol, -p` | `ftp` | ftp/aspera | 下载协议类型｜Download protocol type |
| `--aspera-key, -k` | — |  | Aspera私钥路径｜Path to aspera private key |
| `--method, -m` | `save` | save/run | 执行模式｜Execution mode |
| `--metadata-only, -M` | — |  | 仅下载元数据｜Only download metadata |
| `--fields, -F` | — |  | 自定义元数据字段｜Custom metadata fields |
| `--max-retries, -r` | `3` | int | API请求最大重试次数｜Maximum API request retries |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--accession, -a` | 必填 |  | ENA项目编号｜ENA accession number |
| `--output-dir, -o` | — |  | 输出目录｜Output directory |
| `--create-dir, -d` | — | store_true | 创建专门的输出目录｜Create dedicated output directory |
| `--metadata-format, -f` | `tsv` | tsv/csv/xlsx | 元数据文件格式｜Metadata file format |
| `--protocol, -p` | `ftp` | ftp/aspera | 下载协议类型｜Download protocol type |
| `--aspera-key, -k` | — |  | Aspera私钥路径｜Path to aspera private key |
| `--method, -m` | `save` | save/run | 执行模式｜Execution mode |
| `--metadata-only, -M` | — | store_true | 仅下载元数据，不处理FASTQ文件｜Only download metadata, do not process FASTQ files |
| `--fields, -F` | — |  | 自定义元数据字段｜Custom metadata fields |
| `--max-retries, -r` | `3` | int | API请求最大重试次数｜Maximum API request retries |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 库 requests（调用 ENA API）
- wget（ftp 协议下载，脚本内用 `wget -c`）
- ascp / Aspera（仅 aspera 协议，需自备私钥）
- pandas + openpyxl（仅 `-f xlsx` 时，未安装会自动退回 tsv）

无固定 conda 环境名；wget / ascp 直接通过 PATH 查找。

## 常见问题 | FAQ

**Q1：下载中断了，能续传吗？**
wget 脚本里每条命令都带 `-c`，重跑脚本会自动跳过已下载的部分继续下；aspera 的 ascp 不支持续传。save 模式最省心：脚本没跑完，重新 `bash` 一次即可。

**Q2：aspera 报「密钥文件权限不安全」？**
程序会检查私钥权限必须是 600。执行 `chmod 600 私钥路径` 后再跑即可。

**Q3：xlsx 格式保存不了？**
xlsx 需要 pandas 和 openpyxl。没装时程序会打印警告并自动退回 tsv 格式，装好依赖再重跑才能拿到 xlsx。

**Q4：为什么 csv 文件里字段错位？**
csv 转换是「把制表符换成逗号」的朴素做法，若字段内容本身含逗号会错位。需要可靠表格建议用默认 tsv 或 xlsx。

**Q5：save 和 run 怎么选？**
推荐 save：先生成脚本、人工检查一遍链接和命令，再 `bash` 执行，可后台挂机、可断点续传；run 适合文件少、想一键完成的小项目。
