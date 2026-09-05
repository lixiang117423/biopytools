# ENA 数据下载 | ENA Data Download

一句话理解：**给定 ENA 编号（项目号、run 号均可，也可以是一个每行一个编号的 ID 文件批量下载），自动抓取元数据、找出所有 FASTQ 文件的下载链接，并生成可直接执行的下载脚本**，省去手动去网页上一个个点下载的麻烦。

## 功能概述 | Overview

- 输入一个 ENA accession（项目或 run 编号），或一个每行一个编号的 ID 文件批量处理
- 自动识别编号类型（项目级/运行级/样本级）并写入日志
- 从元数据里提取 FASTQ 下载链接，生成可执行的下载脚本（wget 或 Aspera）
- 支持两种执行方式：save（只生成脚本，自己再跑）和 run（当场直接下载）
- 支持仅下载元数据模式，不碰 FASTQ；断点续传——已有元数据的编号自动跳过
- 元数据可存 tsv / csv / xlsx 三种格式，每个编号附带汇总报告，批量时额外生成总览

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

### accession 编号 或 ID 文件

`-a` 接受两种输入：

1. **单个编号**：项目级 `PRJNA661210`、`PRJEB...`、`ERP...`，run 级 `SRR12345678`、`ERR...`、`DRR...` 均可，工具把它传给 ENA 的 filereport API 查询
2. **ID 文件**：一个已存在的文件路径，**每行一个编号**，空行和 `#` 开头的注释行自动跳过，文件内重复编号自动去重：

```text
# 我的项目清单
PRJNA661210
SRR9969612
ERR204944
```

无需其他本地输入文件，唯一「输入」就是编号串或这份清单。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** `-a` 后面跟一个编号，或一份每行一个编号的清单文件。它决定「去 ENA 查哪些数据」。编号写错或 ENA 里查不到时程序会明确报错退出，不会静默产出空文件。

### 输出设置 | Output settings

**通俗理解|In plain words:** `-o` 指定结果写到哪个目录（不写就用当前目录）；`-d`（`--create-dir`）是「给我建个专用文件夹」开关——开了它且不写 `-o`，单编号会自动建 `<accession>.ena.download` 目录、多编号（ID 文件）建 `ena_batch_download` 目录装所有结果，避免互相覆盖。

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
输入(编号 或 ID文件)
  |
  v
对每个编号循环 | for each accession:
  1. 请求 ENA filereport API 下载元数据 -> <accession>_meta.<format>
     (已有含数据行的元数据则跳过, 断点续传)
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
  4. 生成该编号的汇总报告 -> <accession>_download_summary.txt
  (某个编号失败只跳过它, 不影响其余编号继续)
批量(ID文件)时最后生成总览 -> batch_overview.txt
```

## 输出 | Output

```text
输出目录/
├── PRJNA661210_meta.tsv                          # 元数据表(格式由 -f 决定), 每个编号一份
├── download_PRJNA661210_fastq_by_wget.sh         # 下载脚本(ftp 协议), 每个编号一份
├── download_PRJNA661210_fastq_by_aspera.sh       # 下载脚本(aspera 协议,二选一)
├── PRJNA661210_download_summary.txt              # 该编号的汇总报告
├── batch_overview.txt                            # 批量(ID文件)时的总览: 各编号状态
└── ena_download.log                              # 运行日志
```

- `*_meta.tsv / _meta.csv / _meta.xlsx`：从 ENA API 拉回的元数据表，含 run、样本、物种、fastq_ftp 下载链接等列；文件名以编号开头，批量时互不覆盖
- `download_*_fastq_by_wget.sh`：save 模式生成的下载脚本，内含逐条 `wget -c` 命令（`-c` 支持断点续传），执行 `bash download_xxx.sh` 即可开始下载
- `download_*_fastq_by_aspera.sh`：aspera 协议对应的脚本（含 ascp 命令）
- `<accession>_download_summary.txt`：每个编号一份，记录协议、方式、发现的 FASTQ 文件数、下一步操作提示
- `batch_overview.txt`：仅 ID 文件批量时生成，逐行列出每个编号「成功|ok / 无元数据|no metadata」及成功计数

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 先看每个编号的 `<accession>_download_summary.txt` 里「发现的 FASTQ 文件数量」，确认抓到的东西符合预期；再跑下载脚本，下载完核对文件数是否与元数据行数一致。批量时最后扫一眼 `batch_overview.txt`，哪行不是「成功|ok」就单独补跑哪个编号。

- **FASTQ 文件数**：汇总报告中的「发现的FASTQ文件数量」。为 0 通常说明该 accession 没有公开的 FASTQ（可能是项目元数据问题或选错了编号）
- **下载脚本**：save 模式下报告里会提示「执行以下命令开始下载」。脚本带 `set -e`，某条下载失败会中止，可单独重跑
- **好坏判据**：元数据表有数据行、链接列非空、下载脚本里的命令数等于预期文件数，即为正常；程序明确报「编号在ENA无数据」说明该编号在 ENA 查不到或未公开，先核对编号拼写

## 参数选择建议 | Parameter Guidance

- **只想看看项目有什么数据**：加 `--metadata-only`，先下元数据再决定
- **日常下载**：默认参数（ftp + save），生成脚本后后台跑
- **要下几百 GB 大项目**：换 `--protocol aspera --aspera-key 私钥路径`（先把私钥 `chmod 600`），速度快很多
- **一次处理多个项目**：把编号写进一个清单文件传给 `-a`（每行一个，支持 `#` 注释），再加 `-o` 或 `--create-dir` 指定独立输出目录，避免文件互相覆盖
- **中断后接着下**：直接用同一个 ID 文件重跑，已有完整元数据的编号自动跳过，只补缺的

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--accession, -a` | 必填 |  | ENA编号或ID文件路径(每行一个, 支持#注释)｜ENA accession or ID file (one per line, # comments) |
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
| `--accession, -a` | 必填 |  | ENA编号或ID文件路径(每行一个编号, 支持#注释)｜ENA accession or ID file path (one per line, # comments supported) |
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

**Q6：批量清单（ID 文件）怎么写？**
普通文本文件，每行一个编号，空行和 `#` 开头的注释行会被跳过，重复编号自动去重：
```text
# 我的项目清单
PRJNA661210
SRR9969612
```

**Q7：报「编号在ENA无数据|No data found in ENA」？**
传给 ENA 的编号在库里查不到任何 run——常见原因：编号拼写错误；传了样本级编号（ERS/SAMN 开头）但该样本没有公开 run；数据尚未公开。先核对编号，或到 ENA 网页端搜索确认。批量场景下单个编号失败不影响其余编号继续处理。
