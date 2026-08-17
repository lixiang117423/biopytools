# get-link-from-CNCB - 从 CNCB 批量获取测序数据下载链接 | Batch Download Links from CNCB

一句话理解：**给一个「项目 + 样本」编号清单，自动到 CNCB（GSA/NGDC）的 FTP/HTTPS 上把这些测序数据的下载链接一个个找出来，整理成可执行的下载脚本**。找不到的 ID 还会自动去 ENA、NCBI 兜底再查一遍。

## 功能概述 | Overview

- 输入一份两列表格（项目编号 + Run ID），自动按前缀识别数据来源：SRR→SRA、ERR→ERA、DRR→DRA、CRR→GSA
- 三种查找通道层层兜底：CNCB FTP 镜像 → ENA Portal API → NCBI SDL 数据定位，保证尽量找到链接
- CRR（GSA 原生 Run ID）走 HTTPS 通道，不依赖 FTP；GSA 项目列填 CRA 编号直达，填 PRJCA 等 BioProject 编号也能自动反查
- 自动生成三类产物：链接清单（每行一个 URL）、失败记录、wget -c 下载脚本（含断点续传参数）
- 输出一份人类可读的总结报告，统计成功率、各回退通道命中数、缓存命中率

## 快速开始 | Quick Start

```bash
biopytools get-link-from-CNCB -i projects.txt
```

最小输入：一个两列 Tab 分隔的文本文件 projects.txt（第一列项目编号，第二列 Run ID）。输出默认写到与输入同目录的 <输入名>_links.txt 等文件。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 测序数据 | 一个样本测出来的原始读段文件（如 FASTQ/SRA），通常要下载到本地分析 |
| Run ID | 某次测序运行的编号，像快递单号，本工具就是拿它去「查快递」 |
| 项目编号 | 一组 Run 的集合编号（如 CRA/PRJCA），像「整批订单」的编号 |
| 前缀 | Run ID 前三个字母，用来判断它属于哪家数据库（SRR=SRA，ERR=ENA，DRR=DRA，CRR=GSA） |
| FTP/HTTPS 镜像 | 同一批数据的多份「仓库副本」，本工具从多个仓库找同一份文件 |
| 回退（fallback） | 一个地方查不到，自动换下一个地方查，像「这家店没货换下一家」 |

## 输入 | Input

一个 UTF-8 编码的文本文件，两列用 Tab 分隔（制表符，不是空格），每行一个「项目编号 + Run ID」组合：

```text
# 注释行以 # 开头，会被跳过
CRA010060	CRR123456
PRJCA001234	CRR234567
PRJNA1014406	SRR28526560
```

格式要点：

- 第一列是项目编号，第二列是 Run ID，必须严格两列 Tab 分隔，否则该行报错
- 允许空行和 # 开头的注释行
- 同一个项目下的 Run ID 可以写多行，程序会自动分组、去重、排序
- CRR（GSA 原生）的 Run ID，项目列建议直接填它的 CRA 编号（如 CRA010060）可直达；填 PRJCA 等 BioProject 编号时程序会通过 NGDC 搜索页反查，但依赖网络可达
- Run ID 前缀决定查找通道：SRR/ERR/DRR 走 CNCB FTP 镜像（找不到再回退 ENA/NCBI）；CRR 走 GSA HTTPS 通道

## 参数说明 | Parameters

### 必需与输出参数 | Required and output

**通俗理解|In plain words:** -i 是那份两列清单；-o / --failed / --download-script 决定「找到的链接、没找到的 ID、下载脚本」分别写到哪。不指定输出文件名时，会按输入文件基名自动生成 <输入名>_links.txt、<输入名>_failed.txt、download.sh，一般不用动。

### FTP 连接参数 | FTP connection

**通俗理解|In plain words:** 连 CNCB FTP 服务器的地址、超时和重试次数。默认服务器一般不用改；网络差、频繁超时时可以适当调大 --ftp-timeout 和 --retry-attempts，但调太大只会让失败等更久。

### 回退开关 | Fallback switches

**通俗理解|In plain words:** 控制「CNCB 找不到时要不要去 ENA / NCBI 再查」。默认三层兜底全开（找到链接概率最大）。--no-ena-fallback 会关掉整条回退链（纯 CNCB 模式）；--no-ncbi-fallback 只关 NCBI 一级、保留 ENA。只有当网络只能访问 CNCB、或想严格限定数据来源时才需要关。

### 日志与脚本开关 | Logging and script switches

**通俗理解|In plain words:** -v 打印更详细的调试日志；--log-file 把日志同时写进文件；--no-download-script 不生成下载脚本；--no-executable 生成的脚本不加可执行权限。这些都属于「锦上添花」，正常使用基本不用管。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先读清单 → 连 FTP → 逐个 Run 找文件，找不到的按「ENA → NCBI」顺序兜底，最后统一写结果和脚本。

```text
输入 projects.txt（项目编号 + Run ID 两列）
    │
    ▼
按前缀分类：SRR/ERR/DRR → FTP 通道；CRR → GSA HTTPS 通道
    │
    ├─ CRR：GSA 浏览页精确链接 → 失败回退 autoindex 目录探测
    │
    ├─ SRR/ERR/DRR：CNCB FTP 镜像逐文件模板匹配
    │
    ▼
CNCB 未命中的 INSDC ID → ENA Portal API 回退（批量 + 逐个两阶段）
    │
    ▼
ENA 也未命中的 ID → NCBI SDL 数据定位回退（并发查询）
    │
    ▼
写输出：链接清单 + 失败记录 + wget 下载脚本 + 总结报告
```

## 输出 | Output

```text
<输入文件所在目录>/
├── <输入名>_links.txt            # 找到的下载链接（每行一个 URL，已排序）
├── <输入名>_failed.txt           # 未找到的 ID（两列：项目编号 \t Run ID）
├── download.sh                   # 自动生成的下载脚本（wget -c 断点续传）
└── CNCB_download_report.txt      # 总结报告（统计信息，人类可读）
```

- <输入名>_links.txt：所有成功找到的下载 URL，直接可用于下载。每个 Run 可能对应多个文件（如 R1/R2、.sra 等），都会单独列一行。
- <输入名>_failed.txt：最终仍未找到链接的 ID，格式与输入一致（两列 Tab），可直接作为下次重试的输入。注意：重跑时该文件会被无条件覆盖（包括清空）。
- download.sh：按成功链接生成的 Bash 脚本，每条 wget -c '<url>'；-c 表示支持断点续传，中断后重跑不会重头下。用法：bash download.sh。
- CNCB_download_report.txt：总结报告，含项目数、Run 数、成功/失败数、成功率、总链接数。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看两件事——_links.txt 里有没有你要下载的文件，_failed.txt 里还剩多少没找到。

- 成功率：运行结束日志会打印 成功率|Success Rate。100% 说明全部找到；低于 100% 时看 _failed.txt 是哪些 ID 没找到。
- _failed.txt 非空的原因：ID 拼写错误、项目编号填错（尤其 CRR 应填 CRA 而非 BioProject）、数据未公开/已撤下、或网络不通。
- 回退命中：日志里的 ENA回退命中数、NCBI SDL回退命中数 说明有多少 ID 是 CNCB 找不到、靠兜底通道找到的。若大量命中回退通道，说明 CNCB 镜像不含这些数据，属正常现象。
- 缓存命中率：程序内部会缓存「父目录路径」避免重复扫描 FTP，命中率高说明批量 ID 前缀相近、扫描更省时，仅作参考。

## 参数选择建议 | Parameter Guidance

- 默认参数即可覆盖绝大多数场景，最常用的只有 -i。
- 只想拿链接清单、不生成下载脚本：加 --no-download-script。
- 网络受限、只能访问 CNCB：加 --no-ena-fallback（会同时关掉 NCBI 回退）。
- 想保留 ENA 回退但不要 NCBI 兜底：加 --no-ncbi-fallback。
- 网络不稳：调大 --retry-attempts（如 5）和 --ftp-timeout（如 120）。
- GSA 数据（CRR 前缀）：项目列务必填 CRA 编号，避免每次都要反查、减少失败概率。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入文件｜Input file path (ProjectID and Run ID; GSA项目(CRR前缀)项目列用CRA编号｜use CRA accession for GSA projects with CRR runs) |
| `--output, -o` | — | Path | 输出链接文件｜Output links file path |
| `--failed` | — | Path | 失败记录文件｜Failed records file path |
| `--download-script` | `download.sh` | Path | 下载脚本名｜Download script filename |
| `--ftp-host` | `download2.cncb.ac.cn` |  | FTP服务器｜FTP server host |
| `--ftp-timeout` | `60` | int | FTP连接超时｜FTP connection timeout in seconds |
| `--retry-attempts` | `3` | int | FTP重试次数｜FTP connection retry attempts |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose output mode |
| `--log-file` | — | Path | 日志文件路径｜Log file path |
| `--no-download-script` | — |  | 不生成下载脚本｜Don't generate download script |
| `--no-executable` | — |  | 不设置可执行权限｜Don't make script executable |
| `--no-ena-fallback` | — |  | 关闭ENA/NCBI回退查询(纯CNCB模式)｜Disable ENA/NCBI fallback queries (pure CNCB mode) |
| `--no-ncbi-fallback` | — |  | 仅关闭NCBI SDL回退(保留ENA)｜Disable only the NCBI SDL fallback (keep ENA) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-o, --output` | — |  | 输出文件路径｜Output file path (default: [input]_links.txt) |
| `-f, --failed` | — |  | 失败记录文件路径｜Failed records file path (default: [input]_failed.txt) |
| `--download-script` | — |  | 下载脚本文件名｜Download script filename (default: download.sh) |
| `--ftp-host` | `download2.cncb.ac.cn` |  | FTP服务器地址｜FTP server host (default: download2.cncb.ac.cn) |
| `--ftp-timeout` | `60` | int | FTP连接超时时间｜FTP connection timeout in seconds (default: 60) |
| `--retry-attempts` | `3` | int | FTP连接重试次数｜FTP connection retry attempts (default: 3) |
| `-v, --verbose` | — | store_true | 详细输出模式｜Verbose output mode |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `--no-download-script` | — | store_true | 不生成下载脚本｜Don't generate download script |
| `--no-executable` | — | store_true | 不设置脚本执行权限｜Don't make script executable |
| `--no-ena-fallback` | — | store_true | 关闭ENA/NCBI回退查询(纯CNCB模式)｜Disable ENA/NCBI fallback queries (pure CNCB mode) |
| `--no-ncbi-fallback` | — | store_true | 仅关闭NCBI SDL回退(保留ENA)｜Disable only the NCBI SDL fallback (keep ENA) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3 标准库（ftplib 等）即可完成纯 CNCB FTP 模式
- requests 与 urllib3（用于 ENA/NCBI/GSA 的 HTTP 通道，未安装时这些回退通道优雅降级为空，不报错中断）
- 无需 conda 环境，无需额外生信软件；需要能访问外网（CNCB、EBI、NCBI）

## 常见问题 | FAQ

**Q1：重跑之后 _failed.txt 里的旧失败记录怎么没了？**
这是刻意行为：失败记录每次无条件重写（包括空文件），避免上次的陈旧失败记录残留误导。重跑前请先保存好上次的 _failed.txt，或直接用 _links.txt 判断本次结果。

**Q2：CRR 开头的 Run 为什么提示「项目列请用 CRA 编号」？**
CRR 是 GSA 原生 Run ID，不存放在 CNCB 的 FTP INSDC 镜像里，也无法通过 ENA/NCBI 兜底（INSDC 侧没有 GSA 数据）。项目列填 CRA 编号（如 CRA010060）可直达；填 PRJCA 等 BioProject 编号时需要联网反查，查不到就记为失败。

**Q3：FTP 连不上会直接失败吗？**
不会。SRR/ERR/DRR 的 ID 会先记为失败，再交给 ENA 回退查询；CRR 走 HTTPS 通道、本就不依赖 FTP。所以即使 FTP 不通，只要 ENA 回退开着，仍可能找到大部分链接。

**Q4：为什么同一个 Run 会输出多个链接？**
一个 Run 通常包含多个文件（双端测序的 R1/R2、单文件 .sra 等），程序会把每个真实存在的文件都列出来，属正常现象。下载脚本会逐个 wget。

**Q5：有没有断点续传？**
本工具本身是「查链接」不是「下数据」，不存在中间计算结果可续传；但生成的下载脚本用 wget -c 保证下载环节支持断点续传。
