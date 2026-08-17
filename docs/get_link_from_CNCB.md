# get-link-from-CNCB - 批量获取测序数据下载链接 | Batch Fetch Sequencing Download Links

一句话：给它一份「项目编号 + 测序数据编号」清单，它帮你把每条数据的真实下载网址找出来，并生成一个可以直接跑的下载脚本——CNCB 上找不到的会自动去 ENA、NCBI 逐级兜底。

|This tool takes a list of "project accession + run accession" pairs, resolves each run to its real download URL, and writes a ready-to-run download script. Runs missed by CNCB fall back to ENA and then NCBI automatically.

## 功能概述 | Overview

- 支持 CNCB(GSA) 原生 `CRR` 与 INSDC 镜像 `SRR`/`ERR`/`DRR` 两类 Run ID
- 三级回退链：CNCB FTP 镜像 → ENA Portal API → NCBI SDL 数据定位（可找到仅存 NCBI S3 的 source BAM 等文件）
- GSA 通道优先刮取 NGDC 浏览页精确链接（覆盖 gsa2/gsa3 布局与 tar.gz 归档），失败自动回退目录探测
- 生成 `download.sh` 下载脚本与失败清单，支持断点续传（`wget -c`）
- 项目列填 `PRJCA` 等 BioProject 编号时，自动按 CRR 反查其所属 `CRA`

## 快速开始 | Quick Start

示例|Examples: biopytools get-link-from-CNCB -i projects.txt

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗解释<br>In plain words |
|--------------|--------------------------|
| Run ID（SRR/ERR/DRR/CRR） | 一次测序的「快递单号」；前缀表示数据存放的数据库（SRR=NCBI、ERR=欧洲、DRR=日本、CRR=CNCB 自营）<br>A sequencing run's tracking number; the prefix tells you which archive holds it |
| Project ID（PRJNA/ERP/DRP/CRA） | 一个课题的「项目编号」，一个项目下通常有多条 Run<br>The study-level accession; one project usually contains many runs |
| FTP 镜像 | CNCB 把 INSDC 公共数据复制了一份放在自己服务器上，国内下载快<br>CNCB's local copy of INSDC public data, fast to download in China |
| 回退（fallback） | 第一个数据库没有时，自动去问第二个、第三个，像打客服电话转接<br>When one archive misses, the query is forwarded to the next, like call-center transfers |
| ENA | 欧洲核酸数据库，INSDC 三兄弟之一，网页 API 好查<br>The European archive, one of the three INSDC siblings, with a friendly API |
| NCBI SDL | NCBI 的「数据定位器」，能查出某条 Run 的文件真实放在哪个云存储上<br>NCBI's data locator; it tells you which cloud bucket a run's files actually live in |
| source BAM | 测序中心上传的原始 BAM（如 PacBio subreads.bam），常只在 NCBI 有<br>The submitter's original BAM (e.g. PacBio subreads), often held only by NCBI |
| wget 下载脚本 | 工具生成的 `download.sh`，每行一条 `wget -c`，断网重跑会接着下<br>The generated script with one `wget -c` per URL; safe to re-run after interruption |

## 输入 | Input

制表符（Tab）分隔的两列表格文件，第一列项目编号，第二列 Run ID：

```
CRA010060	CRR705258
PRJNA123456	SRR29936798
PRJNA123456	SRR12145514
```

要求|Requirements:

- CRR 行的项目列**优先填 CRA 编号**（如 `CRA010060`）；填 `PRJCA` 等 BioProject 编号时工具会自动反查，但依赖 NGDC 可达
- SRR/ERR/DRR 行的项目列填对应 BioProject 编号即可，也可留空由镜像自动定位

## 参数说明 | Parameters

**通俗理解|In plain words:** 这个工具参数很少，绝大多数场景只需要 `-i` 指定输入文件；其余参数只在网络环境特殊或想微调输出时才需要动。

### 回退链控制 | Fallback control

**通俗理解|In plain words:** 「回退」就是 CNCB 查不到时自动去 ENA、NCBI 兜底。默认全开；如果你所在网络访问国外慢、或只想要 CNCB 自己的数据，用这两个开关关掉对应层级。一般不用动。

- `--no-ena-fallback`：整条回退链（ENA + NCBI）全部关闭，纯 CNCB 模式
- `--no-ncbi-fallback`：只关 NCBI SDL 这一级，保留 ENA——适合出国流量受限、`*.s3.amazonaws.com` 云存储地址不可达的机器

### 输出控制 | Output control

**通俗理解|In plain words:** 控制生成哪些文件、脚本有没有执行权限。默认生成下载脚本并加执行权限，一般不用动。

### 网络与连接 | Network & connection

**通俗理解|In plain words:** FTP 地址、超时、重试次数。只有服务器地址变更或网络极差导致频繁断线时才需要调大超时/重试。注意 GSA 浏览页查询用的是独立的短超时（15 秒），不受 `--ftp-timeout` 影响。

## 分析流程 | Pipeline

不同前缀的 Run ID 走不同通道，逐级回退，最多能定位到数据所在的真实位置：

| Run 前缀<br>Run prefix | 主通道<br>Primary channel | 回退链<br>Fallback chain |
|---------|--------|--------|
| `CRR`（GSA 原生） | NGDC 浏览页刮取精确下载链接（覆盖 gsa2/gsa3 布局、tar.gz 归档）；项目列填 PRJCA 时先按 CRR 反查 CRA | 浏览页失败 → autoindex 目录探测（gsa2 → gsa3 → gsa） |
| `SRR`/`ERR`/`DRR`（INSDC） | CNCB FTP INSDC 镜像 | FTP 未找到 → ENA Portal API → 仍无链接 → NCBI SDL 数据定位 |

流程步骤|Steps:

1. 读取并校验输入文件（项目/Run 两列 Tab 分隔）
2. 连接 CNCB FTP 服务器（失败不退出：INSDC ID 交给回退链，CRR 走 HTTPS 不依赖 FTP）
3. 逐条解析：CRR 走 GSA HTTPS 通道，INSDC 前缀走 FTP 镜像
4. 未命中的 INSDC ID 批量回退 ENA；ENA 也无链接的再回退 NCBI SDL
5. 写出链接清单、失败清单、下载脚本与总结报告

## 输出 | Output

| 文件<br>File | 内容<br>Content |
|--------------|----------------|
| `<输入名>_links.txt` | 全部成功解析的下载链接（每行一条，已排序） |
| `<输入名>_failed.txt` | 三级回退后仍失败的「项目 + Run」清单（全成功时也会写空文件，清掉上次残留） |
| `download.sh` | 每行一条 `wget -c '<URL>'` 的下载脚本，可直接 `bash download.sh` 或 `nohup` 后台跑 |
| `<输入名>_report.txt` | 总结报告（成功/失败统计） |

## 结果解读 | Interpreting Results

### 链接清单（`_links.txt`）

- `ftp://download.big.ac.cn/...` 与 `https://download.cncb.ac.cn/...` 是 CNCB 直链，国内网络首选
- `ftp.sra.ebi.ac.uk/...` 来自 ENA 回退；`https://sra-pub-*.s3.amazonaws.com/...` 等对象存储直链来自 NCBI SDL 回退，需要机器能访问对应云存储
- 同一条 Run 出现多个链接是正常的（双端数据 f1/r2、多条文件）

### NCBI SDL 回退的取舍（重要）

SDL 返回的同一文件常有多个云镜像，工具**每个文件只保留一个免费位置**，并主动跳过：

- `payRequired: true` 的镜像（需请求方付费签名，`wget` 无法直接下载）
- `.lite` 无质量精简副本（与完整文件数据重复，下了浪费带宽）
- 只有「待取回」（rehydration）标记、暂时没有链接的文件——会在日志里给出 WARNING

若日志出现「文件无可用的免费下载位置」，说明该文件当前在 NCBI 侧不可直接下载，Run 已按未命中处理或在链接上不完整，可稍后重试或到 SRA Run Browser 手动获取。

### 失败清单（`_failed.txt`）

- CRR 行失败最常见原因是项目列不是 CRA 且 NGDC 反查不可达——把 CRA 编号填进项目列即可
- INSDC 行出现在这里 = CNCB 镜像、ENA、NCBI SDL 三处都没有可用链接

### 下载脚本（`download.sh`）

`wget -c` 支持断点续传；中断后直接重跑脚本即可，已下完的文件会接着补齐。

## 参数选择建议 | Parameter Guidance

| 场景<br>Scenario | 建议<br>Suggestion |
|------|------|
| 常规批量下载（国内超算） | 全部默认即可，三级回退全自动 |
| 出网受限、S3/GCS 对象存储不可达 | 加 `--no-ncbi-fallback`，只保留 CNCB+ENA 两级 |
| 完全离线内网 / 只信 CNCB 镜像 | 加 `--no-ena-fallback`（ENA 与 NCBI 一并关闭） |
| NGDC 页面偶发慢但 FTP 正常 | 不用动；浏览页查询自带 15 秒短超时与连续失败熔断，不会拖垮整批 |
| 同一 Run 挂在多个项目下 | 尽量合并到一行；同一 Run 出现多行时其链接也会重复出现 |

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

- Python 标准库（`ftplib` 等）：纯 CNCB FTP 模式零第三方依赖
- `requests`（可选）：GSA HTTPS 通道、ENA 回退、NCBI SDL 回退需要；未安装时这些通道优雅降级为未命中，不影响 FTP 通道

## 常见问题 | FAQ

**Q1: CRR 报「项目列请用CRA编号」？**
项目列填了非 CRA 编号且 NGDC 反查失败。把该 Run 所属的 CRA 编号（如 `CRA010060`）填进项目列。

**Q2: NCBI 回退给的 S3 链接下载 403/超时？**
这些是对象存储直链，需要网络可达且文件本身免费开放。工具已过滤 `payRequired` 镜像；若所在网络访问 `*.s3.amazonaws.com` 受限，用 `--no-ncbi-fallback` 关闭这一级。

**Q3: 为什么有的 Run 只有 `.sra`/`.bam` 没有 fastq 链接？**
NCBI 只存了原始提交文件（如 PacBio source BAM）时，SDL 返回的就是该文件本身；需要 fastq 请用 `sra-tools` 的 `fasterq-dump` 从 SRA 对象转换。

**Q4: GSA 明明有数据却查不到？**
新项目可能位于 gsa3 布局，工具会自动探测 gsa2 → gsa3 → gsa；若仍失败，确认项目列 CRA 编号正确，并查看日志中浏览页与目录探测的具体错误。

**Q5: 中断后重跑会重复下载吗？**
不会。`download.sh` 用 `wget -c` 断点续传；链接清单与失败清单每次运行都会整体重写。
