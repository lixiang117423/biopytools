# bam2fastq - BAM 转 FASTQ | BAM to FASTQ Conversion

一句话理解：**把比对/测序生成的 BAM 文件（如 PacBio HiFi）批量还原成 FASTQ 测序读段**，方便重新比对、纠错或喂给下游工具；支持单文件或整个目录，可多文件并行。

## 功能概述 | Overview

- 用 `bam2fastq` 工具把 BAM 批量转回 FASTQ（读段序列），主要用于三代测序（PacBio）BAM 的下游回退或重处理
- 输入支持**单个 BAM 文件**或**包含多个 BAM 的目录**，目录自动扫描所有 `*.bam`
- 输出路径智能识别：以 `.fq.gz` / `.fastq.gz` 结尾视为"输出到指定文件"，否则视为"输出到目录"
- 自动处理 PBI 索引（缺失时用 `pbindex` 生成），并正确解析软链接
- 支持 `-j` 多文件并行 + `-t` 每文件多线程
- 无断点续传（每次都会重新转换，见 FAQ）

## 快速开始 | Quick Start

~~~bash
biopytools bam2fastq -i ./sample.bam -o ./sample.fq.gz -t 32
~~~

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| BAM | 测序读段比对到参考基因组后的"存档格式"，二进制、体积大 |
| FASTQ | 测序读段的原始文本格式，记录每条读段的序列 + 每个碱基的质量分 |
| PBI 索引 | PacBio BAM 的"目录"，工具靠它快速定位读段；缺了就先补一份 |
| 软链接 | 指向真实文件的"快捷方式"，本工具会自动顺着它找到真实文件 |
| 读段(read group) | 同一次测序运行的一组读段；不同读段组会产出多个 FASTQ，工具会合并 |

## 输入 | Input

- 单个 BAM 文件（`.bam`），或一个包含多个 BAM 的目录
- 支持软链接输入（会自动解析到真实文件）
- 典型来源：PacBio HiFi 比对产物

~~~text
# 单文件
biopytools bam2fastq -i sample.bam -o out.fq.gz

# 目录(批量)
biopytools bam2fastq -i ./bam_dir -o ./fastq_dir -j 4 -t 16
~~~

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i` 是输入（BAM 文件或目录），`-o` 是输出（目录，或带 `.fq.gz`/`.fastq.gz` 结尾的单个文件）。两个都必填。给 `-o` 一个以序列扩展名结尾的路径就表示"只转成一个指定名字的文件"，否则按目录输出、每个 BAM 转成同名 FASTQ。

### 并行与线程 | Threads & jobs

**通俗理解|In plain words:** `-t`（默认 12）是**每个** BAM 转换用的线程数，`-j`（默认 1）是**同时处理几个** BAM 文件。文件多、机器核多时调大 `-j` 收益最明显；总线程约等于 `-t × -j`，别超过机器实际核数。

### 工具路径 | Tool path

**通俗理解|In plain words:** `--bam2fastq-path` 指定 `bam2fastq` 可执行文件位置，默认直接找 PATH 里的 `bam2fastq`。装了 conda 环境但 PATH 里没有时才需要显式指定，一般不用动。

## 分析流程 | Pipeline

~~~text
输入 BAM 文件/目录
    │
    ▼
1. 收集 BAM 列表(目录扫描 *.bam / 单文件)
2. 逐文件: 解析软链接 → 检查/生成 PBI 索引(pbindex)
3. bam2fastq 转换(-o 前缀 -j 线程 -c 6 压缩)
4. 文件模式: 重命名/合并 bam2fastq 输出到用户指定文件名
    │
    ▼
输出 FASTQ(.fastq.gz/.fq.gz) + 日志 bam2fastq_conversion.log
~~~

## 输出 | Output

~~~text
输出目录/
├── sample1.fastq.gz              # 每个 BAM 对应的 FASTQ(目录模式, 名字同 BAM)
├── sample2.fastq.gz
└── bam2fastq_conversion.log      # 运行日志
~~~

- **目录模式**（`-o` 不是文件）：每个 BAM 生成同名 FASTQ（`样本名.fastq.gz`）
- **文件模式**（`-o out.fq.gz`）：输出精确命名为 `out.fq.gz`；若 bam2fastq 因多个读段组生成了多个文件，会自动拼接成一个
- 转换失败的 BAM 会记入日志并在结尾的统计里列出

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看日志结尾的成功/失败统计即可判断转换是否正常；理想情况"成功率 100.00%"。

- 日志末尾有 `Total / Success / Failed / Success rate` 四行，`Failed: 0` 即全部成功
- 失败常见原因：PBI 索引生成失败、BAM 损坏、`bam2fastq` 未安装（见 FAQ）
- FASTQ 大小与读段数量应大致正比，异常偏小可能是 BAM 里读段本就不多

## 参数选择建议 | Parameter Guidance

- **单个大 BAM**：`-t 32` 即可，`-j` 保持 1
- **目录里几十个 BAM**：`-j 4 -t 16`，充分利用多核；注意 `-t × -j` 不超过总核数
- **想精确命名输出**：给 `-o` 传 `xxx.fq.gz` 结尾的路径
- **bam2fastq 不在 PATH**：用 `--bam2fastq-path /path/to/bam2fastq`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入BAM文件或文件夹路径｜Input BAM file or directory path |
| `-o, --output-dir` | 必填 | Path | 输出路径(文件或目录,以.fq.gz/.fastq.gz结尾视为文件)｜Output path (file if ends with .fq.gz/.fastq.gz, else directory) |
| `-t, --threads` | `12` | int | 每个BAM文件转换使用的线程数｜Threads per BAM file conversion |
| `-j, --jobs` | `1` | int | 并行处理的BAM文件数量｜Number of parallel BAM file processing |
| `--bam2fastq-path` | `bam2fastq` | str | bam2fastq可执行文件路径｜bam2fastq executable path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入BAM文件或文件夹路径(包含BAM文件)｜Input BAM file or directory path (containing BAM files) |
| `-o, --output-dir` | 必填 |  | 输出路径(文件或目录,自动识别)｜Output path (file or directory, auto-detected) |
| `-t, --threads` | `64` | int | 每个BAM文件转换使用的线程数 (默认: 64)｜Threads per BAM file conversion (default: 64) |
| `-j, --jobs` | `1` | int | 并行处理的BAM文件数量 (默认: 1)｜Number of parallel BAM file processing (default: 1) |
| `--bam2fastq-path` | `bam2fastq` |  | bam2fastq可执行文件路径 (默认: bam2fastq)｜bam2fastq executable path (default: bam2fastq) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- `bam2fastq`（核心转换工具，通过 `--bam2fastq-path` 指定，默认 PATH 查找；自动检测 conda 环境）
- `pbindex`（PBI 索引工具，来自 pbbam/pbmm2 包；缺失索引时自动调用）
- 安装提示（代码内）：`conda install -c bioconda pbmm2 pbbam`

## 常见问题 | FAQ

**Q1：会不会断点续传？**
不会。本工具每次运行都会重新转换（只有 PBI 索引会复用已存在的）。想避免重复转换，重跑前确认目标 FASTQ 还没生成，或手动避开已完成的文件。

**Q2：为什么输出的 FASTQ 名字和我想的不一样？**
`bam2fastq` 自身输出的文件名不稳定（可能是 `.fastq.gz` 也可能是 `.fq.gz`）。文件模式下本工具会自动找到生成的文件并重命名为你指定的名字；目录模式下直接沿用样本名。多个读段组时输出会被合并成一个文件。

**Q3：目录输入 + 文件输出会怎样？**
目录里有多个 BAM 却指定了单文件输出时，会打印警告并自动回退到"目录输出模式"，每个 BAM 各出一个 FASTQ，不会把多个样本混成一个文件。

**Q4：报"未找到 bam2fastq"？**
先确认已安装：`conda install -c bioconda pbmm2 pbbam`，或用 `--bam2fastq-path` 指向实际路径。工具会自动检测它是否在某个 conda 环境里并正确包装调用。

**Q5：软链接的 BAM 能处理吗？**
能。工具会沿软链接解析到真实文件（相对路径的软链接会从多级目录尝试 `realpath` 解析），避免因链接失效而失败。
