# fastq 完整性检查 | check-reads

给你的测序数据(fastq)做一次"体检"——检查压缩文件有没有损坏、有没有 0 字节的空文件、双端测序的 R1/R2 是不是成对齐全,一次性把所有有问题的文件列出来。

## 功能概述 | Overview { #overview }

- **gz 完整性**:逐文件流式解压校验,发现被截断/损坏的压缩文件(下载中断、转换出错的常见后果)
- **0 字节检查**:找出 fastp 等工具"跑过但没产出"的空文件(如质控参数把全部 reads 滤掉)
- **配对完整性**:自动识别双端命名(`_1/_2`、`_R1/_R2`、`_1.clean` 等),标记缺 R1 或缺 R2 的样本
- **并行校验**:多线程并发解压,88 线程跑 100GB 级别数据也只要几分钟
- **多目录支持**:逗号分隔一次检查多个目录(如 `2nd/clean,3rd/clean`),递归扫描子目录
- 输出人类可读报告 + 日志,有问题时退出码为 1(便于作业调度判断成败)

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools check-reads -i 2nd/clean/ -o check_out/ -t 88
```

检查 `2nd/clean/` 下所有 fastq,报告写到 `check_out/check_reads_report.txt`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗解释<br>In plain words |
|------|------|
| gz 完整性 | 压缩文件就像拉链,中途断掉(下载中断、拷贝出错)就拉不开;校验就是"全拉开试试",拉不开就是损坏 |
| 0 字节文件 | 文件存在但内容为空——质控工具跑了但一条 reads 都没留下 |
| R1/R2 配对 | 双端测序同一条片段的两头,缺了一头下游分析(如比对、组装)会出错 |

## 输入 | Input { #input }

一个或多个 fastq 目录(逗号分隔),递归扫描子目录:

- 支持后缀:`.fq`、`.fastq`,可带 `.gz`
- 命名识别(配对检查用):`SRR123_1.fq.gz`+`_2.fq.gz`、`A_R1.fq`+`A_R2.fq`、`SRR123_1.clean.fq.gz`+`_2.clean.fq.gz` 等
- 单端文件(`C.fq` 等无法配对的)合法,报告里单独统计
- 非 fastq 文件忽略并记录在报告尾部

## 参数说明 | Parameters { #parameters }

**通俗理解|In plain words:** 三个参数都直白——`-i` 告诉它查哪里(可多个目录),`-t` 开多少线程(数据量大就开满,如 88),`-o` 报告放哪。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | fastq 目录(逗号分隔多个,递归扫描)｜FASTQ dir(s), comma-separated, recursive |
| `-o, --output-dir` | `./check_reads_output` |  | 输出目录(默认./check_reads_output)｜Output directory |
| `-t, --threads` | `12` | int | 并行线程数(默认12)｜Parallel threads (default 12) |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | fastq 目录(逗号分隔多个,递归扫描)｜FASTQ dir(s), comma-separated, recursive |
| `-o, --output-dir` | `./check_reads_output` |  | 输出目录(默认./check_reads_output)｜Output directory |
| `-t, --threads` | `12` | int | 并行线程数(默认12)｜Parallel threads (default 12) |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

<!-- END PARAMS:auto -->

## 输出 | Output { #output }

```
output/
├── check_reads_report.txt   # 检查报告(逐文件状态 + 配对 + 汇总)
└── 99_logs/
    └── check_reads.log      # 运行日志
```

报告格式示例:

```
=== check_reads 报告|report | 2026-08-27 09:22:19 ===
-- 文件状态|file status --
OK       /path/SRR1_1.clean.fq.gz
CORRUPT  /path/SRR2_1.clean.fq.gz      ← 损坏
EMPTY    /path/SRR3_1.clean.fq.gz      ← 空文件
-- 配对检查|pairing --
缺R2 missing_R2: SRR4                   ← 只有 R1 没有 R2
-- 汇总|summary --
文件总数|files: 58
OK: 55  CORRUPT: 1  EMPTY: 1  PLAIN: 1
配对样本|paired samples: 28
单端|single-end: 1
结果|RESULT: 发现问题(见上方)|issues found
```

## 结果解读 | Interpreting Results { #interpreting-results }

- **`OK`**:gz 完整,可直接用于下游
- **`CORRUPT`**:压缩文件损坏——重新下载/重新转换该样本(常见于下载中断、磁盘故障)
- **`EMPTY`**:空文件——检查上游质控为什么什么都没留下(如 fastp 参数把 reads 全滤掉)
- **`PLAIN`**:明文 fastq(非 gz),不做压缩校验,只查空
- **`缺R2/缺R1`**:双端样本缺一侧——补齐或确认该样本本来就是单端
- **退出码 0** = 全部通过,**1** = 有问题(作业系统可直接据此判失败)

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

| 场景<br>Scenario | 推荐参数<br>Recommended |
|------|------|
| 常规检查 | `-t 12` 即可 |
| 100GB+ 大目录 | `-t 88` 开满核,几分钟出结果 |
| 多批数据一次查 | `-i dir1/,dir2/` 逗号分隔 |
| 提交作业自动判成败 | 直接看退出码(0 通过/1 有问题) |

**通俗理解|In plain words:** 先小线程快速试一下,数据量大再开满线程;报告文件就是给下游/同事看的体检单。

## 依赖 | Dependencies { #dependencies }

| 依赖<br>Dependency | 说明<br>Note |
|------|------|
| Python 3.10+ | 内置 gzip/concurrent.futures,无第三方依赖 |

## 常见问题 | FAQ { #faq }

- **为什么建议在计算节点跑?** gz 校验要读全文件,100GB+ 数据在登录节点会拖慢别人;开满线程在计算节点几分钟完成
- **`PLAIN` 是什么意思?** 明文 fastq(没有 .gz 后缀),没法做压缩完整性校验,只检查了是否为空
- **`CORRUPT` 的文件还能用吗?** 不能,截断的 gz 解压出的数据不完整,必须重新获取
- **检查会影响原文件吗?** 不会,只读不改
