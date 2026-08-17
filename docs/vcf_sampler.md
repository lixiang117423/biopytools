# vcf-sampler - VCF抽样 | VCF Sampling

一句话理解：**从一个很大的 VCF 里按比例随机抽取一部分 SNP 位点（如 25%），得到一个小得多的 VCF**，用于快速测试流程或降低下游计算量；指定随机种子可让抽样结果完全可复现。

## 功能概述 | Overview { #overview }

- 按染色体统计 SNP 数量
- 按比例随机抽取位点（每条染色体至少保留 1 个）
- 优先用 pysam 加速，未安装时自动退回 Python 逐行处理（较慢）
- 自动为输出 VCF 创建索引（`.tbi`）
- 输入若无索引会自动临时创建（tabix），跑完清理

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vcf-sampler -i input.vcf.gz -o output.vcf.gz
```

按默认 25% 比例从 `input.vcf.gz` 随机抽取 SNP 位点，写到 `output.vcf.gz`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| SNP | 基因组上单个碱基的差异位点，本工具抽的就是这些「点」 |
| 抽样比例 | 想保留多少比例的位点，0.25 就是每 4 个留 1 个 |
| 随机种子 | 抽签用的「随机数起点」，种子相同则每次抽到的位点完全一样 |
| 索引 (.tbi) | 给压缩 VCF 做的「目录」，让程序能按区间快速取数 |
| pysam | 一个 Python 库，读 VCF 比逐行解析快很多 |

## 输入 | Input { #input }

输入一个 VCF 文件，**必须是 gzip 压缩格式（`.vcf.gz`）**。普通未压缩的 `.vcf` 会被校验拒绝。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 指定输入和输出 VCF 路径（输出也必须以 `.vcf.gz` 结尾），每次必填。

### 抽样控制 | Sampling control

**通俗理解|In plain words:** `-r` 决定留多少位点（0 到 1 之间，不含 0），`-s` 决定「怎么抽」。想复现上次的抽样结果就固定同一个 `-s`；想每次都不一样就换种子或留默认。**默认 25% + 种子 1288，跑流程测试很合适。**

### 日志 | Logging

**通俗理解|In plain words:** `-v` 看更多过程信息，`--log-file` 把日志写到文件。**一般不用动。**

## 分析流程 | Pipeline { #pipeline }

```text
输入 VCF(.vcf.gz)
    │
    ▼
步骤1: 统计每条染色体的 SNP 数量(pysam 加速 / Python 兜底)
    │
    ▼
步骤2: 按比例随机选每个染色体里要保留的位点索引(每染色体至少1个)
    │
    ▼
步骤3: 写入抽样结果 -> output.vcf.gz，并创建 .tbi 索引
```

## 输出 | Output { #output }

```text
/path/output.vcf.gz        # 抽样后的 VCF(头信息原样保留)
/path/output.vcf.gz.tbi    # 索引
```

输出路径完全由 `-o` 指定，无固定目录结构。

## 结果解读 | Interpreting Results { #interpreting-results }

- 抽样后的 VCF 与输入结构一致（头信息、样本、基因型全部保留），只是数据行变少。
- 日志末尾会打印「总 SNP 数」「抽取的 SNP 数」「实际抽样比例 vs 目标比例」。实际比例会因四舍五入略有偏差，属正常。
- 每条染色体至少保留 1 个位点：即使某条染色体位点极少，也不会被抽空。

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- **`-r/--sample-rate`**：只想快速验证流程可用 `-r 0.1` 或更小；需要代表性较强的小数据集可用 `-r 0.25` 或 `-r 0.5`。
- **`-s/--random-seed`**：做「可复现」实验时固定种子；多次抽样想取不同子集时换种子。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | 必填 | Path | 输出VCF文件路径｜Output VCF file path |
| `--sample-rate, -r` | `0.25` | float | 抽样比例｜Sampling rate |
| `--random-seed, -s` | `1288` | int | 随机种子｜Random seed |
| `--log-file` | — | Path | 日志文件路径｜Log file path |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `-o, --output` | 必填 |  | 输出VCF文件路径｜Output VCF file path |
| `-r, --sample-rate` | `0.25` | float | 抽样比例｜Sampling rate |
| `-s, --random-seed` | `1288` | int | 随机种子｜Random seed |
| `-v, --verbose` | `0` | count | 详细输出模式｜Verbose mode |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `-V, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3
- pysam（可选，用于加速；未安装时自动退回 Python 逐行处理，会变慢）
- tabix / htslib（输入无索引时自动创建索引用；缺失时退回 Python 模式）

## 常见问题 | FAQ { #faq }

**Q1：输入必须是 .vcf.gz 吗？**
是的。普通 `.vcf` 会被校验拒绝，先用 `bgzip input.vcf`（或 `gzip`）压缩成 `.vcf.gz` 再用。

**Q2：没装 pysam 会怎样？**
能正常跑，但会退回逐行解析，大数据集明显变慢。建议安装 pysam 提速。

**Q3：有断点续传吗？**
没有。每次运行都重新统计、重新抽样并覆盖输出文件。

**Q4：能保证抽样完全可复现吗？**
能。固定 `-s/--random-seed` 和 `-r/--sample-rate`，对同一个输入文件，每次抽到的位点完全一致。
