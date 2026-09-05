# VCF 转 DeepBSA 输入 CSV | VCF to DeepBSA Input CSV { #vcf2deepbsa }

一句话理解：**把 BSA 分析用的 VCF「瘦身」成 DeepBSA 直接认识的 CSV 表**——只保留染色体位置和每个混池的等位基因测序深度(AD)，其余全扔掉，文件更小、DeepBSA 吃得更顺。

> 本模块从 `deepbsa` 的 `vcf2csv` 子命令拆分而来(1.67.0 起)，功能与输出格式完全一致，命令从 `biopytools deepbsa vcf2csv` 变为 `biopytools vcf2deepbsa`。

## 功能概述 | Overview

- 从 VCF 的 FORMAT 字段提取 **AD(Allele Depth,等位基因深度)**，转成 DeepBSA 的无表头 CSV 输入
- 支持 `.vcf` 和 `.vcf.gz`(gzip 压缩)输入
- 输出格式：每行 `CHROM,POS,REF,ALT,池1_REF深度,池1_ALT深度,池2_REF深度,池2_ALT深度,...`
- 坏行不中断：无 AD 字段、AD 缺失(`.`)、数值非法、样本列数不一致的记录**自动跳过并分类计数**(日志可见)，能转多少转多少
- 纯 Python 标准库实现，无第三方依赖，流式处理大 VCF 不占内存

## 快速开始 | Quick Start

```bash
biopytools vcf2deepbsa -i bsa_pools.vcf -o deepbsa_input/
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| VCF | 变异检测的「标准账本」，记录每个位点上测到了什么变异 |
| AD(等位基因深度) | 在某个位点，支持「参考碱基」和「变异碱基」的测序读数各有多少条——像投票的票数 |
| 混池(pool) | 一群个体 DNA 混在一起测序；BSA 通常有高表型池和低表型池两个 |
| BSA(混池分离分析) | 比较两个极端表型混池的「投票差异」，差异最大的区域就是候选基因区 |
| DeepBSA | 一套 BSA 分析工具，同时跑 7 种统计算法找候选区；它吃的是本模块产的 CSV |

## 输入 | Input

- 一个 BSA 混池 VCF 文件(`.vcf` 或 `.vcf.gz`)，由「高表型池 + 低表型池」的变异检测产生
- **硬要求**：样本的 FORMAT 字段必须含 `AD`(如 `GT:AD` 或 `GT:AD:DP`)；没有 AD 的记录会被跳过
- 多个样本(≥2 个混池)都可以，每个样本贡献两列；单样本也能转(但 BSA 至少要两个池才有意义)
- 位点不要求双等位：多等位 ALT(如 `T,G`)只取 AD 的前两个值

```text
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	pool_R	pool_S
chr1	100	.	A	T	50	PASS	.	GT:AD	0/0:10,5	0/1:20,8
chr1	200	.	G	C	50	PASS	.	GT:AD	1/1:0,30	0/0:15,3
```

## 参数说明 | Parameters

**通俗理解|In plain words:** 参数只有三个——`-i` 给输入 VCF，`-o` 给输出目录(不给就落在 `./vcf2deepbsa_output/`)，`--log-level` 调日志详略(一般不用动)。

### 输入输出 | Input & output

**通俗理解|In plain words:** `-i` 是必填的 VCF 路径；`-o` 目录不存在会自动创建，CSV 和日志都放里面。

### 日志 | Logging

**通俗理解|In plain words:** 默认 INFO 级别，统计信息齐全；排查问题时改 `--log-level DEBUG` 才有额外细节，日常无需改。

## 输出 | Output

```text
deepbsa_input/
├── bsa_pools_deepbsa.csv   # DeepBSA 输入 CSV(文件名=输入名_deepbsa.csv)
└── vcf2deepbsa.log         # 运行日志(含转换统计)
```

CSV 内容示例(无表头)：

```text
chr1,100,A,T,10,5,20,8
chr1,200,G,C,0,30,15,3
```

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 打开日志看三个数——总记录数、写入数、跳过数。**写入/总数 越接近 100% 越好**；跳得多时按日志里的分类对症下药(见 FAQ)。

| 日志字段 | 怎么读 |
|------|------|
| 总记录数|Total records | VCF 里的位点总数 |
| 写入记录数|Written records | 成功转成 CSV 的位点数，**核心产物指标** |
| 样本数|Sample count | 每个位点贡献两列 AD 的样本数，应等于混池数 |
| 跳过-无AD字段 | VCF 变异检测时没算 AD(如某些 callers 默认不开)，需重新 call 或换工具 |
| 跳过-AD字段不完整 | AD 只有 1 个值(如 `AD:5`)，属畸形记录 |
| 跳过-AD/POS值非法 | AD 或 POS 是 `.` 或非数字 |
| 跳过-样本列数不一致 | 某行样本数和其他行不一样，通常 VCF 被截断/损坏 |

**好坏判据**：正常 GATK/freebayes 产出的含 AD 的 VCF，跳过率应接近 0；若「无AD字段」占了大多数，说明这个 VCF 根本没有深度信息，转出来也没法跑 BSA，应回上游重新检测。

## 参数选择建议 | Parameter Guidance

| 场景 | 建议 |
|------|------|
| 常规使用 | 全默认，只给 `-i` |
| 输出想放别处 | `-o 指定目录/` |
| 想看哪些行被跳过 | `--log-level DEBUG` 后看日志 |
| 转完直接跑 BSA | `biopytools deepbsa run -i deepbsa_input/bsa_pools_deepbsa.csv -o bsa_results` |

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入 VCF 文件(支持 .gz)｜Input VCF file (.gz supported) |
| `-o, --output-dir` | `./vcf2deepbsa_output` |  | 输出目录(默认./vcf2deepbsa_output)｜Output directory |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-vcf` | 必填 |  | 输入 VCF 文件(支持 .gz)｜Input VCF file (.gz supported) |
| `-o, --output-dir` | `./vcf2deepbsa_output` |  | 输出目录(默认: ./vcf2deepbsa_output)｜Output directory (default: ./vcf2deepbsa_output) |
| `--log-level` | `INFO` |  | 日志级别(默认: INFO)｜Log level (default: INFO) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- 纯 Python 标准库(gzip/csv/argparse)，**无第三方依赖**，基础环境即可运行
- 上游：任何产出带 AD 字段 VCF 的变异检测流程(GATK、freebayes、bcftools 等)
- 下游：DeepBSA(`biopytools deepbsa run/batch` 直接吃本模块的 CSV)

## 常见问题 | FAQ

**Q1：以前用 `biopytools deepbsa vcf2csv`，现在怎么报错了？**
1.67.0 起该子命令拆分为独立模块，改用 `biopytools vcf2deepbsa -i input.vcf -o out_dir/`。注意两点变化：`-o` 从「输出 CSV 文件名」改为「输出目录」；CSV 文件名自动按输入名生成(`xxx.vcf → xxx_deepbsa.csv`)。

**Q2：日志大量「跳过-无AD字段」怎么办？**
说明 VCF 的 FORMAT 里没有 AD。常见原因：变异检测时没开 AD 输出。解决：回上游重新检测(如 bcftools mpileup 默认含 AD)，或用 `bcftools query -f '%CHROM\t%POS[\t%AD]\n'` 先确认。

**Q3：gzip 的 VCF 能直接用吗？**
能，`-i x.vcf.gz` 直接识别，无需先解压。

**Q4：多等位位点(ALT 有多个)会怎样？**
按原 DeepBSA 工具的行为，只取 AD 的前两个值(REF 和第一个 ALT)，REF/ALT 列保持 VCF 原文。建议上游先用 `bcftools view -m2 -M2 -v snps` 过滤成双等位 SNP 再转。

**Q5：输出的 CSV 能给 DeepBSA 的哪些命令用？**
`deepbsa run -i xxx_deepbsa.csv` 和 `deepbsa batch -i xxx_deepbsa.csv` 都行，DeepBSA 对 VCF/CSV 输入同等对待；CSV 更小、解析更快。
