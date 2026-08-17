# busco - 基因组完整性评估 | BUSCO genome completeness assessment

一句话理解：**用一套「几乎每个物种都该有、且只该有一份」的保守基因当尺子，量一量你的基因组/转录组/蛋白集是不是完整、有没有缺胳膊少腿**。
它回答的问题：这个组装/注释「完整度多少、缺了多少、有没有重复组装导致的多拷贝」，是基因组组装和注释质量评估的事实标准。

## 功能概述 | Overview { #overview }

- 基于单拷贝直系同源基因(BUSCO)评估基因组完整性，支持 genome / transcriptome / proteins 三种模式
- 输入支持单文件或目录批量(目录下所有匹配样本逐个评估)
- 结果编译成一张汇总表(支持 txt / csv / xlsx)，并生成总结报告
- 完整度、单拷贝、多拷贝、片段化、缺失五类指标一次给出
- 支持自动选谱系、Augustus / Metaeuk / Miniprot 多种预测器、离线模式
- 断点续传：已完成的样本自动跳过

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools busco -i genome.fa -l eukaryota_odb12
```

最小输入：一个 FASTA(基因组/转录组/蛋白)+ 一个 BUSCO 谱系(如 `eukaryota_odb12`)。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| BUSCO | 一套「基准基因」集合，用来给基因组/注释打完整度分数 |
| 单拷贝直系同源基因 | 进化上保守、几乎所有物种都只有一份的基因，是最公平的「尺子」 |
| 谱系(lineage) | 按物种大类预分好的基准基因集(如 eukaryota 真核、fungi 真菌、metazoa 动物) |
| 完整(Complete) | 基准基因被完整找到，是好事 |
| 单拷贝(Single-copy) | 基准基因只找到一份(最理想) |
| 多拷贝(Duplicated) | 基准基因找到多份，可能是重复组装或近期全基因组复制 |
| 片段化(Fragmented) | 只找到基因的一部分，说明组装有缺口 |
| 缺失(Missing) | 完全没找到，说明这块区域丢了或没注释出来 |
| odb12 | BUSCO 数据集版本号，12 是最新主版本 |

## 输入 | Input { #input }

### 输入文件或目录(-i)

- 单文件：直接给 FASTA
- 目录：默认按 `*.fa` 匹配(用 `--sample-suffix` 改模式)，目录下所有匹配文件逐个评估；匹配不到时会回退尝试常见序列后缀(.faa/.fna/.cds/.pep 等)

### 谱系(-l，必填)

BUSCO 谱系名称，如 `eukaryota_odb12`、`fungi_odb12`、`metazoa_odb12`。选错谱系分数会失真，选与物种最接近的类。

### 分析模式(-m)

- `genome`(默认)：输入是基因组序列，用 Metaeuk 预测再评估
- `transcriptome`：输入是转录组，评估转录组完整性
- `proteins`：输入是蛋白序列，评估蛋白集完整性

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 输入 + 谱系两个必填。谱系是打分基准，**一定要选对大类**，选错整个分数就没意义了。

### 运行参数 | Runtime

**通俗理解|In plain words:** `-o` 输出目录(默认 `./busco_output`)、`-t` 线程数(默认 12，只影响速度)、`-m` 分析模式。`--sample-suffix` 是目录批处理时提取样本名的规则，**一般不用动**。`-f` 强制重跑(默认会跳过已完成样本)。

### 输出格式 | Output format

**通俗理解|In plain words:** `--output-format` 决定汇总表格式，txt(制表符)/csv/xlsx 三选一，**一般用默认 txt 即可**，要交给 Excel 就选 xlsx。

### 谱系选择 | Lineage selection

**通俗理解|In plain words:** 不知道选哪个谱系时，用 `--auto-lineage-euk`(自动选真核谱系)或 `--auto-lineage-prok`(原核)让 BUSCO 自己挑，省心但稍慢。`--datasets-version`(默认 odb12)、`--download-path`(离线数据目录)一般不用动。

### 预测器与高级参数 | Predictors & advanced

**通俗理解|In plain words:** BUSCO 内部用预测器把基因组「翻译」成蛋白再比基准。默认(Metaeuk)对绝大多数情况够用；`--augustus` 可指定用 AUGUSTUS 预测器(需额外配置 species 参数)，`--miniprot` 用 miniprot 处理蛋白输入。`-e`(E 值阈值)、`--limit`(候选数)、`--contig-break` 等是内部阈值，**一般不用动**。

### 运行行为 | Behavior

**通俗理解|In plain words:** `--offline` 离线模式(只用已下载的谱系数据，不联网)、`-r` 重启未完成分析、`-q` 静默、`--scaffold-composition` 生成 scaffold 组成、`--tar` 压缩子目录、`--busco-path` 指定 BUSCO 软件路径。**离线模式默认已开启**，其余按需。

## 分析流程 | Pipeline { #pipeline }

```text
输入: FASTA(单文件或目录批量) + 谱系
    |
    v
步骤1: 获取输入文件(单文件 or 目录批量匹配)
    |
    v
步骤2: 逐样本运行 BUSCO(断点续传, 已完成跳过)
  - 预测器把序列转成蛋白(genome/transcriptome 模式)
  - 与谱系基准基因比对 -> 统计 完整/单拷贝/多拷贝/片段化/缺失
  - 输出 {sample}_busco/ 目录(含 short_summary.*.json)
    |
    v
步骤3: 编译结果 -> 汇总表 + 总结报告 + 版本信息
```

## 输出 | Output { #output }

```text
output_dir/
├── 00_pipeline_info/
│   └── software_versions.yml         # 软件版本
├── {sample}_busco/                   # 每个样本的 BUSCO 原始输出
│   ├── short_summary.specific.<lineage>.{sample}_busco.json   # 结果 JSON
│   ├── short_summary.generic.<lineage>.{sample}_busco.json
│   ├── short_summary.specific.<lineage>.{sample}_busco.txt    # 文本摘要
│   ├── full_table.tsv                # 每个基准基因的详细结果
│   └── missing_busco_list.tsv        # 缺失基因清单
├── busco_summary_results.txt         # 汇总表(所有样本合并, 格式由 --output-format 决定)
├── busco_analysis_summary.txt        # 总结报告(成功率/平均完整度)
└── 99_logs/
    └── busco_analysis.log            # 运行日志
```

## 结果解读 | Interpreting Results { #results }

### 1. 汇总表(`busco_summary_results.txt`)

**通俗理解|In plain words:** 把所有样本的分数排成一张表，每行一个样本，列是完整度/单拷贝/多拷贝/片段化/缺失的百分比和数量。

- **完整度(Complete %)**：越高越好。基因组 >90% 优秀、>80% 合格、<80% 说明组装/注释有明显缺口
- **单拷贝(Single-copy %)**：接近完整度说明没有重复组装问题
- **多拷贝(Duplicated %)**：偏高提示重复组装(杂合组装未去冗余)或近期全基因组复制
- **片段化(Fragmented %)**：偏高提示组装 contig 化严重、有缺口
- **缺失(Missing %)**：高说明大量保守基因没找到，组装或注释质量差

### 2. 总结报告(`busco_analysis_summary.txt`)

**通俗理解|In plain words:** 一页纸的总览——多少样本、成功几个、平均完整度/单拷贝/缺失。快速扫一眼看整体质量。

### 3. 单样本详情(`{sample}_busco/short_summary.*.txt` + `full_table.tsv`)

**通俗理解|In plain words:** 想追查「到底缺了哪些基因」就查这里。`full_table.tsv` 逐基因列出状态(Complete/Duplicated/Fragmented/Missing)和坐标；`missing_busco_list.tsv` 直接列出缺失基因。

## 参数选择建议 | Parameter Guidance { #guidance }

- **基因组完整性**：`-m genome -l <大类>_odb12`，如动物用 metazoa、植物用 embryophyta、真菌用 fungi、普通真核用 eukaryota
- **注释质量(评估 braker.aa)**：`-m proteins -i braker.aa -l <大类>_odb12`
- **不知道选哪个谱系**：用 `--auto-lineage-euk`(真核)或 `--auto-lineage-prok`(原核)自动选
- **目录批量**：`-i dir/ -l ... --sample-suffix '*.fa'` 一次性评估多个基因组
- **要 Excel 汇总**：加 `--output-format xlsx`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入文件或目录｜Input file or directory |
| `--lineage, -l` | 必填 | str | BUSCO数据库谱系｜BUSCO database lineage |
| `--output-dir, -o` | `./busco_output` | Path | 输出目录｜Output directory |
| `--mode, -m` | `genome` | genome/geno/transcriptome/tran/proteins/prot | BUSCO分析模式｜BUSCO analysis mode |
| `--threads, -t` | `12` | int | CPU线程数｜Number of CPU threads |
| `--sample-suffix` | `*.fa` |  | 样本名提取后缀｜Sample name extraction suffix pattern |
| `--output-format` | `txt` | txt/csv/xlsx | 输出文件格式｜Output file format |
| `--force, -f` | — |  | 强制重写现有文件｜Force rewriting existing files |
| `--augustus` | — |  | 使用Augustus基因预测器｜Use Augustus gene predictor |
| `--augustus-parameters` | — | str | Augustus额外参数｜Additional Augustus parameters |
| `--augustus-species` | — | str | Augustus物种名｜Augustus species name |
| `--auto-lineage` | — |  | 自动选择谱系｜Automatically select lineage |
| `--auto-lineage-euk` | — |  | 自动选择真核生物谱系｜Automatically select eukaryote lineage |
| `--auto-lineage-prok` | — |  | 自动选择原核生物谱系｜Automatically select prokaryote lineage |
| `--contig-break` | `10` | int | Contig打断的N数量｜Number of Ns for contig break |
| `--datasets-version` | `odb12` |  | 数据集版本｜Dataset version |
| `--download-path` | — | str | 数据集下载路径｜Dataset download path |
| `--evalue, -e` | `0.001` | float | BLAST E值阈值｜BLAST E-value threshold |
| `--limit` | `3` | int | 候选区域限制｜Candidate region limit |
| `--long` | — |  | 启用Augustus长模式优化｜Enable Augustus long mode optimization |
| `--metaeuk` | — |  | 使用Metaeuk基因预测器｜Use Metaeuk gene predictor |
| `--metaeuk-parameters` | — | str | Metaeuk额外参数｜Additional Metaeuk parameters |
| `--metaeuk-rerun-parameters` | — | str | Metaeuk重新运行参数｜Metaeuk rerun parameters |
| `--miniprot` | — |  | 使用Miniprot基因预测器｜Use Miniprot gene predictor |
| `--skip-bbtools` | — |  | 跳过BBTools统计｜Skip BBTools statistics |
| `--offline` | — |  | 离线模式｜Offline mode |
| `--restart, -r` | — |  | 重启未完成的分析｜Restart incomplete analysis |
| `--quiet, -q` | — |  | 静默模式｜Quiet mode |
| `--scaffold-composition` | — |  | 生成scaffold组成文件｜Generate scaffold composition file |
| `--tar` | — |  | 压缩子目录｜Compress subdirectories |
| `--busco-path` | — |  | BUSCO软件路径｜BUSCO software path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入文件或目录路径｜Input file or directory path |
| `-l, --lineage` | 必填 |  | BUSCO数据库/谱系名称｜BUSCO database lineage name |
| `-o, --output-dir` | `./busco_output` |  | 输出目录｜Output directory |
| `-m, --mode` | `genome` | genome/geno/transcriptome/tran/proteins/prot | BUSCO分析模式｜BUSCO analysis mode |
| `-t, --threads` | `12` | int | CPU线程数｜Number of CPU threads |
| `--sample-suffix` | `*.fa` |  | 样本名称提取后缀模式｜Sample name extraction suffix pattern |
| `--output-format` | `txt` | txt/csv/xlsx | 输出文件格式｜Output file format |
| `-f, --force` | — | store_true | 强制重写现有文件｜Force rewriting of existing files |
| `--augustus` | — | store_true | 使用Augustus基因预测器｜Use Augustus gene predictor |
| `--augustus-parameters` | — |  | Augustus附加参数｜Additional Augustus parameters |
| `--augustus-species` | — |  | Augustus物种名称｜Augustus species name |
| `--auto-lineage` | — | store_true | 自动选择谱系｜Automatically select lineage |
| `--auto-lineage-euk` | — | store_true | 自动选择真核生物谱系｜Automatically select eukaryote lineage |
| `--auto-lineage-prok` | — | store_true | 自动选择原核生物谱系｜Automatically select prokaryote lineage |
| `--contig-break` | `10` | int | Contig断点N数量｜Number of Ns for contig break |
| `--datasets-version` | `odb12` |  | 数据集版本｜Dataset version |
| `--download-path` | — |  | 数据集下载路径｜Dataset download path |
| `-e, --evalue` | `0.001` | float | BLAST E值阈值｜BLAST E-value threshold |
| `--limit` | `3` | int | 候选区域数量限制｜Candidate region limit |
| `--long` | — | store_true | 启用Augustus长模式优化｜Enable Augustus long mode optimization |
| `--metaeuk` | — | store_true | 使用Metaeuk基因预测器｜Use Metaeuk gene predictor |
| `--metaeuk-parameters` | — |  | Metaeuk附加参数｜Additional Metaeuk parameters |
| `--metaeuk-rerun-parameters` | — |  | Metaeuk重运行参数｜Metaeuk rerun parameters |
| `--miniprot` | — | store_true | 使用Miniprot基因预测器｜Use Miniprot gene predictor |
| `--skip-bbtools` | — | store_true | 跳过BBTools统计｜Skip BBTools statistics |
| `--offline` | — | store_true | 离线模式｜Offline mode |
| `-r, --restart` | — | store_true | 重启未完成的分析｜Restart incomplete analysis |
| `-q, --quiet` | — | store_true | 静默模式｜Quiet mode |
| `--scaffold-composition` | — | store_true | 生成scaffold组成文件｜Generate scaffold composition file |
| `--tar` | — | store_true | 压缩子目录｜Compress subdirectories |
| `--busco-path` | — |  | BUSCO软件路径｜BUSCO software path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- **BUSCO**：conda 环境 `busco`(默认 `~/miniforge3/envs/busco/bin/busco`，可用 `--busco-path` 或环境变量 `BUSCO_PATH` 覆盖)
- **BUSCO 谱系数据集**：odb12，需预先下载(离线模式)或联网自动下载
- BUSCO 内部依赖(由 BUSCO 环境提供)：Metaeuk、AUGUSTUS、HMMER、BBTools 等
- **pandas / openpyxl**：结果编译(xlsx 格式需 openpyxl)

## 常见问题 | FAQ { #faq }

**Q1：中断后重跑会从头再来吗？**
不会。每个样本按 `{sample}_busco/short_summary.*.json` 存在性判断，已完成的样本自动跳过。要强制重跑用 `-f`。

**Q2：谱系怎么选？**
选与物种最接近的大类：动物 metazoa_odb12、植物 embryophyta_odb12、真菌 fungi_odb12、不确定就 eukaryota_odb12。选错谱系分数会整体失真。

**Q3：多拷贝比例高说明什么？**
可能有两个原因：一是组装没去冗余(两条单倍型都装进去了)，二是物种近期发生过全基因组复制。前者是组装问题，后者是生物学事实，需要结合物种背景判断。

**Q4：完整度多少算合格？**
经验值：基因组 >90% 优秀，>80% 合格，<80% 说明有明显缺口。但也要看谱系和物种，高杂合/高重复物种本身就更难。

**Q5：`--offline` 是什么，需要联网吗？**
离线模式只用已下载到本地的谱系数据，不联网下载。程序默认开启离线；首次使用前需先把谱系数据下到本地(或用 `--auto-lineage` 联网自动下载)。

**Q6：目录批量为什么有些样本没被评估？**
默认按 `*.fa` 匹配。如果文件名不是 .fa 结尾，用 `--sample-suffix` 指定正确模式(如 `'*.fna'`)；模式里必须含通配符 `*`。
