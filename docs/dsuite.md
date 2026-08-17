# D统计基因渗入分析 | D-statistic Introgression Analysis (Dsuite)

一句话理解：**用「ABBA-BABA」这个统计检验，判断三个群体(两个亲缘群体 + 一个外群)之间是否发生过基因渗入或杂交**，并估计渗入的比例。

## 功能概述 | Overview

- 封装 Dsuite 的 `Dtrios` 程序，计算 D 统计量（ABBA-BABA 检验）
- 输入一个 VCF 和一份 SETS 分组文件，自动对全部可能的三元组(P1/P2/P3)计算
- 用 Jackknife 估计标准误与 Z 分数，判断 D 值是否显著
- 默认只保留双等位 SNP（D 统计量的标准做法），可切换为 indel / both
- 输出 `BBAA`（D 统计量 + f4-ratio）、`Dmin`（最小 D）、`tree`（按树排列）三种结果
- 无断点续传：每次运行直接重算并覆盖同名输出（换参数重跑需换 `-p` 前缀，见 FAQ）

## 快速开始 | Quick Start

```bash
biopytools dsuite -i variants.vcf.gz -s sets.txt -o output_dir
```

最小输入：一个 VCF（.vcf/.vcf.gz）+ 一份两列 SETS 分组文件（样本 ID、群体 ID）。

## 零基础概念速览 | Concepts in plain words

不熟悉生信术语的话，先花两分钟看这张表，后面的参数说明都会用到：

| 术语 | 通俗理解 |
|------|----------|
| 基因渗入(introgression) | 一个物种的基因「混进」了另一个物种，通常由杂交引起，像一滴墨水溶进另一杯水 |
| ABBA-BABA | 用四个群体(P1、P2、P3、外群)在两个位点上的基因型模式，统计「谁跟谁更像」 |
| D 统计量 | 渗入信号的打分，范围 -1 到 1；0=无渗入，偏离 0 越远越可能有渗入 |
| 外群(Outgroup) | 一个确定「最早分出去」的群体，用来给模式定方向，像参照物/祖辈 |
| f4-ratio | 渗入比例的估计，0-1，越大表示渗入越多 |
| Z 分数 | 显著性检验分数，绝对值大于 3 一般视为显著（近似 p < 0.01） |
| Jackknife | 「留一重抽样」估计误差的方法，用来给 D 值算标准误 |

## 输入 | Input

### VCF 文件

标准 VCF 格式（支持 `.vcf` / `.vcf.gz`），样本名必须与 SETS 文件一致。

### SETS 分组文件

两列制表符分隔：`样本ID` 与 `群体ID`。至少一个样本标为 `Outgroup`（外群）；想忽略某样本，群体写 `xxx`：

```text
Ind1	Species1
Ind2	Species1
Ind3	Species2
Ind4	Species2
Ind5	Species3
Ind6	Outgroup
```

- 每一行的样本名必须能在 VCF 头部找到，找不到的样本会被忽略
- `Outgroup` 用于给 ABBA-BABA 模式「定方向」，没有外群无法计算 D 值

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 说明 |
|------|------|
| `-i, --input` | 输入 VCF 文件（.vcf / .vcf.gz） |
| `-s, --sets` | SETS 分组文件（样本 ID + 群体 ID） |
| `-o, --output-dir` | 输出目录 |

### 等位基因过滤 | Allele filtering

**通俗理解|In plain words:** 决定「什么样的变异位点参与计算」。D 统计量本质是把每个位点分成「两种形态」来计数(ABBA/BABA)，所以标准做法是只保留双等位位点(一个位点只有两种等位基因)；多等位位点(一个位点 3 种及以上等位基因)很难归入两态，默认被排除。**绝大多数项目用默认值即可，几乎不需要动。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--min-alleles` | 2 | 最小等位基因数，低于此值的位点排除 |
| `--max-alleles` | 2 | 最大等位基因数，高于此值的位点排除 |
| `--variant-type` | snps | 变异类型：snps / indels / both / none |

### 其他参数 | Other options

**通俗理解|In plain words:** 这一组管「输出文件名」和「底层工具位置」。`-p` 换输出文件前缀用；`--dsuite-bin`/`--bcftools` 指定底层软件的实际路径；`--collect-stats` 只影响日志详略，不影响计算结果。**一般只需动 `-p` 和 `--dsuite-bin`。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-p, --prefix` | dsuite | 输出文件前缀 |
| `--dsuite-bin` | ~/software/Dsuite/Build/Dsuite | Dsuite 可执行文件路径 |
| `--bcftools` | bcftools | bcftools 命令路径 |
| `--collect-stats` | 关 | 收集并打印 VCF 样本数/变异数统计（仅日志详略） |

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先把 VCF 按等位基因数/类型过滤成干净的位点，再喂给 Dsuite 的 Dtrios 做 ABBA-BABA 检验。整体是「过滤 → 计算」两步，一条管道完成。

```text
输入 VCF + SETS 文件
    |
    v
bcftools 过滤(等位基因数 / 变异类型) —— 可选: 收集 VCF 统计信息
    |
    v
Dsuite Dtrios 计算(所有 P1/P2/P3 三元组的 D 统计量)
    |
    v
输出 BBAA / Dmin / tree 三张结果表 + 运行日志
```

## 输出 | Output

```text
output_dir/
├── dsuite_BBAA.txt                 # D 统计量 + f4-ratio(核心结果)
├── dsuite_Dmin.txt                 # 每个三元组的最小 D 值
├── dsuite_tree.txt                 # 按树结构排列的结果
└── dsuite_analysis_<时间戳>.log    # 运行日志
```

`dsuite` 为默认前缀，可用 `-p` 修改。日志文件名带时间戳，多次运行不会互相覆盖。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看 `dsuite_BBAA.txt` 每一行：P1/P2/P3 是三个群体，D 值偏离 0 越远、Z 分数绝对值越大，越说明 P3 与 P1 或 P2 之间发生过渗入。

- `Dstatistic`：-1 到 1。正/负号表示渗入方向(P3 更接近 P1 还是 P2)，绝对值越大信号越强
- `Z-score`：Jackknife 估计的显著性，绝对值大于 3 通常视为显著
- `p-value`：对应的显著性水平
- `f4-ratio`：渗入比例估计，0-1，越大渗入越多
- `BBAA / ABBA / BABA`：三种基因型模式的计数，D 值由 ABBA、BABA 之差算出
- `Dmin`（在 `dsuite_Dmin.txt`）：每个三元组里 D 的最小值，用于更保守地判断

## 参数选择建议 | Parameter Guidance

- `--min-alleles/--max-alleles`：默认都是 2（只保留双等位 SNP），是 D 统计量的标准做法，**一般不用动**
- `--variant-type`：默认 snps；要分析 indel 再改 `indels` 或 `both`
- `-p/--prefix`：换输出文件名前缀用，默认 dsuite 即可；想保留多组参数的结果，换不同前缀
- `--collect-stats`：想先看看 VCF 有多少样本/变异位点时加上，只影响日志详略

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入VCF文件｜Input VCF file path |
| `--sets, -s` | 必填 |  | SETS分组文件｜SETS file path |
| `--output-dir, -o` | 必填 | Path | 输出目录｜Output directory |
| `--prefix, -p` | `dsuite` |  | 输出文件前缀｜Output file prefix |
| `--dsuite-bin` | `~/software/Dsuite/Build/Dsuite` |  | Dsuite可执行文件路径｜Dsuite binary path |
| `--min-alleles` | `2` | int | 最小等位基因数｜Min number of alleles |
| `--max-alleles` | `2` | int | 最大等位基因数｜Max number of alleles |
| `--variant-type` | `snps` | snps/indels/both/none | 变异类型｜Variant type |
| `--bcftools` | `bcftools` |  | bcftools命令路径｜bcftools command path |
| `--collect-stats` | `False` |  | 收集VCF统计信息｜Whether to collect VCF statistics |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | [FILE] 输入VCF文件路径｜Input VCF file path |
| `-s, --sets` | 必填 |  | [FILE] SETS分组文件路径｜SETS file path |
| `-o, --output-dir` | 必填 |  | [DIR] 输出目录｜Output directory |
| `-p, --prefix` | `dsuite` |  | [STR] 输出文件前缀｜Output file prefix |
| `--dsuite-bin` | `~/software/Dsuite/Build/Dsuite` |  | [FILE] Dsuite可执行文件路径｜Dsuite binary path |
| `--min-alleles` | `2` | int | [INT] 最小等位基因数｜Min number of alleles |
| `--max-alleles` | `2` | int | [INT] 最大等位基因数｜Max number of alleles |
| `--variant-type` | `snps` | snps/indels/both/none | [STR] 变异类型｜Variant type |
| `--bcftools` | `bcftools` |  | [CMD] bcftools命令路径｜bcftools command path |
| `--collect-stats` | — | store_true | [FLAG] 是否收集VCF统计信息｜Whether to collect VCF statistics |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Dsuite（默认 `~/software/Dsuite/Build/Dsuite`，需自行编译安装，用 `--dsuite-bin` 指定实际路径）
- bcftools（按等位基因数/变异类型过滤 VCF）

## 常见问题 | FAQ

**Q1：为什么样本里必须有一个 Outgroup？**
外群用来给「基因型模式」定方向(判断谁是祖先状态)。没有 Outgroup，ABBA-BABA 的 D 值无从计算。

**Q2：VCF 样本名和 SETS 文件对不上会怎样？**
SETS 里写了但 VCF 里没有的样本会被忽略；VCF 里有但 SETS 没写的样本不会被纳入任何三元组。两边样本名务必完全一致。

**Q3：D 值显著就一定是基因渗入吗？**
不一定。群体结构、测序错误、祖先群体大小变化等也可能让 D 显著偏离 0。显著 D 是「提示有渗入」，还需结合生物学背景(杂交记录、地理分布)判断。

**Q4：换过滤参数重跑要删旧文件吗？**
Dsuite 模块无断点续传，每次运行直接重算并覆盖同名输出。换 `--variant-type` 等参数重跑同一 `-o`/`-p` 即可；想保留多组结果请换 `-p` 前缀。

**Q5：提示「Dsuite not found」怎么办？**
程序会检查 `--dsuite-bin` 指向的文件是否存在。默认路径不可用时，用 `--dsuite-bin /绝对路径/Dsuite` 传入编译好的 Dsuite 可执行文件的实际绝对路径。