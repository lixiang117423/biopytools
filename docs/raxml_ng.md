# RAxML-NG - 最大似然建树 | RAxML-NG Maximum Likelihood Phylogeny

一句话理解：**RAxML 的升级版，一条命令自动完成「模型选择 + bootstrap + 最大似然建树 + 打支持值」，输入 FASTA 或 PHYLIP 比对即可。**

## 功能概述 | Overview

- 三种分析模式：`all`(默认,搜索 + bootstrap + 支持值一次做完)、`search`(只找树)、`support`(给已有树打支持值)
- 进化模型不指定时自动选择(自适应)
- bootstrap 支持固定重复数(如 1000)或 `autoMRE{N}` 自动判断收敛
- 支持多种分支支持值类型(fbp / tbe / sh / ebg / rbs 等)
- 原生断点续传：不传 `--redo` 时自动从上次检查点继续

## 快速开始 | Quick Start

```bash
biopytools raxml-ng -i alignment.fasta -o tree_results
```

最小输入：一个 FASTA(或 PHYLIP)比对(`-i`)和一个输出目录(`-o`)。输出前缀默认用输入文件名。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| 多序列比对(MSA) | 把多条同源序列「对齐排好」的表格，每列是同一「位点」 |
| 最大似然(ML) | 在众多候选树里挑「最可能产生眼前数据」的那棵 |
| 替换模型 | 描述碱基如何互相「变来变去」的规则；让程序自动选最省心 |
| bootstrap | 把位点随机抽样重算很多次，看每个分叉「站得住」多少次 |
| 支持值 | 分叉上的可信度百分比，越高越可信 |
| 外群(outgroup) | 明确「最远亲」的样本，给树定方向 |

## 输入 | Input

### 比对文件

支持 **FASTA** 和 **PHYLIP** 两种格式。FASTA 示例：

```text
>seq1
ACGTACGTACGT...
>seq2
ACGTACGTACGT...
```

PHYLIP 示例(首行为「条数 长度」)：

```text
3 100
seq1   ACGTACGT...
seq2   ACGTACGT...
seq3   ACGTACGT...
```

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 就两个必填项：`-i` 是输入比对，`-o` 是输出目录。输出目录下会自动建 `99_logs` 和 `00_pipeline_info` 两个子目录。

### 核心参数 | Core parameters

**通俗理解|In plain words:** `-p/--prefix` 是输出文件前缀(默认用输入文件名)，所有结果文件都叫 `<前缀>.raxml.*`。`--mode` 决定做多少事：`all` 一次拿全(最常用)，`search` 只找树，`support` 只给已有树打支持值。`-m/--model` 是进化模型，**留空让程序自动选即可，一般不用指定**。`-t/--threads` 是线程数(默认 12)。

### Bootstrap | Bootstrap

**通俗理解|In plain words:** `-b/--bs-trees` 在 `all` 模式下是「bootstrap 重复次数」(默认 1000)；也可写 `autoMRE{N}` 让程序自动判断要跑多少次。`--bs-metric` 是支持值的计算方式(fbp 是默认最常用)，**一般不用动**。

### 树选项 | Tree options

**通俗理解|In plain words:** `--tree` 在 `search` 模式下是「起始树」，在 `support` 模式下是「要打支持值的参考树」。`--outgroup` 指定外群(逗号分隔多个)给树定方向。**没有特殊需求都不用加。**

### 随机与续传 | Random seed and resume

**通俗理解|In plain words:** `--seed` 固定随机种子(让结果可复现，需要时再设)。`--redo` 强制重跑并覆盖已有结果——**不传它就是断点续传**。`--raxml-ng-path` 指定软件路径(默认已配好)。

## 分析流程 | Pipeline

```text
比对文件(FASTA/PHYLIP)
    │
    ▼
步骤1: 检查 RAxML-NG 依赖并取版本号
    │
    ▼
步骤2: 按 mode 执行 RAxML-NG(--all / --search / --support)
    │      (不传 --redo 时从上次检查点续传)
    ▼
步骤3: 写 00_pipeline_info/software_versions.yml
```

`all` 模式等价于 `search`(找最优树) + `bootstrap`(生成 bootstrap 树) + `support`(给最优树标支持值)三合一。

## 输出 | Output

```text
tree_results/
├── <prefix>.raxml.bestTree        # 最佳 ML 树(Newick)
├── <prefix>.raxml.support         # 带支持值的最优树
├── <prefix>.raxml.bootstraps      # bootstrap 树集合
├── <prefix>.raxml.mlTrees         # ML 搜索过程中的候选树
├── <prefix>.raxml.startTree       # 起始树(简约法/随机)
├── <prefix>.raxml.log             # RAxML-NG 运行日志
├── 00_pipeline_info/
│   └── software_versions.yml      # 软件版本与参数记录
└── 99_logs/
    └── raxml_ng_analysis.log      # 本工具运行日志
```

其中 `<prefix>` 是 `-p` 指定的前缀(默认输入文件名)。`support` 模式会额外要求/生成 `<prefix>.raxml.bootstraps`。

## 结果解读 | Interpreting Results

### 1. 最佳树(bestTree)

**通俗理解|In plain words:** 最终要用的树，Newick 格式。用 FigTree、iTOL 打开即见树形图。

### 2. 支持值树(support)

**通俗理解|In plain words:** 带支持值的最优树——每个分叉上的数字就是可信度百分比，**越高越可信**(>70 较可靠、>95 非常可靠)。

### 3. bootstrap 树集合(bootstraps)

**通俗理解|In plain words:** 一堆 bootstrap 重复产生的树，普通用户一般不用看，它是用来算支持值的中间产物。

### 4. 运行日志(raxml.log)

**通俗理解|In plain words:** RAxML-NG 的详细日志，跑挂了先看这个文件，里面会写明错误原因。

## 参数选择建议 | Parameter Guidance

- **绝大多数场景**：默认 `all` 模式 + 模型留空自动选，一次拿全结果
- **想省时间**：bootstrap 用 `-b autoMRE{200}`(最少跑 200 次、自动判断收敛)代替固定 1000
- **已有最优树、只想补支持值**：用 `--mode support --tree <树文件> -b <bootstrap树文件>`
- **要固定随机种子复现结果**：加 `--seed <整数>`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入比对文件｜Input alignment file |
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--prefix, -p` | — |  | 输出文件前缀(默认输入文件名)｜Output prefix (default: input filename) |
| `--mode` | `all` | all/search/support | 分析模式｜Analysis mode |
| `--model, -m` | — |  | 进化模型(不指定则自适应)｜Evolutionary model (auto if not specified) |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--bs-trees, -b` | `1000` |  | Bootstrap重复数或autoMRE{N}｜Bootstrap replicates or autoMRE{N} |
| `--bs-metric` | `fbp` |  | 支持值类型: fbp/tbe/sh/ebg/rbs｜Branch support metric |
| `--tree` | — | Path | 起始树(search)/参考树(support)｜Starting tree (search) / reference tree (support) |
| `--outgroup` | — |  | 外群名称(逗号分隔)｜Outgroup taxon names (comma-separated) |
| `--seed` | — | int | 随机种子｜Random seed |
| `--redo` | — |  | 覆盖已有结果,忽略checkpoint｜Overwrite results, ignore checkpoint |
| `--raxml-ng-path` | — |  | RAxML-NG软件路径｜RAxML-NG program path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入比对文件 (FASTA/PHYLIP)｜Input alignment file (FASTA/PHYLIP) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-p, --prefix` | — |  | 输出文件前缀 (默认输入文件名)｜Output prefix (default: input filename) |
| `--mode` | `all` | all/search/support | 分析模式 (默认: all)｜Analysis mode (default: all) |
| `-m, --model` | — |  | 进化模型 (不指定则自适应)｜Evolutionary model (auto if not specified) |
| `-t, --threads` | `12` | int | 线程数 (默认: 12)｜Number of threads (default: 12) |
| `-b, --bs-trees` | `1000` |  | Bootstrap重复数或autoMRE{N} (默认: 1000)｜Bootstrap replicates or autoMRE{N} |
| `--bs-metric` | `fbp` |  | 支持值类型: fbp/tbe/sh/ebg/rbs (默认: fbp)｜Branch support metric |
| `--tree` | — |  | 起始树(search)/参考树(support)｜Starting tree (search) / reference tree (support) |
| `--outgroup` | — |  | 外群名称 (逗号分隔)｜Outgroup taxon names (comma-separated) |
| `--seed` | — | int | 随机种子｜Random seed |
| `--redo` | — | store_true | 覆盖已有结果,忽略checkpoint｜Overwrite results, ignore checkpoint |
| `--raxml-ng-path` | — |  | RAxML-NG程序路径｜RAxML-NG program path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

| 软件<br>Software | 说明<br>Description |
|------|------|
| raxml-ng | RAxML-NG 静态二进制，默认路径 `~/software/RAxML/raxml-ng_v2.0.2_linux_x86_64/raxml-ng`，可用环境变量 `RAXML_NG_PATH` 覆盖 |

静态二进制不在 conda 环境内时直接调用；若位于 conda 环境则自动用 `conda run` 包装。

## 常见问题 | FAQ

### 1. 会断点续传吗？怎么强制重跑？

会。RAxML-NG 本身支持检查点续传：**不传 `--redo` 时，已完成的步骤自动跳过**。换参数或想从头重跑时，加 `--redo` 覆盖已有结果。

### 2. support 模式报「bootstrap 树文件不存在」？

`support` 模式下 `-b/--bs-trees` 的含义变了：它**必须是已经存在的 bootstrap 树文件**(而不是重复次数)。请先跑一次 bootstrap 得到 `<前缀>.raxml.bootstraps`，再把它传给 `--bs-trees`。

### 3. 模型要不要自己指定？

一般不用。不指定 `-m` 时程序会自动选最合适的替换模型；只有你明确知道该用什么模型时才指定。

### 4. 输出文件在哪？

所有 `<前缀>.raxml.*` 文件都写在 `-o` 指定的输出目录里(RAxML-NG 会把它们写到运行工作目录)。日志在 `99_logs/`，版本信息在 `00_pipeline_info/`。
