# RAxML - 最大似然法建树 | RAxML Maximum Likelihood Phylogeny

一句话理解：**从一份多序列比对或群体 VCF 出发，用「最大似然法」自动找一棵最合理的进化树，并可选做 bootstrap 给每个分叉打上「可信度分数」。**

## 功能概述 | Overview

- 经典 RAxML 最大似然(ML)建树，输入 PHYLIP / FASTA 比对或 **VCF 变异文件(自动转换为 PHYLIP 再建树)**
- 输入格式自动检测(后缀 + 文件头)，也可用 `--input-format` 显式指定
- 默认 GTRGAMMA 替换模型(核酸)、25 个速率类别，适合绝大多数 DNA 数据
- 支持多种算法：`-f d`(默认,快速爬山法)到 `-f a`(rapid bootstrap)等
- 支持 bootstrap 收敛自动停止(autoFC / autoMR / autoMRE 等)
- 输出最佳树、带支持值树、bootstrap 树集合等标准 RAxML 文件
- 多线程加速(需要 PTHREADS 版 RAxML)

## 快速开始 | Quick Start

```bash
biopytools raxml -s variants.vcf -n my_tree
```

最小输入：一个比对文件或 VCF 文件(`-s`，格式自动检测)，一个输出后缀名(`-n`)。结果默认写到 `./raxml_output/`。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| 多序列比对(MSA) | 把多条同源序列「对齐排好」的表格，每一列是同一个「位点」，像把多份手稿对齐到同一行 |
| 系统发育树 | 描述样本之间亲缘关系的「家谱图」，分叉越近亲缘越近 |
| 最大似然(ML) | 一种「找最合理答案」的统计方法：在众多候选树里挑那棵「最可能产生眼前数据」的树 |
| 替换模型 | 描述碱基之间如何互相「变来变去」的一套规则(GTRGAMMA 是核酸最常用的) |
| bootstrap | 「抽查重算」评估可信度：把比对里的位点随机抽样重算很多次，看每个分叉出现多少次 |
| 支持值 | 每个分叉上 bootstrap 的百分比，越高说明这个分叉越可信(如 95 表示 95% 的重复都支持) |
| 外群(outgroup) | 明确「最远亲」的样本，用来给树定方向(定根)，通常选关系最远的物种 |

## 输入 | Input

### VCF 变异文件(自动转换 | auto-converted)

**群体 SNP 的 VCF 可以直接输入**(`.vcf` / `.vcf.gz`)：程序会先调用内置的 vcf2phylip 转换器把基因型变成 PHYLIP 矩阵(双等位 SNP、杂合子按 IUPAC 模糊码编码，如 A/G 杂合记为 `R`)，再自动接着建树。转换产物为 `输出目录/<输出名>_converted.min<m>.phy`，重跑时若该文件已存在会**跳过转换直接建树**。

- 位点会按「最少检出样本数」过滤(默认 4，可用 `--min-samples-locus` 调整)：一个位点至少要有 m 个样本有基因型才保留
- `--resolve-iupac` 可把杂合子随机解析成单一碱基(默认保留 IUPAC 模糊码，RAxML 原生支持)
- InDel 和多等位位点会被跳过，只保留双等位 SNP

### PHYLIP 比对文件

RAxML 经典版要求 **PHYLIP(宽松 relaxed)格式**，首行是两个整数：`序列条数 序列长度`，之后每行一条序列。

```text
5 60
seq1   ACGTACGTACGTACGTACGT...
seq2   ACGTACGTACGTACGTACGT...
seq3   ACGTACGTACGTACGTACGT...
seq4   ACGTACGTACGTACGTACGT...
seq5   ACGTACGTACGTACGTACGT...
```

### FASTA 比对文件

对齐好的 **FASTA 也可以直接输入**（`>` 开头），RAxML 原生支持，程序直接透传使用。

### 格式自动检测

程序按「文件后缀 → 文件首行」自动判断：`.vcf`/`.vcf.gz` 按 VCF 处理；首行 `>` 按 FASTA；首行两个整数按 PHYLIP。判断不符合预期时可用 `--input-format phylip|fasta|vcf` 显式指定。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 就两个必填项：`-s` 告诉程序「比对文件或 VCF 在哪」(格式自动检测)，`-n` 给输出文件起个「后缀名」(不是路径，只是一个名字片段，所有结果文件都会带这个后缀)。

### 输入格式 | Input format

**通俗理解|In plain words:** 这组参数只在「输入是 VCF」或「格式识别不对」时才需要关心。`--input-format` 强制指定输入格式(默认 auto 自动检测，**一般不用动**)；`--min-samples-locus` 是 VCF 转换时「一个位点至少几个样本有基因型才保留」的门槛，样本少(如 <10)可降到 2~3，否则默认 4 即可；`--resolve-iupac` 把杂合子随机拆成纯合，一般保持关闭(RAxML 能直接处理 IUPAC 码)。

### 替换模型 | Substitution model

**通俗理解|In plain words:** 告诉程序「碱基是怎么演化的」。DNA 数据默认 `GTRGAMMA` 基本不用动；只有分析蛋白质序列时才需要换成 `PROTGAMMA...` 一类。`--categories`(速率类别数)和 `--likelihood-epsilon`(优化精度)是细节参数，**一般不用动**。

### 算法与运行次数 | Algorithm and runs

**通俗理解|In plain words:** `-f/--algorithm` 决定「怎么找树」：默认 `d` 是快速爬山法(只找一棵最好树、不做 bootstrap)；想要支持值要换 `a`(rapid bootstrap)。`-N/--runs` 是跑几遍(或 bootstrap 重复次数)。

### Bootstrap | Bootstrap

**通俗理解|In plain words:** 这组参数控制「可信度评估」。**关键点：默认算法 `-f d` 是不做 bootstrap 的**——要拿到支持值，用 `-f a -x <种子> -N <重复数>`(`-f a` 忘给 `-x` 时程序会自动生成一个)。给了 `-x` 就不要再给 `-b`(两种 bootstrap 种子互斥)。`-I/--bootstrap-convergence` 是「自动判断 bootstrap 次数是否够」的省时开关，程序会把它翻译成 RAxML 的 `-# autoXXX` 写法，一般不用动。

### 树选项 | Tree options

**通俗理解|In plain words:** `-t/--starting-tree` 给一个「起始树」让程序从它开始优化；`-g/--constraint-tree` 强制结果树必须符合某棵「约束树」的拓扑；`-o/--outgroup` 指定外群(逗号分隔多个)给树定方向。**没有特殊需求都不用加。**

### 性能与内存 | Performance and memory

**通俗理解|In plain words:** `-T/--threads` 是多线程数(需 PTHREADS 版)。比对特别大内存吃紧时开 `-U/--memory-saving`(内存节省模式,稍慢)。`--disable-pattern-compression` 关闭「相同列压缩」(节省内存的优化)，**一般不用关**。

### 输出与工具 | Output and tool path

**通俗理解|In plain words:** `-w/--output-dir` 指定结果目录。`--raxml-path` 指定 RAxML 程序路径(默认 `raxmlHPC-PTHREADS`，自动解析，缺失时回退 PATH)。`--no-seq-check` 跳过输入检查、`--silent` 静默输出，**一般不用动**。

## 分析流程 | Pipeline

```text
比对文件(PHYLIP/FASTA)或 VCF
    │
    ▼
步骤1: 检测输入格式(VCF 则先转换为 PHYLIP 矩阵, 已转换过则跳过) + 统计序列条数/长度
    │
    ▼
步骤2: 复制序列/起始树/约束树到输出目录
    │
    ▼
步骤3: 构建并执行 RAxML 命令(模型/算法/bootstrap/线程)
    │
    ▼
步骤4: 收集输出文件 + 生成总结报告(<输出名>_summary.txt)
```

## 输出 | Output

```text
raxml_output/
├── RAxML_bestTree.my_tree              # 最佳 ML 树(Newick)
├── RAxML_bipartitions.my_tree          # 带支持值的树(bootstrap 时生成)
├── RAxML_bipartitionsBranchLabels.my_tree  # 分支标签带支持值的树
├── RAxML_bootstrap.my_tree             # bootstrap 树集合(bootstrap 时生成)
├── RAxML_info.my_tree                  # 分析信息(似然值/参数/耗时)
├── RAxML_log.my_tree                   # RAxML 运行日志
├── RAxML_parsimonyTree.my_tree         # 简约法起始树
├── RAxML_result.my_tree                # 结果树
├── my_tree_summary.txt                 # 总结报告(本工具生成)
├── raxml_analysis.log                  # 本工具的运行日志
├── alignment.phy                       # 复制到工作目录的输入比对
└── my_tree_converted.min4.phy          # VCF 输入时自动生成的 PHYLIP 矩阵(仅 VCF 输入)
```

以上文件名中的 `my_tree` 即 `-n` 指定的输出后缀名。

## 结果解读 | Interpreting Results

### 1. 最佳树(RAxML_bestTree)

**通俗理解|In plain words:** 这是最终要用的树，Newick 格式(`(A,(B,C));` 这样的括号串)。用 FigTree、iTOL 等打开即可看到树形图。

### 2. 支持值树(RAxML_bipartitions)

**通俗理解|In plain words:** 只有做了 bootstrap 才会有这个文件。它在每个分叉上标了支持值(0-100 的百分比)，**支持值越高分叉越可信**：一般认为 >70 较可靠、>95 非常可靠。

### 3. 分析信息(RAxML_info)

**通俗理解|In plain words:** 记录本次分析的「体检报告」：最终似然值、模型参数、GAMMA 形状参数 alpha、运行耗时等。似然值本身不直接比较好坏，主要用来核对参数是否符合预期、跑得是否正常。

### 4. 总结报告(my_tree_summary.txt)

本工具额外生成的人类可读总结，汇总了参数、输出文件清单和关键结果，适合归档和论文 Methods 记录。

## 参数选择建议 | Parameter Guidance

- **核酸数据**：`-m GTRGAMMA`(默认)即可；**蛋白质数据**改用 `-m PROTGAMMAWAG` 等蛋白模型
- **群体 VCF 建树**：直接 `-s xxx.vcf.gz`，程序自动转 PHYLIP；样本少时加 `--min-samples-locus 2` 避免位点被过滤太多
- **只要一棵树**：默认 `-f d`(最快)；**要支持值**：`-f a -x <种子> -N 1000`(1000 次 rapid bootstrap 是常用精度)
- **样本多、想让 bootstrap 自动省时**：加 `-I autoMRE` 让程序自己判断 bootstrap 次数是否足够
- **大比对内存吃紧**：加 `-U` 开内存节省模式
- **想固定树的某些关系**：用 `-g` 约束树；**给树定方向**：用 `-o` 外群

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--sequence-file, -s` | 必填 |  | 输入序列文件(PHYLIP/FASTA/VCF, 默认自动检测)｜Input sequence file (PHYLIP/FASTA/VCF, auto-detected by default) |
| `--output-name, -n` | 必填 | str | 输出文件名称｜Output file name |
| `--input-format` | `auto` | auto/phylip/fasta/vcf | 输入格式(auto=自动检测, VCF会先转换为PHYLIP再建树)｜Input format (auto=auto-detect; VCF is converted to PHYLIP before tree building) |
| `--min-samples-locus` | `4` | int | VCF转PHYLIP时每位点最少检出样本数｜Min called samples per locus for VCF to PHYLIP |
| `--resolve-iupac` | `False` |  | VCF转换时随机解析杂合子为单一碱基｜Randomly resolve heterozygotes to single bases in VCF conversion |
| `--model, -m` | `GTRGAMMA` | str | 替换模型｜Substitution model (GTRGAMMA, PROTGAMMAWAG, etc.) |
| `--categories, -c` | `25` | int | 速率异质性类别数｜Number of rate heterogeneity categories |
| `--likelihood-epsilon, -e` | `0.1` | float | 似然优化精度｜Likelihood optimization precision |
| `--algorithm, -f` | `d` | str | 算法类型｜Algorithm type (d=rapid hill-climbing, a=rapid bootstrap, etc.) |
| `--parsimony-seed, -p` | — | int | 简约法随机种子｜Parsimony random seed |
| `--runs, -N` | `1` | str | 运行次数或bootstrap次数｜Number of runs or bootstrap replicates |
| `--bootstrap-seed, -b` | — | int | Bootstrap随机种子｜Bootstrap random seed |
| `--rapid-bootstrap-seed, -x` | — | int | 快速bootstrap随机种子｜Rapid bootstrap random seed |
| `--bootstrap-convergence, -I` | — | autoFC/autoMR/autoMRE/autoMRE_IGN | Bootstrap收敛标准｜Bootstrap convergence criterion |
| `--bootstop-threshold, -B` | `0.03` | float | Bootstrap停止阈值｜Bootstrap stop threshold |
| `--bootstop-perms` | `100` | int | Bootstrap停止检验次数｜Bootstrap stop test permutations |
| `--print-bootstrap-trees, -k` | — |  | 输出带分支长度的bootstrap树｜Print bootstrap trees with branch lengths |
| `--starting-tree, -t` | — | Path | 起始树文件｜Starting tree file |
| `--constraint-tree, -g` | — | Path | 约束树文件｜Constraint tree file |
| `--outgroup, -o` | — | str | 外群名称(逗号分隔多个)｜Outgroup name(s) (comma-separated) |
| `--threads, -T` | `12` | int | 线程数｜Number of threads |
| `--memory-saving, -U` | — |  | 启用内存节省模式｜Enable memory saving mode |
| `--ml-search-convergence, -D` | — |  | 启用ML搜索收敛标准｜Enable ML search convergence criterion |
| `--random-starting-tree, -d` | — |  | 使用随机起始树｜Use random starting tree |
| `--disable-rate-heterogeneity, -V` | — |  | 禁用速率异质性模型｜Disable rate heterogeneity model |
| `--gamma-median, -u` | — |  | 使用GAMMA模型中位数｜Use median for GAMMA model |
| `--disable-pattern-compression, -H` | — |  | 禁用模式压缩｜Disable pattern compression |
| `--output-dir, -w` | `./raxml_output` | Path | 输出目录｜Output directory |
| `--raxml-path` | `raxmlHPC-PTHREADS` | str | RAxML程序路径｜RAxML program path |
| `--no-seq-check` | — |  | 跳过序列检查｜Skip sequence checking |
| `--silent` | — |  | 静默模式｜Silent mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-s, --sequence-file` | 必填 |  | 输入序列文件(PHYLIP/FASTA/VCF, 默认自动检测)｜Input sequence file (PHYLIP/FASTA/VCF, auto-detected by default) |
| `-n, --output-name` | 必填 |  | 输出文件名称｜Output file name |
| `--input-format` | `auto` | auto/phylip/fasta/vcf | 输入格式(默认auto自动检测, VCF会先转换为PHYLIP再建树)｜Input format (default auto-detect; VCF is converted to PHYLIP before tree building) |
| `--min-samples-locus` | `4` | int | VCF转PHYLIP时每位点最少检出样本数｜Min called samples per locus for VCF to PHYLIP |
| `--resolve-iupac` | — | store_true | VCF转换时随机解析杂合子为单一碱基｜Randomly resolve heterozygotes to single bases in VCF conversion |
| `-m, --model` | `GTRGAMMA` |  | 替换模型｜Substitution model (GTRGAMMA, PROTGAMMAWAG, etc.) |
| `-c, --categories` | `25` | int | 速率异质性类别数｜Number of rate heterogeneity categories |
| `-e, --likelihood-epsilon` | `0.1` | float | 似然优化精度｜Likelihood optimization precision |
| `-f, --algorithm` | `d` |  | 算法类型｜Algorithm type (d=rapid hill-climbing, a=rapid bootstrap, etc.) |
| `-p, --parsimony-seed` | — | int | 简约法随机种子｜Parsimony random seed |
| `-#, -N, --runs` | `1` |  | 运行次数或bootstrap次数｜Number of runs or bootstrap replicates |
| `-b, --bootstrap-seed` | — | int | Bootstrap随机种子｜Bootstrap random seed |
| `-x, --rapid-bootstrap-seed` | — | int | 快速bootstrap随机种子｜Rapid bootstrap random seed |
| `-I, --bootstrap-convergence` | — | autoFC/autoMR/autoMRE/autoMRE_IGN | Bootstrap收敛标准｜Bootstrap convergence criterion |
| `-B, --bootstop-threshold` | `0.03` | float | Bootstrap停止阈值｜Bootstrap stop threshold |
| `--bootstop-perms` | `100` | int | Bootstrap停止检验次数｜Bootstrap stop test permutations |
| `-k, --print-bootstrap-trees` | — | store_true | 输出带分支长度的bootstrap树｜Print bootstrap trees with branch lengths |
| `-t, --starting-tree` | — |  | 起始树文件｜Starting tree file |
| `-g, --constraint-tree` | — |  | 约束树文件｜Constraint tree file |
| `-o, --outgroup` | — |  | 外群名称 (逗号分隔多个)｜Outgroup name(s) (comma-separated) |
| `-T, --threads` | `88` | int | 线程数｜Number of threads |
| `-U, --memory-saving` | — | store_true | 启用内存节省模式｜Enable memory saving mode |
| `-D, --ml-search-convergence` | — | store_true | 启用ML搜索收敛标准｜Enable ML search convergence criterion |
| `-d, --random-starting-tree` | — | store_true | 使用随机起始树｜Use random starting tree |
| `-V, --disable-rate-heterogeneity` | — | store_true | 禁用速率异质性模型｜Disable rate heterogeneity model |
| `-u, --gamma-median` | — | store_true | 使用GAMMA模型中位数｜Use median for GAMMA model |
| `-H, --disable-pattern-compression` | — | store_true | 禁用模式压缩｜Disable pattern compression |
| `-w, --output-dir` | `./raxml_output` |  | 输出目录｜Output directory |
| `--raxml-path` | — |  | RAxML程序路径(默认域环境自动解析)｜RAxML program path (default: auto domain env) |
| `--no-seq-check` | — | store_true | 跳过序列检查｜Skip sequence checking |
| `--silent` | — | store_true | 静默模式｜Silent mode |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

| 软件<br>Software | 说明<br>Description |
|------|------|
| raxmlHPC-PTHREADS | RAxML 经典版(多线程 PTHREADS 编译)，经 conda run 自动检测包装，可用 `--raxml-path` 或环境变量 RAXML_PATH 覆盖；未安装于 conda 环境时回退 PATH 直接调用 |

无对应功能域环境(回退 PATH 直接调用 RAxML)。

## 常见问题 | FAQ

### 1. 我给 FASTA 文件可以吗？VCF 可以直接喂吗？

都可以。FASTA 是 RAxML 原生支持的输入，直接透传使用；VCF(`.vcf`/`.vcf.gz`)会自动转换为 PHYLIP 矩阵后再建树，转换产物 `<输出名>_converted.min<m>.phy` 留在输出目录里可复用。

### 2. 跑完没有 RAxML_bipartitions 文件？

因为默认算法 `-f d` 是快速爬山法，**不做 bootstrap**，自然没有支持值树。要支持值请用 `-f a -x <种子> -N <重复数>`(忘给 `-x` 程序会自动生成)。

### 3. 会断点续传吗？

**部分支持。** VCF→PHYLIP 转换有断点续传：转换产物已存在时跳过转换直接建树。RAxML 建树本身每次运行都会重新执行并覆盖同名文件；大数据建议在计算节点一次性跑完。

### 4. 报「找不到 RAxML」？

默认程序名是 `raxmlHPC-PTHREADS`，请确认它已安装且在 PATH 中，或通过 `--raxml-path` 指定完整路径。

### 5. `-n` 是输出路径吗？

不是。`-n` 只是一个「后缀名」，所有结果文件都命名为 `RAxML_xxx.<这个名字>` 并放在 `-w` 指定的目录里。
