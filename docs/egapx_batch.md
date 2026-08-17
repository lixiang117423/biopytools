# egapx-batch - EGAPx 批量运行配置生成 | EGAPx batch config generator

一句话理解：**它自己不跑注释，而是把一个大基因组按染色体拆开，给每条染色体生成一份 EGAPx 注释的配置文件和运行脚本，让你能一条条并行提交**。
它解决的问题：EGAPx(NCBI 的真核基因组注释流程)对超大基因组或单条超长染色体不友好，拆成按染色体的多个小任务就能并行跑、也更容易容错。

## 功能概述 | Overview { #overview }

- 按染色体(或序列)拆分基因组，为每条序列生成独立的 EGAPx YAML 配置 + 运行脚本
- 生成一个总的提交列表脚本(`all_jobs_submit.list.sh`)，支持串行 / GNU parallel / xargs 三种方式并行执行
- 支持短读、长读测序数据作为注释证据，自动生成 EGAPx 格式的 reads 列表
- 在输出目录上层自动创建 EGAPx 运行所需的软链接(ui/nf/egapx_config)与 SIF 镜像配置
- 可自定义 locus 标签前缀、报告名、染色体前缀过滤、NCBI 物种分类 ID
- 注意：本模块只**生成**配置，真正的注释由用户随后执行生成的脚本完成

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools egapx-batch -g genome.fa -o output_dir
```

最小输入：一个基因组 FASTA + 一个输出目录。生成后按提示执行 `bash all_jobs_submit.list.sh`(或 parallel 并行)开始真正注释。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| EGAPx | NCBI 出品的真核基因组自动注释流程，一条龙完成基因预测到提交报告 |
| 染色体拆分 | 把一个大基因组按每条序列拆成小文件，让注释任务变小、能并行 |
| YAML 配置 | 给 EGAPx 的任务说明书，写清基因组、证据、物种 ID、前缀等 |
| locus 标签 | 给每个预测基因起的唯一编号前缀，如 `PSOJA_` 这种 |
| taxid | NCBI 物种分类编号，告诉注释软件「这是哪个物种」 |
| Singularity(.sif) | 打包好的软件镜像，EGAPx 用它保证依赖一致 |
| 并行提交 | 同时跑多个染色体的任务，省时间 |

## 输入 | Input { #input }

### 基因组(-g，必填)

FASTA 格式，支持 `.fa` / `.fa.gz` / `.fasta` / `.fasta.gz`。按序列(通常是一条染色体)拆分，序列头里的特殊字符(冒号、竖线、斜杠、反斜杠)会自动替换成下划线 `_`。

### 测序数据(可选)

- `--short-reads` / `--long-reads`：短读/长读数据，文件或目录，支持 `.fq` / `.fq.gz` / `.fastq` / `.fastq.gz`。目录会扫描并按样本名(自动识别 `_1/_2`、`_R1/_R2` 等双端后缀)分组生成 EGAPx 格式 reads 列表

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 基因组 + 输出目录。输出目录会按染色体名建子目录，每个子目录放该染色体的配置和脚本。

### EGAPx 路径 | EGAPx paths

**通俗理解|In plain words:** `-e`(EGAPx 安装目录)、`--local-cache`(本地缓存)、`--sif`(Singularity 镜像)三个路径默认指向 `~/software/EGAPX_v.0.4.1-alpha/` 下的标准位置。**一般不用动**，除非你的 EGAPx 装在别处。

### 拆分与过滤 | Split & filter

**通俗理解|In plain words:** 默认按序列拆分(`--no-split` 关闭拆分、整个基因组一个任务)。`-p`/`--chr-prefix` 只保留序列 ID 以某前缀开头的序列(如只想注释核染色体、跳过 organelle/未定位 scaffold)。**一般不用动，默认全拆**。

### 命名与物种 | Naming & taxonomy

**通俗理解|In plain words:** `--locus-prefix` 是预测基因编号的前缀(最终会拼成 `<prefix>_<chr>`)，`--report-name` 是报告名，`--taxid` 是 NCBI 物种分类 ID。**taxid 默认 71234 是占位值，务必改成自己物种的 NCBI taxonomy ID**，否则注释报告里的物种信息是错的。

## 分析流程 | Pipeline { #pipeline }

```text
输入: 基因组 FASTA [+ 短读/长读数据]
    |
    v
步骤1: 拆分基因组(默认) -> 每条序列一个 .fa 文件(临时目录 temp_split)
    |
    v
步骤2: 为每条序列生成配置与脚本
  - {chr}.fa: 该序列的 FASTA
  - {chr}.yaml: EGAPx 配置(genome/taxid/locus前缀/reads/线程)
  - egapx_{chr}.sh: 运行脚本(激活 conda + singularity 跑 EGAPx)
  - 追加到 all_jobs_submit.list.sh
    |
    v
步骤3: 清理临时目录 + 生成 EGAPx 软链接 + 打印执行方式
  -> 用户执行 bash all_jobs_submit.list.sh 开始注释
```

## 输出 | Output { #output }

```text
output_dir/
├── egapx_batch.log                   # 生成日志
├── <chr1>/                           # 每条染色体一个子目录
│   ├── <chr1>.fa                     # 该染色体序列
│   ├── <chr1>.yaml                   # EGAPx 配置
│   ├── egapx_<chr1>.sh               # 运行脚本(可执行)
│   ├── work/                         # EGAPx 工作目录
│   └── output/                       # EGAPx 输出目录
├── <chr2>/ ...                       # 其余染色体同结构
├── short_reads_list.txt              # 短读列表(给了 --short-reads 才有)
└── long_reads_list.txt               # 长读列表(给了 --long-reads 才有)

output_dir 的上一级目录:
├── all_jobs_submit.list.sh           # 总提交列表(逐行 bash 每条染色体的脚本)
├── ui/ nf/ egapx_config/             # EGAPx 运行所需软链接(自动创建)
└── (egapx_config/singularity.config) # 自定义 SIF 镜像配置(给了 --sif 时)
```

## 结果解读 | Interpreting Results { #results }

### 1. 生成物是什么

**通俗理解|In plain words:** 本模块的「结果」是配置和脚本，不是生物学注释。每个 `<chr>.yaml` 是该染色体的任务说明书，`egapx_<chr>.sh` 是执行它的命令。

### 2. 执行方式(日志末尾会打印)

- 串行：`bash all_jobs_submit.list.sh`(一条条跑，最慢但最稳)
- 并行：`cat all_jobs_submit.list.sh | parallel -j 4` 或 `xargs -P 4`

### 3. 真正的注释结果在哪

执行脚本后，注释结果在各 `<chr>/output/` 里；全部跑完后把各染色体的 GFF 合并，即得到全基因组注释。

## 参数选择建议 | Parameter Guidance { #guidance }

- **常规用法**：`-g genome.fa -o out --taxid <你的物种ID>`，其余默认
- **有 RNA-seq 证据**：加 `--short-reads <dir>`(有长读再 `--long-reads`)
- **只想注释核染色体**：`-p chr` 只保留前缀为 chr 的序列(视具体命名)
- **基因组小、不想拆**：`--no-split` 整体一个任务
- **统一基因命名**：`--locus-prefix 物种缩写`，预测基因会编号成 `物种缩写_染色体_序号`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--genome, -g` | 必填 |  | 基因组FASTA文件路径｜Genome FASTA file path |
| `--output, -o` | 必填 | Path | 输出目录路径｜Output directory path |
| `--egapx, -e` | `~/software/EGAPX_v.0.4.1-alpha/egapx` |  | EGAPx安装路径｜EGAPx installation path |
| `--local-cache` | `~/software/EGAPX_v.0.4.1-alpha/local_cache` |  | EGAPx本地缓存路径｜EGAPx local cache path |
| `--sif` | `~/software/EGAPX_v.0.4.1-alpha/egapx/egapx_0.4.1-alpha.sif` |  | Singularity镜像路径｜Singularity image path |
| `--no-split` | `False` |  | 不按染色体拆分基因组｜Do not split genome by chromosome |
| `--chr-prefix, -p` | — |  | 染色体前缀过滤｜Chromosome prefix filter |
| `--locus-prefix` | `` |  | locus标签前缀｜Locus tag prefix |
| `--report-name` | `EGAPx` |  | 报告名称｜Report name |
| `--short-reads` | `` |  | 短读长测序数据(目录或文件)｜Short reads (directory or file) |
| `--long-reads` | `` |  | 长读长测序数据(目录或文件)｜Long reads (directory or file) |
| `--taxid` | `71234` |  | 物种分类ID｜Species taxonomy ID |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | [FILE] 基因组FASTA文件路径｜Genome FASTA file path |
| `-o, --output` | 必填 |  | [DIR] 输出目录路径｜Output directory path |
| `-e, --egapx` | `~/software/EGAPX_v.0.4.1-alpha/egapx` |  | [PATH] EGAPx安装路径｜EGAPx installation path |
| `--local-cache` | `~/software/EGAPX_v.0.4.1-alpha/local_cache` |  | [PATH] EGAPx本地缓存路径｜EGAPx local cache path |
| `--sif` | `~/software/EGAPX_v.0.4.1-alpha/egapx/egapx_0.4.1-alpha.sif` |  | [FILE] Singularity镜像路径｜Singularity image path |
| `--no-split` | `False` | store_true | 不按染色体拆分基因组｜Do not split genome by chromosome |
| `--taxid` | `71234` |  | [INT] 物种分类ID｜Species taxonomy ID |
| `-p, --chr-prefix` | — |  | [STR] 染色体前缀过滤｜Chromosome prefix filter |
| `--locus-prefix` | `` |  | [STR] locus标签前缀｜Locus tag prefix |
| `--report-name` | `EGAPx` |  | [STR] 报告名称｜Report name |
| `--short-reads` | `` |  | [FILE] 短读长测序数据文件路径｜Short reads file path |
| `--long-reads` | `` |  | [FILE] 长读长测序数据文件路径｜Long reads file path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- **awk**：基因组拆分用(系统自带)
- **EGAPx**：默认 `~/software/EGAPX_v.0.4.1-alpha/egapx`(含 ui/nf/egapx_config)，**需预先安装**
- **Singularity** 与 EGAPx SIF 镜像：默认 `~/software/EGAPX_v.0.4.1-alpha/egapx/egapx_0.4.1-alpha.sif`
- **conda(EGAPx 环境)**：生成的脚本会激活 `EGAPx_v.0.4.0-alpha` 环境与 JAVA_HOME
- **python3**：EGAPx 运行时调用(脚本内 `python3 ui/egapx.py`)

## 常见问题 | FAQ { #faq }

**Q1：这个工具会直接跑注释吗？**
不会。它只生成配置和脚本，真正的注释要你随后执行 `all_jobs_submit.list.sh`(或 parallel)才会开始。

**Q2：生成的 YAML 里 threads 写的是 88，对吗？**
对。生成器给每条染色体任务默认写 `threads: 88`(针对大染色体)。如需调整，改对应 `<chr>.yaml` 后再跑脚本。

**Q3：taxid 默认 71234 是什么？**
是生成器内置的占位物种 ID。务必用 `--taxid` 改成自己物种的 NCBI taxonomy ID，否则报告里的物种信息不对。

**Q4：有断点续传吗？**
本模块是轻量生成器，没有断点续传概念——重跑就是重新生成配置。真正注释阶段由 EGAPx 自己在 `<chr>/work/` 里管理续传。

**Q5：为什么输出目录上层多了 ui/nf/egapx_config 软链接？**
EGAPx 运行时需要这些目录，生成器自动把它们软链接到输出目录的上一级，供脚本调用。这是正常现象，别删。

**Q6：拆分后序列名变了怎么办？**
序列 ID 里的特殊字符(冒号、竖线、斜杠、反斜杠)会被替换成下划线 `_`，避免下游出错。若某序列不想拆分，用 `-p <前缀>` 过滤或 `--no-split`。
