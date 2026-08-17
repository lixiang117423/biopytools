# LAI - LTR 组装指数 | LTR Assembly Index (LAI)

一句话理解：**输入一个基因组组装，自动识别全长 LTR 逆转座子并算出一个 0~30 的「组装完整度分数(LAI)」**，分数越高说明组装越完整、连续性越好，是评估基因组组装质量的常用指标。

## 功能概述 | Overview { #overview }

- 双工具识别 LTR 候选：LTRharvest（GenomeTools）+ LTR_FINDER_parallel，合并去冗余
- LTR_retriever 从候选中筛出「完整 LTR 逆转座子(intact LTR-RT)」，并做 RepeatMasker 注释
- 用 LAI 脚本按滑动窗口计算全基因组与各染色体的 LAI 值
- 4 种运行模式：`full`(完整) / `harvest`(仅候选识别) / `retrieve`(仅筛选) / `calculate`(仅算 LAI)
- 断点续传：各步骤按输出文件存在性自动跳过（`--skip-completed` 默认开）

## 快速开始 | Quick Start { #quickstart }

```bash
biopytools lai -i genome.fa -o lai_output
```

最小输入：一个基因组 FASTA。结果（含 `{genome}.out.LAI`）写入 `lai_output/`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| LTR 逆转座子 | 最大的一类转座子，两端有重复的「长末端重复」；在植物基因组里占大头 |
| 全长/完整 LTR-RT | 首尾、内部结构都齐全的 LTR 逆转座子；组装越完整，能拼出的全长拷贝越多 |
| LAI(组装指数) | 把「全长 LTR 占总 LTR 的比例」折算成 0~30 的分数，像给组装「完整性打分」 |
| 组装连续性 | 基因组是不是被拼成又长又完整的片段；断掉的部分会把 LTR 也截断 |
| 候选识别(harvest) | 用工具粗筛「可能是 LTR」的片段，先捞一网，不计真假 |
| 筛选(retrieve) | 从候选里挑出真正结构完整的 LTR，丢掉碎片和假阳性 |

## 输入 | Input { #input }

### 基因组 FASTA

必需，标准 FASTA 格式（建议用接近染色体的组装，LAI 才有意义）：

```text
>chr1
ATGCGTACGTACGTAGCTA...
>chr2
TTGCAAGCTAGCATCGATC...
```

## 参数说明 | Parameters { #parameters }

### 运行模式 | Run mode

**通俗理解|In plain words:** `-m/--mode` 决定「跑到哪一步」——默认 `full` 从头跑到底；只想重新识别候选用 `harvest`，只想重新筛选用 `retrieve`，只想重新算分数用 `calculate`。配合断点续传，可以只重跑出错的那一步。**一般用默认 `full`。**

本组参数：`-m/--mode`（默认 full，可选 full/harvest/retrieve/calculate）。

### 运行与续传 | Runtime & resume

**通俗理解|In plain words:** `-t` 线程数（LTR_FINDER/LTR_retriever/LAI 的 BLAST 会吃核）；`--skip-completed` 默认开，已有结果就跳过对应步骤——**换参数重跑前要先删旧产物**，否则会复用。

本组参数：`-t/--threads`（默认 12）、`--skip-completed/--no-skip-completed`（默认开）。

### conda 环境路径 | Conda environments

**通俗理解|In plain words:** 三个软件分属三个 conda 环境，默认路径已配好。**装好后一般不用动**；只有环境装在别处时才需改这三个参数。

本组参数：`--conda-harvest`（默认 `~/miniforge3/envs/ltr_harvest_parallel_v.1.2`）、`--conda-finder`（默认 `~/miniforge3/envs/ltr_finder_parallel_v.1.3`）、`--conda-retriever`（默认 `~/miniforge3/envs/ltr_retriever_v.3.0.1`）。

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 先两个工具各捞一网 LTR 候选 → 合并 → 精选出完整 LTR → 用滑动窗口算 LAI 分数。

```text
基因组 FASTA
    |
    v
步骤1-2: LTRharvest 建索引 + 识别候选  -> {genome}.harvest.scn
          LTR_FINDER_parallel 识别候选  -> {genome}.finder.combine.scn
    |
    v
合并候选                          -> {genome}.rawLTR.scn
    |
    v
步骤3: LTR_retriever 筛选完整 LTR + RepeatMasker 注释
          -> {genome}.pass.list / {genome}.out / {genome}.out.gff
    |
    v
步骤4: LAI 滑动窗口计算               -> {genome}.out.LAI
```

## 输出 | Output { #output }

```text
<输出目录>/
├── {genome}.harvest.scn          # LTRharvest 候选结果
├── {genome}.finder.combine.scn   # LTR_FINDER 候选结果
├── {genome}.rawLTR.scn           # 合并后的候选
├── {genome}.pass.list            # 完整 LTR-RT 列表(核心)
├── {genome}.out                  # RepeatMasker 注释(核心)
├── {genome}.out.gff              # 注释 GFF 格式
├── {genome}.out.LAI              # LAI 结果(核心, 全基因组+各染色体)
├── {genome}.LTRlib.fa            # LTR 库(如 LTR_retriever 生成)
├── tmp/                          # 临时文件
└── lai.log                       # 运行日志
```

其中 `{genome}` 是基因组文件名（含 `.fa` 后缀，如 `genome.fa`）。

## 结果解读 | Interpreting Results { #interpreting }

### LAI 结果文件（`{genome}.out.LAI`）

**通俗理解|In plain words:** 核心输出。文件里按「染色体/窗口」列出各段 LAI，最后一行是 `whole_genome` 的全基因组总分。

```text
Chr  From     To       Intact  Total   Raw_LAI  LAI
chr1 0        3000000  120     150     0.80     24.00
...
whole_genome  0        ...     ...     ...      18.50
```

- **全基因组 LAI 判据**（经验参考）：< 10 为草稿级组装；10-20 为参考级；> 20 为金标准级
- **Raw_LAI / LAI**：前者是「完整 LTR 占比」，后者经总 LTR 含量校正后的最终分数；看 LAI 列即可
- 各染色体分数差异大，提示某些染色体/区段组装明显较差

### 完整 LTR 列表（`{genome}.pass.list`）

完整 LTR-RT 的数量与总长度：数量越多、越完整，说明组装连续性好；数量异常少提示组装碎片化严重。

## 参数选择建议 | Parameter Guidance { #guidance }

- **`--mode`**：首次跑用默认 `full`；中途某步失败，用 `harvest`/`retrieve`/`calculate` 只重跑对应步骤
- **`--threads`**：LAI 计算里的 BLAST 是主要耗时点，建议 12-32；超大基因组可给更多
- **`--skip-completed`**：保持默认开；换参数或换基因组重跑同一目录时，先清空输出目录避免复用旧结果
- **快速预览**：只关心分数、不想要完整流程时可 `-m calculate`（需已有 pass.list 与 .out）

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 | Path | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-m, --mode` | `full` | full/harvest/retrieve/calculate | 运行模式: full(完整流程), harvest(仅候选识别), retrieve(仅筛选), calculate(仅LAI计算)｜Run mode |
| `--skip-completed/--no-skip-completed` | `True` |  | 跳过已完成的步骤｜Skip completed steps |
| `--conda-harvest` | `~/miniforge3/envs/ltr_harvest_parallel_v.1.2` |  | LTR_harvest conda环境路径｜LTR_harvest conda environment path |
| `--conda-finder` | `~/miniforge3/envs/ltr_finder_parallel_v.1.3` |  | LTR_finder conda环境路径｜LTR_finder conda environment path |
| `--conda-retriever` | `~/miniforge3/envs/ltr_retriever_v.3.0.1` |  | LTR_retriever conda环境路径｜LTR_retriever conda environment path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-t, --threads` | `64` | int | 线程数｜Number of threads |
| `-m, --mode` | `full` | full/harvest/retrieve/calculate | 运行模式｜Run mode (default: full) |
| `--conda-env` | `~/miniforge3/envs/ltr_retriever_v.3.0.1` |  | Conda环境路径｜Conda environment path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- LTRharvest（GenomeTools `gt`），conda 环境 `ltr_harvest_parallel_v.1.2`
- LTR_FINDER_parallel，conda 环境 `ltr_finder_parallel_v.1.3`
- LTR_retriever + LAI，conda 环境 `ltr_retriever_v.3.0.1`
- Perl（运行 LTR_retriever 与 LAI 脚本）
- 依赖按 `--mode` 检查：harvest 只查前两个，retrieve 只查 LTR_retriever，calculate 只查 LAI

## 常见问题 | FAQ { #faq }

**Q1：换参数重跑，结果没变？**
断点续传按输出文件存在性判断（harvest.scn / finder.combine.scn / rawLTR.scn / pass.list+out / out.LAI）。换参数或换基因组重跑前，先删除对应旧产物（或清空输出目录）。

**Q2：`calculate` 模式报找不到 pass.list / .out？**
这两个文件由 `retrieve` 步骤生成。先跑 `-m retrieve`（或完整 `full`）拿到它们，再 `-m calculate`。

**Q3：LAI 分数很低是不是组装就废了？**
LAI 主要反映「重复序列区的组装连续性」，低分通常说明组装碎片化、尤其是重复区没拼好；但它是众多评估指标之一，应结合 N50、BUSCO 等一起看。

**Q4：报 conda 环境不存在？**
用 `--conda-harvest` / `--conda-finder` / `--conda-retriever` 指定实际环境路径（模块会自动从环境 bin 目录找 `gt`、`LTR_FINDER_parallel`、`LTR_retriever`、`LAI`）。

**Q5：LTR_FINDER 报输出文件未生成？**
LTR_FINDER_parallel 需在输出目录内运行并生成 `{genome}.finder.combine.scn`；检查输出目录可写、以及该环境里 `LTR_FINDER_parallel` 是否可用。
