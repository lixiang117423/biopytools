# EDTA - 转座子注释 | EDTA Transposable Element Annotation

一句话理解：**给一个基因组 FASTA，自动把里面的「转座子」(会跳跃复制的 DNA 片段) 找出来、分类、并标注位置**，产出一个非冗余 TE 库和全基因组注释文件，供下游基因预测前屏蔽重复序列。

## 功能概述 | Overview { #overview }

- 封装 EDTA（Extensive de-novo TE Annotator），单基因组转座子鉴定、分类与注释一条龙
- 内置 Rice / Maize 物种预设（用对应 CDS 做误检过滤），其余物种用 `others`
- 步骤可拆分：`all`(完整) / `filter`(筛选) / `final`(汇总) / `anno`(注释)，方便中途重启
- 支持喂入 CDS、已筛选 TE 库(curatedlib)、RepeatModeler 库(rmlib) 等辅助输入提高精度
- `--anno 1` 产出全基因组 TE 注释(GFF3 + 统计汇总表)；`--evaluate 1` 可评估注释一致性
- 依赖 EDTA 自身做断点续传（靠 `--step`/`--overwrite` 控制），换参数重跑需覆盖旧产物

## 快速开始 | Quick Start { #quickstart }

```bash
biopytools edta -i plant.fa --anno 1 -t 24
```

最小输入：一个基因组 FASTA。产物会写在**基因组 FASTA 所在目录**（不是 `-o` 目录），详见「输出」。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 转座子/TE | 基因组里会「复制并跳到别处」的 DNA 片段，像会自动繁殖的段落；基因预测前通常要先把它们标记出来 |
| LTR 逆转座子 | 最大的一类转座子，两端有重复的「长末端重复」，像一段话首尾各抄一遍作为识别标记 |
| TIR / LINE / SINE / Helitron | 其他几类转座子，各有各的「跳法」；EDTA 会分别用不同方法找它们 |
| TE 库(TElib) | 把找到的转座子「去重汇总」成一份代表序列文件，可复用于屏蔽其它基因组 |
| 注释(annotation) | 把 TE 在基因组上的精确坐标、方向、类型写成 GFF3 表格 |
| 屏蔽(masking) | 把重复区「糊掉/标记」使基因预测软件不去那里找基因 |
| 分歧度(divergence) | 一个转座子拷贝与「祖先序列」差了多少；越高越老、越退化 |
| CDS | 已知基因的编码序列，用来把「混进 TE 库里的真基因」过滤掉 |

## 输入 | Input { #input }

### 基因组 FASTA

必需，标准 FASTA 格式（可含多序列）。示例：

```text
>chr1
ATGCGTACGTACGTAGCTA...
>chr2
TTGCAAGCTAGCATCGATC...
```

### 可选输入文件

- `--cds`：CDS 序列 FASTA（推荐非水稻/玉米物种提供，用于去除假 TE）
- `--curatedlib`：已人工筛选的可信 TE 库 FASTA（跳过从头鉴定，直接用它）
- `--rmlib`：RepeatModeler 产生的库 FASTA
- `--rmout`：已有 RepeatMasker 输出(.out)，跳过 RepeatMasker 步骤
- `--exclude`：BED 格式，标记不参与注释的区域（如线粒体/叶绿体）

## 参数说明 | Parameters { #parameters }

### 物种与步骤 | Species & step

**通俗理解|In plain words:** `species` 决定「用哪套 CDS 去过滤假 TE」——Rice/Maize 用内置 CDS，其它物种选 `others` 并建议配合 `--cds`。`step` 决定「跑到哪一步」，是 EDTA 断点续传的开关；**一般保持默认 `all`**。

本组参数：`--species`（默认 others，可选 Rice/Maize）、`--step`（默认 all，可选 all/filter/final/anno）、`--overwrite`（默认 0，置 1 覆盖已有产物）。

### 辅助库与过滤 | Auxiliary libraries

**通俗理解|In plain words:** 这组参数是「纠错与省事」用的——提供已知 CDS/库文件能显著降低假阳性；`--sensitive 1` 让程序额外跑一遍 RepeatModeler（更全但更慢）；`--exclude` 把不想碰的区域划掉。**没这些文件时留空即可，不影响运行。**

本组参数：`--cds`、`--curatedlib`、`--rmlib`、`--rmout`、`--sensitive`（默认 0）、`--exclude`。

### 注释与评估 | Annotation & evaluation

**通俗理解|In plain words:** `--anno 1` 是「要不要最后的全基因组注释」——要拿到 GFF3 坐标就必须开；`--evaluate 1` 额外做一轮注释一致性自检（更慢）。`--force 1` 强制继续过滤步骤（自动填充缺失的 TE 类型文件），只在程序因缺某些 TE 类型报错、想让它自动补空文件续跑时才用。**平时 `--anno 1` 即可，其余一般不用动。**

本组参数：`--anno`（默认 0）、`--evaluate`（默认 0）、`--force`（默认 0）。

### 分歧度与突变率 | Divergence & mutation rate

**通俗理解|In plain words:** `maxdiv` 是「TE 与祖先差多少才算 TE」的上限（0-100），调大=纳入更老的退化拷贝（更全但更脏），调小=只留年轻完整的（更干净但更少）。`u` 是中性突变率，用来把分歧度换算成「年龄」。**默认 40 / 1.3e-8 适用大多数植物，一般不动。**

本组参数：`--maxdiv`（默认 40）、`--u`（默认 1.3e-8）。

### 运行与环境 | Runtime & environment

**通俗理解|In plain words:** `-t` 是线程数（越大越快，RepeatMasker/RepeatModeler 等会吃核）；`--debug 1` 保留中间文件（排查问题用，占磁盘）；`-o` 是输出目录但**只放日志**；`--edta-path` 指向 EDTA 安装目录（含 EDTA.pl 的目录），找不到 EDTA 时用它指定。

本组参数：`-t/--threads`（默认 12）、`--debug`（默认 0）、`-o/--output-dir`（默认 ./edta_output）、`--edta-path`。

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** EDTA 内部按「找候选 → 筛假阳性 → 汇总成库 → 全基因组注释」四步走，`--step` 就是让用户从任意一步切入。

```text
基因组 FASTA
    |
    v
raw    : 各方法分别找候选(LTR/TIR/Helitron/LINE/SINE)  -> {genome}.mod.EDTA.raw/
    |
    v
filter : 用 CDS 等过滤假阳性、去冗余
    |
    v
final  : 汇总成非冗余 TE 库                        -> {genome}.mod.EDTA.TElib.fa
    |
    v
anno   : RepeatMasker 全基因组注释(需 --anno 1)   -> {genome}.mod.EDTA.TEanno.gff3
```

## 输出 | Output { #output }

EDTA 原生产物运行结束后自动归集到 `-o` 输出目录的 `01_edta_raw/`（by-step 结构，`{genome}` 为基因组文件名剥 `.fa` 后缀，如 `plant`）：

```text
-o 输出目录/
├── 00_pipeline_info/
│   └── software_versions.yml         # EDTA 与模块版本
├── 01_edta_raw/                      # EDTA 原生产物(运行结束自动归集到这里)
│   ├── {genome}.mod.EDTA.raw/        # 各类型原始候选(中间产物, --debug 0 可能被清理)
│   │   ├── {genome}.mod.LTR.raw.fa
│   │   ├── {genome}.mod.TIR.raw.fa
│   │   ├── {genome}.mod.Helitron.raw.fa
│   │   └── ...
│   ├── {genome}.mod.EDTA.TElib.fa    # 非冗余 TE 库(核心产物)
│   ├── {genome}.mod.EDTA.TElib.novel.fa # 新鉴定 TE(不属于已知库的部分)
│   ├── {genome}.mod.EDTA.TEanno.gff3 # 全基因组 TE 注释(--anno 1)
│   └── {genome}.mod.EDTA.TEanno.sum  # 注释统计汇总(--anno 1)
└── 99_logs/
    └── edta_analysis.log             # 运行日志
```

注意：运行过程中产物先写在输出目录根（EDTA 以运行目录为写入位置），成功（或 filter 容错续跑成功）后自动归集到 `01_edta_raw/`；归集幂等，重跑安全。

## 结果解读 | Interpreting Results { #interpreting }

### 1. TE 库（`{genome}.mod.EDTA.TElib.fa`）

**通俗理解|In plain words:** 这就是「这个基因组里转座子的代表名单」——去重后的每条序列代表一类转座子，可复用于屏蔽其它基因组或喂给 RepeatMasker。

- 序列 ID 通常带类型标签（如 LTR/Copia、LTR/Gypsy、TIR/...），条数越多代表 TE 多样性越高
- `novel.fa` 是其中「新鉴定」的部分，说明这个物种有参考库里没有的转座子

### 2. 注释表（`TEanno.gff3` 与 `TEanno.sum`）

**通俗理解|In plain words:** GFF3 是「每个转座子在基因组上的坐标+类型」清单；sum 是「按类型统计的数量/长度/占比」汇总表。

- `TEanno.sum` 会按 Order/Superfamily 统计 TE 数量与总长度，**看总 TE 占基因组比例**是核心指标（多数植物 30%-80%）
- 异常偏低（如 <10%）提示可能漏检（如没提供 CDS、或 `--maxdiv` 太小）；异常偏高需排查是否把基因区也标成 TE

### 3. 一致性评估（`--evaluate 1` 时）

评估不同方法对同一 TE 的注释一致性，一致性高说明注释可信；分歧大提示某些区域可能是假阳性。

## 参数选择建议 | Parameter Guidance { #guidance }

- **`--species`**：水稻选 `Rice`、玉米选 `Maize`（内置 CDS 过滤）；其它物种选默认 `others` 并加 `--cds your_cds.fa`
- **`--anno`**：只要需要 GFF3 坐标就置 1；只想拿 TE 库可先 `--step final` 快速出库
- **`--sensitive 1`**：从头注释更全面，但会额外跑 RepeatModeler，耗时明显增加；默认关闭对多数项目已足够
- **`--maxdiv`**：年轻完整 TE 为主可降到 20-30；想尽量找全（含古老退化拷贝）可保持 40 或略升
- **`--overwrite 1`**：换参数重跑同一基因组时置 1，否则 EDTA 复用旧产物
- **`--threads`**：EDTA 内部多步会吃核，建议 24-64 加速

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--species` | `others` | Rice/Maize/others | 物种类型｜Species type |
| `--step` | `all` | all/filter/final/anno | 运行步骤｜Step to run |
| `--overwrite` | `0` | IntRange | 覆盖已有结果｜Overwrite existing results |
| `--cds` | — |  | CDS序列文件｜CDS sequences file |
| `--curatedlib` | — |  | 筛选TE库｜Curated TE library |
| `--rmlib` | — |  | RepeatModeler库｜RepeatModeler library |
| `--sensitive` | `0` | IntRange | 敏感模式RepeatModeler建库｜Sensitive-mode RepeatModeler library |
| `--anno` | `0` | IntRange | 执行全基因组注释｜Perform whole-genome annotation |
| `--rmout` | — |  | RepeatMasker输出文件｜RepeatMasker output file |
| `--maxdiv` | `40` | int | 最大分歧度｜Maximum divergence |
| `--evaluate` | `0` | IntRange | 评估注释一致性｜Evaluate annotation consistency |
| `--exclude` | — |  | 排除区域BED文件｜Exclude regions BED file |
| `--force` | `0` | IntRange | 强制继续过滤步骤(自动填充缺失TE文件)｜Force filtering (auto-fill missing TE files) |
| `--u` | `1.3e-08` | float | 中性突变率｜Neutral mutation rate |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--debug` | `0` | IntRange | 保留中间文件｜Retain intermediate files |
| `-o, --output-dir` | `./edta_output` | Path | 输出目录｜Output directory |
| `--edta-path` | — | Path | EDTA安装路径｜EDTA installation path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --genome` | 必填 |  | 基因组FASTA文件｜Genome FASTA file |
| `--species` | `others` | Rice/Maize/others | 物种类型｜Species type |
| `--step` | `all` | all/filter/final/anno | 运行步骤｜Step to run |
| `--overwrite` | `0` | 0/1 | 覆盖已有结果｜Overwrite existing results |
| `--cds` | — |  | CDS序列文件｜CDS sequences file |
| `--curatedlib` | — |  | 筛选TE库｜Curated TE library |
| `--rmlib` | — |  | RepeatModeler库｜RepeatModeler library |
| `--sensitive` | `0` | 0/1 | 敏感模式RepeatModeler建库(同源阈值放宽)｜Sensitive-mode RepeatModeler library |
| `--anno` | `0` | 0/1 | 执行全基因组注释｜Perform whole-genome annotation |
| `--rmout` | — |  | RepeatMasker输出文件｜RepeatMasker output file |
| `--maxdiv` | `40` | int | 最大分歧度｜Maximum divergence |
| `--evaluate` | `0` | 0/1 | 评估注释一致性｜Evaluate annotation consistency |
| `--exclude` | — |  | 排除区域BED文件｜Exclude regions BED file |
| `--force` | `0` | 0/1 | 强制继续过滤步骤(自动填充缺失TE文件)｜Force filtering (auto-fill missing TE files) |
| `--u` | `1.3e-08` | float | 中性突变率｜Neutral mutation rate |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--debug` | `0` | 0/1 | 保留中间文件｜Retain intermediate files |
| `-o, --output-dir` | `./edta_output` |  | 输出目录｜Output directory |
| `--edta-path` | — |  | EDTA安装路径｜EDTA installation path |
| `-i, --genome-list` | 必填 |  | 基因组列表文件｜Genome list file |
| `-c, --cds` | — |  | CDS序列文件｜CDS sequences file |
| `-l, --curatedlib` | — |  | 筛选TE库｜Curated TE library |
| `-f, --fl-copy` | `3` | int | 全长拷贝数阈值｜Full-length copy number cutoff |
| `-a, --anno` | `1` | 0/1 | 执行全基因组注释｜Perform whole-genome annotation |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Perl + EDTA（默认查找 `~/miniforge3/envs/EDTA_v.2.2.2/share/EDTA/EDTA.pl`，也可 `--edta-path` 指定；conda 环境名 `EDTA_v.2.2.2`）
- EDTA 内部会调用：RepeatMasker、RepeatModeler、LTR_retriever、CD-HIT、BLAST 等（随 EDTA 环境提供）
- 无 conda 环境自动检测：`get_conda_env` 会按 PATH 自动包装 `conda run`

## 常见问题 | FAQ { #faq }

**Q1：结果输出到哪了？`-o` 目录里为什么只有日志？**
EDTA 默认把产物写在**基因组 FASTA 所在目录**（如 `plant.fa.mod.EDTA.TElib.fa`）。`-o` 目录只存放 `edta_analysis.log` 运行日志。

**Q2：换参数重跑，结果没变？**
EDTA 有断点续传，复用旧中间文件。换过滤参数重跑时加 `--overwrite 1` 覆盖；或先用 `--step filter`/`--step final` 指定切入步骤。

**Q3：报「EDTA 路径不存在」或找不到 EDTA.pl？**
用 `--edta-path /path/to/EDTA` 指定 EDTA 安装目录（该目录下应有 `EDTA.pl`）。

**Q4：运行中报某 TE 类型（如 SINE/LINE/TIR/Helitron）not found？**
说明该物种缺少这类转座子的完整结构。包装器会自动创建空文件并用 `--force 1` 继续跑 filter 步骤，属正常容错。

**Q5：非模式物种没给 `--cds`，注释假阳性很多怎么办？**
提供该物种的 CDS 文件（`--cds`），EDTA 会用它把混进 TE 库的真基因过滤掉，显著降低假阳性。

## 已知限制与预检 | Known limits & prechecks

- **序列 ID ≤13 字符**：EDTA/RepeatMasker 对序列 ID 有 13 字符硬限（如 `HiC_scaffold_100` 会撞限），模块在 `validate` 阶段预检并直接报错提示改名（`biopytools chr_rename` / `assembly_qc`），避免数小时运行中途失败
- **输出位置**：EDTA 原生结果（`{genome}.mod.EDTA.*`）生成在模块输出目录（EDTA 以运行目录为写入位置），不写在基因组旁
- **环境**：默认走 `edta_v.2.3.0` 域解析（env_map 注册 `EDTA.pl`/`panEDTA.sh`）；可用 `--edta-path` 或 `EDTA_PATH`/`PANEDTA_PATH` 环境变量显式指定安装目录
