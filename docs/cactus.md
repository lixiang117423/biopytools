# Cactus 泛基因组构建 | Cactus Pangenome Construction

一句话理解：**把多个基因组的序列「叠成一张公共地图」（泛基因组图），让不同个体之间的差异——哪里多了、哪里少了、哪里倒了——全部显式地画在图上**，并输出 GFA/HAL/GBZ/ODGI 等多种标准格式，供下游变异分析或可视化继续使用。

## 功能概述 | Overview

- 封装 **Minigraph-Cactus**（`cactus-pangenome`）流程，通过 **Singularity 容器**运行，免去手装一整套 Cactus/Toil 依赖的麻烦
- 输入一个两列的「序列清单」文件（样本名 + FASTA 路径），程序自动校验、自动挂载目录、自动构建
- 默认输出四种常用格式：HAL（永远生成）、GFA（压缩）、GBZ、ODGI；也可按需追加 vg / vcf / xg
- 支持断点续传：所有期望输出文件已存在时自动跳过（见 FAQ）
- 运行成功后默认清理 Toil jobstore 中间目录，节省磁盘（`--no-cleanup` 可保留用于调试）
- 全流程日志中英文对照，并记录完整执行命令，便于论文 Methods 复现

## 快速开始 | Quick Start

```bash
biopytools cactus -s seqfile.txt -o output/ -r ref_genome
```

最小输入：一个两列序列清单 `seqfile.txt`（第一行必须是参考基因组）+ 一个输出目录 + 参考基因组名称。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 泛基因组图 | 不是一条线，而是一张「路网」：所有个体共有的路段只画一条，个体之间不同的路段画成岔路，谁有谁没有一目了然 |
| 参考基因组 | 这张图的「主干道」，其余基因组的差异都以它为基准来对照；本工具要求它排在序列清单第一行 |
| FASTA | 存储 DNA 序列的最通用文本格式，一行 `>` 开头是名字，下面是 A/T/C/G 字母序列 |
| GFA | 描述泛基因组图结构的标准文本格式（路段 + 连接关系），很多下游工具读它 |
| HAL | Cactus 的原生「分层对齐」格式，保留了多基因组两两之间的比对信息，是后续分析的主数据 |
| GBZ / ODGI | 两种压缩/索引的图格式，分别用于高效存储和图上操作（如计算路径、排序） |
| jobstore | Cactus 底层 Toil 的工作「记账本」，记录每一步做了没有、做到哪，跑完即可删除 |
| Singularity 容器 | 一种「软件打包箱」，把 Cactus 及其所有依赖打成一个 .sif 文件，跨机器开箱即用 |

## 输入 | Input

### 序列清单文件（seqfile）

两列格式：**样本名 + FASTA 路径**，制表符或空格分隔。**第一行必须是参考基因组**，且样本名必须与 `--reference` 一致（minigraph-cactus 的硬性要求）。

```text
ref_genome   /path/to/ref.fa
sample2      /path/to/sample2.fa
sample3      /path/to/sample3.fa
```

- 至少需要 2 个基因组（参考 + 查询）
- 支持相对路径（相对于 seqfile 所在目录解析）和绝对路径
- 各基因组文件可以散落在不同目录：程序会自动把它们所在的目录逐个挂载进容器

### 基因组 FASTA 文件

标准 FASTA 格式，程序会逐个校验文件是否存在。单个基因组可以含多条序列（如分染色体/分 scaffold 的多条记录）。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 这三样缺一不可：告诉程序「哪些基因组、放到哪、以谁为基准」。

- `-s, --seqfile`：序列清单文件（两列：样本名 + FASTA 路径），第一行是参考
- `-o, --output`：输出目录
- `-r, --reference`：参考基因组名称，**必须与 seqfile 第一行的样本名完全相同**（拼错会直接报错，见 FAQ）

### 容器与镜像 | Container & image

**通俗理解|In plain words:** Cactus 靠 Singularity 容器运行。这两条路径指向「装箱工具」和「箱子本体」，**只要服务器上按默认位置装好了就一般不用动**，只有你自己另装了别的版本时才需要指。

- `--singularity`：Singularity 可执行文件路径，默认 `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity`
- `--cactus-sif`：Cactus 的 SIF 镜像路径，默认 `~/software/singularity/cactus_v3.1.4.sif`

### 输出与中间产物 | Output & intermediates

**通俗理解|In plain words:** 控制「结果文件叫什么、中间账本要不要留」。`--out-name` 决定所有输出文件的前缀；`--jobstore` 是中间目录名；默认跑完就删掉 jobstore 省空间，只有排查失败原因时才用 `--no-cleanup` 留住它。

- `--out-name`：输出文件前缀，默认 `cactus_output`
- `--jobstore`：Toil jobstore 目录名，默认 `cactus-jobstore`（位于输出目录下）
- `--no-cleanup`：保留 jobstore 不删除（调试用）
- `--formats`：输出格式，可多选，默认 `gfa gbz odgi`；HAL 永远是默认输出（`{out_name}.full.hal`），无需在列表里指定

### 性能参数 | Performance

**通俗理解|In plain words:** 泛基因组构建很吃算力。`--threads` 是同时用多少核，`--max-memory` 是允许每个任务最多用多少内存。**一般按机器规模给：核多就调大 threads，内存大就调大 max-memory**；内存给太小会跑得很慢甚至失败。

- `--threads, -t`：CPU 核心数，默认 12
- `--max-memory, -m`：最大内存（如 `100G`、`50G`），默认 `100G`

### 目录绑定与日志 | Bind & logging

**通俗理解|In plain words:** 容器默认看不到宿主机的所有目录，`--bind` 是额外「开门」给容器看某些目录（程序已自动绑定输入、输出、临时目录，**一般不用手动加**）。日志参数控制输出详细程度。

- `--bind`：额外绑定目录到容器，可多次使用
- `--log-level`：DEBUG / INFO / WARN / ERROR，默认 INFO
- `--quiet`：静默模式

## 分析流程 | Pipeline

```text
序列清单 seqfile.txt
    |
    v
校验环境（Singularity / SIF 镜像 / 版本）
    |
    v
解析 seqfile -> 转成绝对路径版本 + 自动挂载各基因组目录
    |
    v
检查是否已完成（所有期望输出文件都在 -> 跳过）
    |
    v
Singularity exec cactus-pangenome ...（构建泛基因组图）
    |
    v
输出 HAL / GFA / GBZ / ODGI 等，成功后清理 jobstore
```

## 输出 | Output

```text
output/
├── cactus_output.full.hal      # HAL 分层对齐（默认永远生成，主数据）
├── cactus_output.gfa.gz        # GFA 图（压缩，若选了 gfa）
├── cactus_output.gbz           # GBZ 图（若选了 gbz）
├── cactus_output.full.og       # ODGI 图（若选了 odgi，注意后缀是 .full.og 不是 .odgi）
├── cactus_output.log           # 运行日志（中英文对照 + 完整命令）
└── work/                       # 临时工作目录（TMPDIR，供 cactus 写临时文件）
```

文件前缀由 `--out-name` 决定（默认 `cactus_output`）。实际文件名与所选 `--formats` 的对应关系：`hal -> {out_name}.full.hal`、`gfa -> {out_name}.gfa.gz`、`gbz -> {out_name}.gbz`、`odgi -> {out_name}.full.og`。

## 结果解读 | Interpreting Results

- **`.full.hal`**：Cactus 的核心结果，保存了所有输入基因组之间的分层比对关系。文件存在且非空即代表构建成功；后续可用 `halStats`、`hal2maf` 等工具提取特定区间或导出多序列比对。
- **`.gfa.gz`**：泛基因组图结构，解压后是标准 GFA 文本，可被 vg、Bandage、odgi 等工具读取和可视化。若图在 Bandage 里能正常展开成有环有支的网络，说明构建合理。
- **`.gbz` / `.full.og`**：面向高效图存储/图操作的格式，供下游 vg/odgi 生态使用。
- **`cactus_output.log`**：末尾若打印「分析成功完成」并逐个列出 `[OK]` 输出文件，即为成功；若出现 `[MISSING]` 或「Cactus 分析失败」，按日志定位。
- **好坏判据**：输出文件非空、GFA 里同时包含 S（片段）/L（连接）/P（路径）三种记录行，且日志末尾无 ERROR，即为一次正常构建。图里「分叉/气泡」越多，通常代表样本间结构变异越多，属正常现象。

## 参数选择建议 | Parameter Guidance

- **默认 `--formats gfa gbz odgi`**：已覆盖最常见的下游用途；HAL 无需指定、永远生成，除非确定用不到（它也不占额外选择成本）。
- **`--threads` / `--max-memory`**：按集群规模给；建议 threads 8–32、内存不低于 50G，基因组越多/越大越要往上调。
- **`--no-cleanup`**：第一次跑新数据或排查失败时加上，便于查看 jobstore 里的报错；正常生产运行保持默认清理。
- **`--bind`**：仅当基因组文件在程序没自动挂载到的特殊位置（如挂载盘、软链目标）时才需要手动加。
- **`--singularity` / `--cactus-sif`**：默认路径装好了就不动；只在自装新版本时覆盖。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--seqfile, -s` | 必填 |  | 序列文件｜Sequence file (两列格式：样本名 + 路径｜Two columns: sample_name + path) |
| `--output, -o` | 必填 |  | 输出目录｜Output directory |
| `--reference, -r` | 必填 |  | 参考基因组名称｜Reference genome name (必须与seqfile第一个基因组名称匹配｜Must match first genome in seqfile) |
| `--singularity` | `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity` |  | Singularity可执行文件路径｜Singularity executable path |
| `--cactus-sif` | `~/software/singularity/cactus_v3.1.4.sif` |  | Cactus SIF镜像路径｜Cactus SIF image path |
| `--jobstore` | `cactus-jobstore` |  | Toil jobstore目录名称｜Toil jobstore directory name |
| `--out-name` | `cactus_output` |  | 输出文件前缀｜Output file prefix |
| `--no-cleanup` | — |  | 保留jobstore不删除｜Keep jobstore without deleting |
| `--formats` | `['gfa', 'gbz', 'odgi']` | gfa/gbz/odgi/hal/vg/vcf/xg | 输出格式｜Output formats (可多选｜can select multiple) |
| `--threads, -t` | `12` | int | CPU核心数｜Number of CPU cores |
| `--max-memory, -m` | `100G` |  | 最大内存｜Maximum memory |
| `--bind` | — |  | 绑定目录到容器｜Bind directory to container (可多次使用｜can be used multiple times) |
| `--log-level` | `INFO` | DEBUG/INFO/WARN/ERROR | 日志级别｜Log level |
| `--quiet` | — |  | 静默模式｜Quiet mode |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **Singularity**（容器运行器；默认 `~/miniforge3/envs/singularity_v.3.8.7/bin/singularity`）
- **Cactus SIF 镜像**（内置 Minigraph-Cactus v3.1.4 全依赖；默认 `~/software/singularity/cactus_v3.1.4.sif`）
- 无需额外安装 Cactus / Toil / minigraph 等软件，全部随 SIF 镜像提供
- Python 3（运行封装脚本本身）

## 常见问题 | FAQ

**Q1：报「参考基因组不是 seqfile 的第一条样本」？**
minigraph-cactus 硬性要求参考基因组排在 seqfile 第一行，且 `--reference` 名称与第一行样本名完全一致（含大小写、空格）。两者之一对不上就会报错，改其中一处即可。

**Q2：换参数重跑，为什么直接跳过不重算？**
断点续传按「所有期望输出文件是否存在」判断。若改了 `--formats`（例如上次只要 gfa，这次想加 gbz），旧输出目录下已有上次的部分文件，程序会认为「未全部完成」而补跑；但**若想用不同参数彻底重算，应先清空输出目录或改 `--out-name`**，否则会复用旧结果。

**Q3：基因组文件散在不同目录，容器里读不到怎么办？**
程序会自动解析 seqfile 并把每个基因组所在目录挂载进容器，一般无需手动处理。只有在自动挂载失效（如软链指向别处）时才用 `--bind` 手动补挂。

**Q4：jobstore 是什么？要不要留？**
它是 Toil 的工作记账目录，跑完就没用了。程序默认成功后删除以省磁盘；失败时若想排查，加 `--no-cleanup` 保留。

**Q5：为什么会看到「TMPDIR=...」开头的命令？**
这是程序刻意把临时目录指到输出目录下的 `work/`，修复 Cactus 直接写系统 `/tmp` 导致爆盘的已知问题，属正常行为，不用干预。
