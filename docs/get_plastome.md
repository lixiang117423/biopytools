# 叶绿体基因组组装 | Plastome Assembly

一句话理解：**把一份二代全基因组测序（WGS）数据丢进去，自动从中「捞」出叶绿体基因组序列**。输入一个存放 clean reads 的目录，输出每个样品的叶绿体基因组 FASTA，省去手工从海量核基因组 reads 里翻找叶绿体的麻烦。

## 功能概述 | Overview

- 基于 GetOrganelle，从全基因组测序（WGS）reads 中自动组装细胞器基因组（默认叶绿体）
- 默认批量模式：自动扫描输入目录里所有成对 reads 样品，逐个组装
- 支持 7 种细胞器类型：被子植物叶绿体（默认）、线粒体、动物线粒体、真菌线粒体等
- 自动整理结果：把 GetOrganelle 主产物重命名为 `{sample}.plastome.fasta`，序列每行 60 碱基
- 无断点续传：重跑会重新执行组装（见 FAQ）

## 快速开始 | Quick Start

```bash
biopytools get-plastome -i fastq_folder -o plastome_output
```

最小输入：一个目录，里面按 `{sample}_1.clean.fq.gz` / `{sample}_2.clean.fq.gz` 命名存放成对 clean reads。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 叶绿体基因组 | 叶绿体自带的一套小 DNA「小册子」，独立于核基因组，被子植物通常约 120–170 kb |
| 全基因组测序(WGS) | 把细胞里所有 DNA（核 + 叶绿体 + 线粒体）混在一起测序，像一锅大杂烩 |
| GetOrganelle | 专门从这锅「大杂烩」里只挑细胞器序列（叶绿体/线粒体）的工具 |
| k-mer | 把序列按固定长度切成小碎片，像拼图的小块；程序用不同大小的块反复尝试拼装 |
| contig | 拼装出来的连续序列片段 |
| 组装图(graph) | 记录碎片之间如何互相连接的「拼装说明书」，GFA 格式 |

## 输入 | Input

### reads 目录

一个目录，包含一个或多个样品的成对二代 clean reads，默认按以下后缀识别：

```text
sampleA_1.clean.fq.gz     # R1
sampleA_2.clean.fq.gz     # R2
sampleB_1.clean.fq.gz
sampleB_2.clean.fq.gz
...
```

- 默认 R1 后缀 `_1.clean.fq.gz`、R2 后缀 `_2.clean.fq.gz`（如需改后缀，用模块直调参数 `--read1-suffix`/`--read2-suffix`，见参数速查）
- 支持 gzip 压缩（.fq.gz）；也支持单端 reads（仅模块直调）

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** `-i` 指定「测序数据放在哪个文件夹」，程序会自动找出里面所有样品；`-o` 指定结果写到哪。

### 细胞器类型 | Organelle type

**通俗理解|In plain words:** 告诉程序你要找的是哪种细胞器的「小册子」。默认 `embplant_pt` 是被子植物叶绿体，绝大多数植物项目直接用默认即可；动物换 `animal_mt`，真菌换 `fungus_mt`，非编码区换 `*_nr`。

### 组装参数 | Assembly parameters

**通俗理解|In plain words:** `--max-rounds` 是「反复扩展拼装的轮数上限」，`--kmer-list` 是「用哪些大小的碎片去拼」。这两个是 GetOrganelle 的核心参数，**默认值对绝大多数被子植物叶绿体已经调好，一般不用动**。数据特别难拼（如高重复区）时才可能需要加大轮数或调整 k-mer 列表。

### 软件路径与日志 | Software path & logging

**通俗理解|In plain words:** `--getorganelle-path` 指向 GetOrganelle 脚本，一般不用动；`-v`/`--log-file` 用于排查问题时输出更多信息。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先扫描目录找出所有样品，再对每个样品分别跑 GetOrganelle，最后把主要产物整理成带样品名的干净 FASTA。

```text
输入目录(fastq)
    │
    ▼
扫描样品: 按 _1/_2.clean.fq.gz 后缀识别成对 reads
    │
    ▼
对每个样品: GetOrganelle 组装(-F 细胞器类型 -R 轮数 -k kmer列表 -t 线程)
    │
    ▼
整理结果: graph1.1.path_sequence.fasta → {sample}.plastome.fasta(每行60碱基)
```

## 输出 | Output

```text
plastome_output/
└── {sample}/
    ├── {sample}.plastome.fasta        # 最终叶绿体序列(重命名+格式化)
    ├── {sample}.graph.gfa             # 组装图(重命名自 selected_graph.gfa)
    ├── *.path_sequence.fasta          # GetOrganelle 原始路径序列(如 graph1.1/others)
    ├── *.selected_graph.gfa           # GetOrganelle 选出的组装图
    └── *.log.txt                      # GetOrganelle 运行日志
```

- `{sample}.plastome.fasta`：最关心的最终序列，优先选择 `graph1.1` 路径（通常是最优解）重命名而来
- `{sample}.graph.gfa`：对应的组装图，可导入 Bandage 等工具可视化

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 打开 `{sample}.plastome.fasta`，看它是不是「一本完整、长度合理的小册子」。

- **长度**：被子植物叶绿体通常约 120–170 kb。远小于此（如几十 kb）往往是没拼全；远大于此可能混入了核基因组或线粒体污染
- **序列数**：正常的环状叶绿体应拼成 1 条（或 1–2 条含重复区的情况）；如果拆成很多条 contig，说明数据或参数不够好
- **完整性**：可进一步用 BUSCO（embryophyta_odb10）或比对近缘物种的叶绿体参考来验证
- 日志里若找不到 `*.path_sequence.fasta`，说明组装失败，查看 `*.log.txt` 排查

## 参数选择建议 | Parameter Guidance

- `--organelle-type`：植物默认 `embplant_pt`；动物线粒体 `animal_mt`；叶绿体和线粒体各跑一次即可都拿到
- `--max-rounds` / `--kmer-list`：**默认即可**，只有拼装结果差时再尝试调大轮数或改 k-mer 列表
- `--threads`：默认 12，样品多或数据大时可适当调高
- 后缀不符：若文件名不是 `_1.clean.fq.gz` 风格，用模块直调参数 `--read1-suffix`/`--read2-suffix` 指定（见参数速查「模块直调参数」）

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入目录(包含reads文件)｜Input directory containing reads files |
| `--output-dir, -o` | `./plastome_output` | Path | 输出目录｜Output directory |
| `--prefix, -p` | — |  | 输出前缀｜Output prefix |
| `--organelle-type` | `embplant_pt` | embplant_pt/embplant_mt/embplant_nr/other_pt/animal_mt/fungus_mt/fungus_nr | 细胞器类型｜Organelle type |
| `--max-rounds, -R` | `15` | int | 最大扩展轮数｜Maximum extension rounds |
| `--kmer-list, -k` | `21,45,65,85,105` |  | Kmer列表(逗号分隔)｜Kmer list comma-separated |
| `--threads, -t` | `12` | int | 线程数｜Threads |
| `--getorganelle-path` | `~/miniforge3/envs/asm/bin/get_organelle_from_reads.py` |  | GetOrganelle脚本路径｜GetOrganelle script path |
| `--verbose, -v` | — |  | 详细输出模式｜Verbose output mode |
| `--log-file` | — |  | 日志文件路径｜Log file path |
| `--single-mode` | — |  | 单样品模式｜Single sample mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入目录(包含reads文件)｜Input directory containing reads files |
| `-o, --output-dir` | `./plastome_output` |  | 输出目录｜Output directory |
| `-p, --prefix` | — |  | 输出前缀｜Output prefix |
| `--organelle-type` | `embplant_pt` | embplant_pt/embplant_mt/embplant_nr/other_pt/animal_mt/fungus_mt/fungus_nr | Organelle类型｜Organelle type |
| `-R, --max-rounds` | `15` | int | 最大扩展轮数｜Maximum extension rounds |
| `-k, --kmer-list` | `21,45,65,85,105` |  | Kmer列表(逗号分隔)｜Kmer list comma-separated |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--read1-suffix` | `_1.clean.fq.gz` |  | R1文件后缀模式｜R1 file suffix pattern (default: %(default)s) |
| `--read2-suffix` | `_2.clean.fq.gz` |  | R2文件后缀模式｜R2 file suffix pattern (default: %(default)s) |
| `--getorganelle-path` | `~/miniforge3/envs/asm/bin/get_organelle_from_reads.py` |  | GetOrganelle脚本路径｜GetOrganelle script path |
| `-v, --verbose` | — | store_true | 详细输出模式｜Verbose output mode |
| `--log-file` | — |  | 日志文件路径｜Log file path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- GetOrganelle（脚本默认路径 `~/miniforge3/envs/asm/bin/get_organelle_from_reads.py`，conda 环境 `asm`）
- 该脚本运行依赖 GetOrganelle 自身及其配套工具（bowtie2、SPAdes、blast 等，随 GetOrganelle 环境安装）

## 常见问题 | FAQ

**Q1：输入目录里一个样品都没识别到？**
检查文件名是否严格符合 `{sample}_1.clean.fq.gz` / `{sample}_2.clean.fq.gz`。若后缀不同，用模块直调参数 `--read1-suffix` / `--read2-suffix` 指定（这两个参数 click 包装器未暴露）。

**Q2：支持断点续传吗？**
不支持。重跑会重新执行 GetOrganelle 组装。若中途失败，建议检查 `*.log.txt` 后重跑。

**Q3：为什么结果里有多条序列（多个 contig）？**
正常叶绿体应拼成接近 1 条。多条通常是：数据量不足、重复区（如 IR 反向重复）导致的拆分、或混入了核基因组/线粒体。可尝试增加 `--max-rounds` 或提供更高质量的 clean reads。

**Q4：路径含中文会有问题吗？**
程序内部把当前工作目录作为运行目录、并对路径做了标准化处理，以规避中文路径问题；若仍报错，可尝试把数据放到纯英文路径下。

