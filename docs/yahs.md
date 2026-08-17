# YaHS Hi-C 染色体挂载 | YaHS Hi-C Scaffolding

一句话理解：**用 Hi-C 测序数据（揭示 DNA 在细胞核里「谁挨着谁」）把 contig 挂载、排序、定向成染色体级 scaffolds**。
输入一条 contig 级组装 FASTA 和 Hi-C 双端 reads（R1/R2），输出染色体级 scaffold 序列、AGP、Hi-C 热图和组装质量评估。

## 功能概述 | Overview

- 完整六步流程：建索引 → Hi-C 比对 → YaHS 挂载 → 生成 .hic 热图 → 生成 JBAT 文件 → 质量评估
- **支持断点续传**：已完成的步骤自动跳过，中途失败可安全重跑（`--force-rerun` 强制全重跑）
- 支持单步运行（`-s/--step 1~6`），方便分步调试
- 内置资源自动调整：线程数和内存可按机器自动适配
- 内置比对/BAM 处理中间步骤（bwa mem → 排序 → fixmate → 去重 → 索引）

## 快速开始 | Quick Start

```bash
biopytools yahs -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz
```

最小输入：一条组装 FASTA + Hi-C R1/R2 两个 fastq（.gz 也支持），输出目录默认 `./yahs_output`。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| Hi-C | 一种实验/测序方法，能测出 DNA 在细胞核里哪些片段「离得近」，从而推断序列排列 |
| contig | 一段连续无缺口的序列，是组装的最小单元 |
| scaffold | 多段 contig 用缺口(gap)连成的更大单元 |
| 挂载(scaffolding) | 用 Hi-C 信号把 contig 排成染色体、定方向的过程 |
| 限制性酶切位点 | Hi-C 实验里酶切断 DNA 的识别序列，常见如 GATC |
| MAPQ | 比对质量分，越高越可信 |
| .hic 文件 | Juicer 生态的热图文件，用 Juicebox 打开看 Hi-C 矩阵 |
| JBAT | Juicebox 的自动/手动校正格式，用于人工纠错 |
| AGP | 描述「哪段序列由哪些 contig 拼起来」的文本文件 |
| N50 | 组装连续性的常用指标：从长到短累加，加到总长一半时那条序列的长度 |

## 输入 | Input

### 参考基因组 FASTA

contig 级组装（待挂载对象）。内部会复制/硬链接到输出目录再建 BWA + samtools 索引。

### Hi-C reads

`-1` 和 `-2` 两个文件（Hi-C 双端 reads），支持 .fq/.fastq/.gz 压缩格式。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 三个都要给：组装 FASTA 和 Hi-C 两端 reads。缺一不可。

### 资源参数 | Resources

**通俗理解|In plain words:** `-t/--threads` 是线程数；`--java-ram` 给 Java（juicer_tools）内存，`--sam-ram` 给 samtools 排序内存。**注意：这两个内存参数留默认时会被自动调整**——Java 取机器总内存的 60%，samtools 排序取「总内存一半」和 300G 的较大值（即至少 300G）。机器内存不够时要主动调小，否则可能申请不到内存而失败（见 FAQ）。

### YaHS 核心参数 | YaHS core parameters

**通俗理解|In plain words:** `-e/--enzyme` 是 Hi-C 实验用的酶切识别序列，常见默认 `GATC`，**必须和实验实际用的酶一致**，否则热图和挂载都会错。`--min-len` 过滤掉过短的 contig（默认 10000bp），`--min-mapq` 过滤低质量比对（默认 30），`--no-contig-ec`/`--no-scaffold-ec` 分别跳过 contig/scaffold 错误校正，`--telo-motif` 指定端粒序列用于判定染色体末端。这些**一般用默认即可**，只有确认数据特殊时才动。

### 工具路径 | Tool paths

**通俗理解|In plain words:** 各软件的可执行文件路径，都有默认值（YaHS/juicer 在 hic 环境、bwa/samtools 在 align 环境、juicer_tools.jar 在 software 目录）。**除非你的安装位置不同，否则一般不用改**。

### 执行控制 | Execution control

**通俗理解|In plain words:** `-s/--step` 只跑指定某一步（1 索引 / 2 比对 / 3 挂载 / 4 热图 / 5 JBAT / 6 评估）；`--force-rerun` 无视断点续传强制全重跑；`--keep-temp` 保留中间临时文件（默认清理，省磁盘）。

## 分析流程 | Pipeline

```text
组装 FASTA + Hi-C R1/R2
    │
    ▼
步骤1: 构建索引(bwa index + samtools faidx)         → 01_indexing/
    │
    ▼
步骤2: Hi-C 比对(bwa mem -5SP → SAM转BAM → 按名排序
       → fixmate → 按坐标排序 → markdup → 索引)      → 02_mapping/aligned_sorted_dedup.bam
    │
    ▼
步骤3: YaHS 挂载                                    → 03_scaffolding/yahs_out_scaffolds_final.fa + .agp
    │
    ▼
步骤4: 生成标准 .hic 热图(juicer pre + juicer_tools) → 04_hic_standard/yahs_out_final.hic
    │
    ▼
步骤5: 生成 JBAT 文件(供 Juicebox 人工校正)          → 05_jbat/out_JBAT.hic + .assembly
    │
    ▼
步骤6: 组装质量评估(N50/N90/L50/L90)                 → 06_assessment/assembly_metrics.txt
```

## 输出 | Output

```text
yahs_output/
├── 00_pipeline_info/
│   └── software_versions.yml          # 软件版本信息
├── 01_indexing/                       # 基因组索引(bwa + samtools)
│   ├── {ref}.bwt / .pac / .ann / .amb / .sa
│   └── {ref}.fai
├── 02_mapping/
│   └── aligned_sorted_dedup.bam(+.bai)  # 去重后的比对
├── 03_scaffolding/
│   ├── yahs_out_scaffolds_final.fa      # 最终染色体级 scaffold(核心产物)
│   ├── yahs_out_scaffolds_final.agp     # 挂载关系
│   └── yahs_out.bin                     # YaHS 二进制中间文件
├── 04_hic_standard/
│   └── yahs_out_final.hic               # 标准 Hi-C 热图
├── 05_jbat/
│   ├── out_JBAT.hic                     # JBAT 热图
│   ├── out_JBAT.assembly                # 供 Juicebox 校正
│   └── out_JBAT.liftover.agp
├── 06_assessment/
│   └── assembly_metrics.txt             # N50/N90/L50/L90 统计
└── 99_logs/
    └── yahs_pipeline.log                # 运行日志
```

- `03_scaffolding/yahs_out_scaffolds_final.fa`：最终挂载结果，下游分析的主输入
- `03_scaffolding/yahs_out_scaffolds_final.agp`：每条 scaffold 由哪些 contig 拼成，复核挂载的依据
- `04_hic_standard/yahs_out_final.hic`：用 Juicebox 打开查看 Hi-C 矩阵
- `06_assessment/assembly_metrics.txt`：组装连续性指标，看挂载前后提升多少

## 结果解读 | Interpreting Results

- **scaffold 数应接近染色体数**：挂载后序列条数应明显减少，接近预期染色体条数（可能略多，含未定位 scaffold）
- **N50 应显著提升**：挂载前是 contig N50，挂载后是 scaffold N50，通常提升数倍到数十倍
- **.hic 热图呈「对角线方块」**：好的挂载在 Juicebox 里沿对角线是整齐的方块；出现大片「十字交叉」或噪点，说明存在错误挂载，可用 JBAT 文件人工校正
- **AGP 里 gap 多而大**：说明某些染色体 contig 之间缺口多，组装连续性还有提升空间
- **端粒（若指定 `--telo-motif`）出现在 scaffold 两端**：说明染色体挂载完整、方向正确

## 参数选择建议 | Parameter Guidance

- **标准 Hi-C 默认酶**：`-e GATC`（DpnII/MboI 类）
- **酶不是 GATC**：按实验记录改成实际识别序列（如 HindIII 用 `AAGCTT`）
- **内存有限**：主动设 `--java-ram 32G --sam-ram 100G`，避免默认自动调整申请过大内存
- **只补跑挂载**：`-s 3`
- **换参数重跑**：加 `--force-rerun`，或删掉对应步骤的 checkpoint 文件

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-r, --ref` | 必填 |  | 参考基因组FASTA｜Reference genome FASTA |
| `-1, --hic-r1` | 必填 |  | Hi-C R1文件｜Hi-C R1 file |
| `-2, --hic-r2` | 必填 |  | Hi-C R2文件｜Hi-C R2 file |
| `-o, --output-dir` | `./yahs_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `--java-ram` | `300G` |  | Java内存｜Java memory |
| `--sam-ram` | `300G` |  | Samtools排序内存｜Samtools sort memory |
| `-e, --enzyme` | `GATC` |  | 限制性酶切位点｜Restriction enzyme sequence |
| `--min-len` | `10000` | int | 最小contig长度｜Minimum contig length |
| `--min-mapq` | `30` | int | 最小MAPQ值｜Minimum MAPQ value |
| `--no-contig-ec` | — |  | 跳过contig错误校正｜Skip contig error correction |
| `--no-scaffold-ec` | — |  | 跳过scaffold错误校正｜Skip scaffold error correction |
| `--resolutions` | — |  | 分辨率列表(逗号分隔)｜Resolution list (comma-separated) |
| `--rounds` | `1` | int | 每分辨率运行轮数｜Rounds per resolution |
| `--telo-motif` | — |  | 端粒序列模体｜Telomeric sequence motif |
| `--yahs-bin` | `~/miniforge3/envs/hic/bin/yahs` |  | YaHS可执行文件｜YaHS executable |
| `--juicer-bin` | `~/miniforge3/envs/hic/bin/juicer` |  | juicer可执行文件｜juicer executable |
| `--juicer-jar` | `~/software/juicer/scripts/juicer_tools.jar` |  | juicer_tools.jar文件｜juicer_tools.jar file |
| `--bwa-bin` | `~/miniforge3/envs/align/bin/bwa` |  | BWA可执行文件｜BWA executable |
| `--samtools-bin` | `~/miniforge3/envs/align/bin/samtools` |  | samtools可执行文件｜samtools executable |
| `--java-cmd` | `java` |  | Java可执行文件｜Java executable |
| `-s, --step` | — | 1/2/3/4/5/6 | 运行指定步骤｜Run specified step |
| `--force-rerun` | — |  | 强制重新运行｜Force rerun all steps |
| `--keep-temp` | — |  | 保留临时文件｜Keep temporary files |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-r, --ref` | 必填 |  | 参考基因组FASTA文件｜Reference genome FASTA file |
| `-1, --hic-r1` | 必填 |  | Hi-C R1测序文件｜Hi-C R1 sequencing file |
| `-2, --hic-r2` | 必填 |  | Hi-C R2测序文件｜Hi-C R2 sequencing file |
| `-o, --output-dir` | `./yahs_output` |  | 输出目录｜Output directory (default: ./yahs_output) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `--java-ram` | `32G` |  | Java内存｜Java memory (default: 32G) |
| `--sam-ram` | `4G` |  | Samtools排序内存｜Samtools sort memory (default: 4G) |
| `-e, --enzyme` | `GATC` |  | 限制性酶切位点｜Restriction enzyme sequence (default: GATC) |
| `--min-len` | `10000` | int | 最小contig长度｜Minimum contig length (default: 10000) |
| `--min-mapq` | `30` | int | 最小MAPQ值｜Minimum MAPQ value (default: 30) |
| `--no-contig-ec` | — | store_true | 跳过contig错误校正｜Skip contig error correction |
| `--no-scaffold-ec` | — | store_true | 跳过scaffold错误校正｜Skip scaffold error correction |
| `--resolutions` | — |  | 分辨率列表(逗号分隔)｜Resolution list (comma-separated) |
| `--rounds` | `1` | int | 每分辨率运行轮数｜Rounds per resolution (default: 1) |
| `--telo-motif` | — |  | 端粒序列模体｜Telomeric sequence motif |
| `--yahs-bin` | `~/miniforge3/envs/hic/bin/yahs` |  | YaHS可执行文件路径｜YaHS executable path (default: ~/miniforge3/envs/hic/bin/yahs) |
| `--juicer-bin` | `~/miniforge3/envs/hic/bin/juicer` |  | juicer可执行文件路径｜juicer executable path (default: ~/miniforge3/envs/hic/bin/juicer) |
| `--juicer-jar` | `~/software/juicer/scripts/juicer_tools.jar` |  | juicer_tools.jar文件路径｜juicer_tools.jar file path (default: ~/software/juicer/scripts/juicer_tools.jar) |
| `--bwa-bin` | `~/miniforge3/envs/align/bin/bwa` |  | BWA可执行文件路径｜BWA executable path (default: ~/miniforge3/envs/align/bin/bwa) |
| `--samtools-bin` | `~/miniforge3/envs/align/bin/samtools` |  | samtools可执行文件路径｜samtools executable path (default: ~/miniforge3/envs/align/bin/samtools) |
| `--java-cmd` | `java` |  | Java可执行文件路径｜Java executable path (default: java) |
| `-s, --step` | — | 1/2/3/4/5/6 | 运行指定步骤(1-6)｜Run specified step (1-6) |
| `--force-rerun` | — | store_true | 强制重新运行所有步骤｜Force rerun all steps |
| `--keep-temp` | — | store_true | 保留临时文件｜Keep temporary files |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- conda 环境 `hic`：YaHS（`~/miniforge3/envs/hic/bin/yahs`）、juicer（`~/miniforge3/envs/hic/bin/juicer`）
- conda 环境 `align`：BWA（`~/miniforge3/envs/align/bin/bwa`）、samtools（`~/miniforge3/envs/align/bin/samtools`）
- `juicer_tools.jar`（默认 `~/software/juicer/scripts/juicer_tools.jar`），用于生成 .hic 文件
- Java（`--java-cmd`，默认 `java`）

## 常见问题 | FAQ

**Q1：有断点续传吗？**
有。每个步骤以关键输出文件作为 checkpoint，存在即跳过；中途失败重跑会自动从断点继续。想强制全部重跑加 `--force-rerun`。

**Q2：换参数重跑结果没变？**
因为断点续传按 checkpoint 文件存在性判断，不会因参数变化自动重跑。换参数请加 `--force-rerun`，或先删除对应步骤的产物。

**Q3：为什么内存参数「默认 300G」这么高？**
CLI 帮助里的 `300G` 其实是个哨兵值。留默认时内部会**自动检测机器内存**：Java 取总内存 60%，samtools 排序取 max(总内存一半, 300G)。机器内存小时建议显式设小，如 `--java-ram 32G --sam-ram 100G`，否则可能申请不到内存而失败。

**Q4：磁盘占用为什么很大？**
步骤 2 的 BWA 先输出**文本 SAM**（巨大）再转 BAM，且排序会产生临时文件。请预留足够磁盘；`--keep-temp` 会保留中间文件，默认会清理。

**Q5：热图/JBAT 步骤没生成 .hic 文件？**
步骤 4/5 需要 `juicer_tools.jar` 和 Java。若 jar 缺失会跳过并提示；确认 `--juicer-jar` 路径正确、`java` 可用。

**Q6：挂载结果和酶对不上有关系吗？**
关系很大。`-e/--enzyme` 必须和 Hi-C 建库实际用的酶一致，否则 Hi-C 信号位置错乱，挂载和热图都会错。
