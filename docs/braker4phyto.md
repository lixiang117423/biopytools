# braker4phyto - 疫霉菌基因组注释(默认不屏蔽重复) | Phytophthora annotation (no repeat masking by default)

一句话理解：**和 annorefine 完全同一条流程，唯一区别是「默认不屏蔽重复序列」——专门给疫霉这类卵菌用，因为它们的效应子基因常常藏在重复区里，屏蔽了反而会丢掉**。
给基因组 + 近缘蛋白(可加 RNA-seq)，一条命令跑完 BRAKER 注释 + 同源查漏补缺，输出整合后的 GFF3。

## 功能概述 | Overview { #overview }

- 与 annorefine 共用同一套端到端流程(BRAKER 注释 + 同源查漏补缺 + 合并拆分 + 质控)
- **唯一差异**：重复序列屏蔽默认关闭(annorefine 默认开启)——疫霉效应子多位于 repeat 区，不 mask 才能保住
- 同样支持小蛋白回收通道、三重质控、断点续传等 annorefine 全部能力
- 适合对象：疫霉属、霜霉属等卵菌，以及一切「重复区里藏真基因」的物种

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools braker4phyto -g genome.fa -s psojae -p prot.fa -o out/
```

最小输入与 annorefine 相同：未屏蔽基因组 + 物种名 + 近缘蛋白。有 RNA-seq 时加 `--rnaseq-dirs` 效果更好。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 卵菌 | 一类长得像真菌但其实是不同进化分支的病原微生物(疫霉、霜霉等) |
| 疫霉 | 卵菌的模式代表，能引起马铃薯晚疫病等重大植物病害 |
| 效应子 | 病原菌分泌出来攻击寄主的小蛋白，常藏在重复区、表达量低、基因短 |
| 重复序列屏蔽 | 把转座子等高重复区域盖住，避免干扰基因预测 |
| 查漏补缺 | 用近缘蛋白把预测软件漏掉的基因补回来 |
| GFF3 | 基因组注释标准格式，记录基因/外显子坐标结构 |

## 输入 | Input { #input }

与 annorefine 完全相同：

- **基因组(`-g`)**：未屏蔽(mask)的原始基因组 FASTA
- **物种名(`-s`)**：输出命名用
- **近缘蛋白(`-p`)**：FASTA 文件或目录
- **RNA-seq(`--rnaseq-dirs`，可选)**：逗号分隔的二代数据目录
- **三代转录本(`--isoseq`，可选)**：文件或目录

## 参数说明 | Parameters { #parameters }

### 与 annorefine 的差异 | Difference from annorefine

**通俗理解|In plain words:** 参数完全一样，只有 `--skip-repeat` 的默认值不同：本命令**默认跳过重复屏蔽**(`--skip-repeat` 默认开)。如果你反而想做屏蔽，用 `--no-skip-repeat` 显式开启——此时行为和 annorefine 默认一致。

### 其余参数 | Other parameters

**通俗理解|In plain words:** 真菌模式(`--fungus`)默认开(疫霉适用)、证据与通用参数、查漏补缺判据、生物学质控、小蛋白回收通道等全部继承 annorefine，含义与默认值一致。详见 annorefine 文档的「参数说明」。做效应子/小蛋白时，建议开启 `--recover-small-proteins`，必要时加 `--no-small-exclude-te`(效应子常在 TE 区)。

## 分析流程 | Pipeline { #pipeline }

与 annorefine 相同，唯一区别在第 1 步重复屏蔽被跳过(默认)或按 `--no-skip-repeat` 执行：

```text
阶段1: BRAKER 注释(重复屏蔽默认跳过, 其余同 annorefine)
  1) [默认跳过] 重复序列屏蔽
  2) 转录组比对(三代 minimap2 / 二代 HISAT2)
  3) BRAKER3 基因结构预测 -> braker.gtf + braker.gff3
    |
    v
阶段2: 同源查漏补缺(miniprot 扫描 + 漏检/合并分析 + 三重质控 + 合并)
  -> 最终整合 GFF3(04_merged/<prefix>.merged.gff3)
```

## 输出 | Output { #output }

输出目录结构与 annorefine 完全一致：

```text
out/
├── 02_long_reads/                  # 三代转录本比对(给了 isoseq 才有)
├── 03_short_reads/                 # 二代 RNA-seq 比对(给了 rnaseq-dirs 才有)
├── 04_braker_annotation/           # BRAKER 预测: braker.gtf / braker.gff3 / braker.aa
├── logs/                           # BRAKER 阶段日志
└── 05_gap_filling/                 # 查漏补缺核心产物
    ├── 01_evidence_scan/<prefix>.miniprot.gff3
    ├── 02_gap_analysis/<prefix>.gap_report.tsv
    ├── 03_gap_filled/<prefix>.gap_filled.gff3
    ├── 04_merged/<prefix>.merged.gff3   # 最终整合结果
    └── 99_logs/annorefine.log
```

注意：因为默认跳过重复屏蔽，**没有 `01_repeat_masking/` 目录**(除非用 `--no-skip-repeat` 显式开启)。

## 结果解读 | Interpreting Results { #results }

与 annorefine 完全一致：

- **最终结果**：`05_gap_filling/04_merged/<prefix>.merged.gff3`，基因 ID 前缀 `<prefix>_gap_N` / `<prefix>_small_gap_N`
- **证据报告**：`05_gap_filling/02_gap_analysis/<prefix>.gap_report.tsv`，逐基因列出蛋白相似度、表达深度、TE 重叠等证据
- 重点看 `te_overlap_pct` 列：本命令不屏蔽重复，意味着候选基因更可能落在 TE 区，**靠表达证据和完整 ORF 把关**，而不是靠屏蔽排除

## 参数选择建议 | Parameter Guidance { #guidance }

- **疫霉/卵菌标准注释**：直接最小四参数 + `--rnaseq-dirs`，其余默认
- **专门做效应子/小蛋白**：加 `--recover-small-proteins` + `--no-small-exclude-te`，必要时 `--small-min-expression-depth 0.1`(效应子低表达)
- **想恢复屏蔽**：加 `--no-skip-repeat`(等价于用 annorefine 默认行为)
- **换参数重跑**：先删 `05_gap_filling` 下对应产物，否则断点续传复用旧结果

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --genome` | 必填 |  | 未mask原始基因组｜Unmasked genome |
| `-s, --species` | 必填 |  | 物种名｜Species name |
| `-p, --prot-seq` | 必填 |  | 近缘蛋白(文件或目录)｜Protein file/dir |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output dir |
| `--rnaseq-dirs` | — |  | 二代RNA-seq目录(逗号分隔)｜RNA-seq dirs |
| `--isoseq` | — |  | 三代转录本(文件或目录)｜Iso-seq file/dir |
| `-t, --threads` | `12` | int | 线程数｜Threads |
| `--fungus/--no-fungus` | `True` |  | 真菌模式(疫霉适用)｜Fungus mode |
| `--no-singularity` | — |  | 不用Singularity｜No singularity |
| `--skip-repeat/--no-skip-repeat` | `True` |  | 跳过repeat屏蔽(默认开, 效应子在repeat区; --no-skip-repeat开启屏蔽)｜Skip repeat masking (default on; --no-skip-repeat to enable) |
| `--skip-repeat-filter` | — |  | 跳过repeat库过滤｜Skip repeat filter |
| `--skip-rescue/--no-skip-rescue` | `True` |  | 证据还原(默认关)｜Rescue (default off) |
| `--split-min-copy-coverage` | `80` | float | 保守合并判据:完整拷贝覆盖率%｜Split copy coverage |
| `--no-split` | — |  | 关闭合并拆分｜Disable split |
| `--repeat-out` | — |  | RepeatMasker .out(filling真TE排除)｜RepeatMasker out |
| `--exclude-te-gap` | — |  | 质控排除TE区gap(默认不排)｜exclude TE-overlap gaps |
| `--no-real-orf` | — |  | 关闭真实完整ORF检查(默认开)｜disable real-ORF check (default on) |
| `--no-coord-zero-overlap` | — |  | 关闭gap坐标零重叠(默认开)｜disable coord-zero-overlap (default on) |
| `--no-unique-reads` | — |  | 关闭唯一比对过滤(默认开)｜disable unique-read filter (default on) |
| `--min-unique-mapq` | `20` | int | 唯一比对MAPQ兜底阈值｜unique MAPQ fallback |
| `--min-expression-depth` | `1.0` | float | 唯一reads平均深度下限｜min unique-read depth |
| `--min-coverage-breadth` | `50.0` | float | CDS覆盖广度%下限｜min coverage breadth |
| `--no-gap-fill` | — |  | 关闭纯漏检填补(只保留合并拆分)｜disable pure gap-fill (split only) |
| `--recover-small-proteins` | — |  | 开启小蛋白回收通道(默认关, 通用)｜enable small-protein lane (default off) |
| `--small-max-cds-len` | `450` | int | 小蛋白CDS上限bp｜small max CDS len |
| `--small-min-identity` | `50.0` | float | 小蛋白放宽identity%(有表达时)｜small min identity |
| `--small-min-coverage` | `50.0` | float | 小蛋白放宽coverage%(有表达时)｜small min coverage |
| `--small-min-expression-depth` | `1.0` | float | 小蛋白表达深度下限(effector低表达可调低如0.1)｜small min expression depth |
| `--small-min-coverage-breadth` | `60.0` | float | 小蛋白CDS覆盖广度%下限｜small min coverage breadth |
| `--no-small-exclude-te` | — |  | 关闭小蛋白TE区排除(effector常在TE区可关)｜disable small-protein TE exclusion |
| `--small-strong-homology-identity` | `95.0` | float | 强同源直通identity%阈值(≥此值绕过TE/表达过滤)｜strong-homology bypass identity |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

与 annorefine 完全一致：

- **BRAKER3**(GeneMark-ETP、AUGUSTUS、TSEBRA)：Singularity 镜像 `~/software/singularity/braker3_devel.sif`
- **RepeatModeler / RepeatMasker**：conda 环境 `repeat`(仅 `--no-skip-repeat` 时用到)
- **minimap2**：`align`；**HISAT2**：`rna`；**miniprot / StringTie**：`annot`；**samtools**
- **Singularity**：`~/miniforge3/envs/singularity_v.3.8.7/bin/singularity`

## 常见问题 | FAQ { #faq }

**Q1：braker4phyto 和 annorefine 到底什么区别？**
同一套代码，唯一区别是重复屏蔽默认值：braker4phyto 默认不屏蔽(保住重复区的效应子)，annorefine 默认屏蔽。其余行为完全一致。

**Q2：为什么疫霉要默认不屏蔽重复？**
疫霉的效应子基因大量位于重复序列区。屏蔽会把重复区盖住，这些真基因就被一起屏蔽掉、预测不出来了。不屏蔽才能保住它们。

**Q3：中断后重跑会从头再来吗？**
不会，断点续传机制与 annorefine 相同(BRAKER 按 braker.gtf、证据扫描按 miniprot.gff3 判断)。

**Q4：我想做屏蔽怎么办？**
加 `--no-skip-repeat` 显式开启屏蔽，或直接用 `biopytools annorefine`。

**Q5：做效应子时要注意什么？**
开启 `--recover-small-proteins` 找回小蛋白；效应子常在 TE 区且低表达，可用 `--no-small-exclude-te` 关闭 TE 排除、`--small-min-expression-depth 0.1` 放宽表达下限。
