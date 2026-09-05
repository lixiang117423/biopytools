# genome2tree:基因组目录直接构建物种进化树 | Alignment-free species tree from genomes

一句话:把一筐基因组(组装 fasta,每物种一个文件)丢进一个目录,一条命令
直接产出考虑"不完全谱系分选"的物种进化树——不需要比对、不需要注释、不需要找基因。

## 功能概述 | Overview

- 输入一个目录,每样本一个 fasta(组装)文件,支持 .gz 压缩
- 底层为 ASTER 包的 waster:k-mer 找 SNP → CASTER coalescent 模型建树
- 输出 Newick 物种树,枝上标注 local bootstrap 支持度(>95 为好)
- 可选:指定外群出有根树;同种多个体合并;在固定拓扑上补算枝长
- 断点续传:已完成的步骤重跑自动跳过

## 快速开始 | Quick Start

```bash
biopytools genome2tree -i genome_dir/ -o results/
```

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗解释<br>In plain words |
|---|---|
| 物种树<br>Species tree | 反映物种分化历史的家谱;单个基因的"家谱"(基因树)可能因杂交、水平转移与物种历史不一致 |
| 不完全谱系分选<br>ILS | 祖先物种的遗传变异在分家时"没分干净",后代的基因树会随机地与真实物种历史矛盾;waster 的模型专门处理这种情况 |
| 支持度<br>Bootstrap support | 每个分支的可信度打分,0-100;>95 表示这个分支非常可靠 |
| 外群<br>Outgroup | 明确最"外围"的物种,用它把树根定下来,树才有方向 |
| K-mer<br>K-mer | 把序列切成定长小片段当"指纹"用;waster 固定用 9bp 的 k-mer 找 SNP,免去全基因组比对 |

## 输入 | Input

一个目录,每个样本一个文件,文件名(去后缀)即样本名:

```
genome_dir/
├── Human.fa          # 组装基因组
├── Chimp.fasta
├── Gorilla.fna.gz    # 支持 .gz(模块自动解压)
└── Orangutan.fas     # 只认 fasta 系后缀,其他文件自动忽略
```

- 认的后缀:`.fa/.fasta/.fna/.fas`(大小写不敏感,可加 `.gz`);fastq、txt 等其他文件一律忽略,不报错
- 同种多个体:准备映射文件(`个体文件名<TAB>物种名` 两列)传 `--samples-map`

## 参数说明 | Parameters

(本节参数表由 scripts/gen_docs_params.py 自动生成)

**通俗理解|In plain words — 基本参数:** `-i` 给目录、`-o` 给输出位置,`-t` 给线程,一般只用这三个。默认 12 线程;waster 并行效率很好,核心多就往高了给。

**通俗理解|In plain words — 进阶参数:** `--root` 指定外群(树会变成有根树,画图好看);`--samples-map` 同种多个体才需要;`--branch-length` 想要枝长(分化时间尺度的长度)才打开,会明显变慢,默认关着即可。

## 输出 | Output

```
results/
├── 00_pipeline_info/software_versions.yml   # 工具版本与参数记录
├── 01_input/input.tsv                       # 最终生效的样本→文件清单
│   ├── samples_map.tsv                      # (有 --samples-map 时)个体→物种映射
│   └── uncompressed/                        # .gz 输入的解压副本
├── 02_waster/waster_species_tree.nw        # 物种树(Newick,含支持度)|species tree (key output)
│   ├── waster_species_tree.nw.snps.fa      # 内部 SNP 矩阵(调试用)
│   └── waster.log
├── 03_branch_length/                        # (--branch-length 时)
│   └── waster_branchlength_species_tree.nw # 带枝长的树|tree with branch lengths
└── 99_logs/genome2tree.log                  # 模块全量日志
```

## 结果解读 | Interpreting Results

- `waster_species_tree.nw`:Newick 格式,枝上方括号数字即 local bootstrap
  - >95:分支很可靠;70-95:可用但需谨慎;<70:该分支关系存疑
- 好坏判据:支持度普遍高=数据量足够;低支持度多半是样本太少(<4 警告)或物种间太远/太近
- `03_branch_length` 的枝长为 coalescent/substitution 单位,不是年代;仅用于相对分化深浅

## 参数选择建议 | Parameter Guidance

- 物种数 4-10:直接跑,不用调参
- 物种数 <4:waster 是四聚体方法,不建议;建议凑样本或改用别的工具
- 基因组只有低覆盖组装:照样能跑,这正是 waster 的强项(reads 需先自行组装为 fasta)
- 想要发表用有根树:加 `--root 外群名`(外群最好选明确外围的近缘种)

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 基因组目录(每样本一个 fasta,可.gz)｜Genome dir (one fasta per sample, .gz ok) |
| `-o, --output-dir` | `./genome2tree_output` |  | 输出目录(默认./genome2tree_output)｜Output directory |
| `-t, --threads` | `12` | int | 线程数(默认12)｜Threads (default 12) |
| `--root` | `` |  | 外群物种名(出有根树)｜Outgroup species name (rooted tree) |
| `--branch-length` | `False` |  | 追加 waster_branchlength 枝长计算｜Also compute branch lengths |
| `--samples-map` | `` |  | 个体→物种映射文件(个体stem<TAB>物种名)｜individual-to-species map |
| `--waster-path` | — |  | waster 路径(默认~/software/ASTER/bin/waster)｜waster binary path |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 基因组目录(每样本一个 fasta,可.gz)｜Genome dir (one fasta per sample, .gz ok) |
| `-o, --output-dir` | `./genome2tree_output` |  | 输出目录(默认./genome2tree_output)｜Output directory |
| `-t, --threads` | `12` | int | 线程数(默认12)｜Threads (default 12) |
| `--root` | `` |  | 外群物种名(出有根树)｜Outgroup species name |
| `--branch-length` | — | store_true | 追加 waster_branchlength 枝长计算｜Also compute branch lengths on the fixed topology |
| `--samples-map` | `` |  | 个体→物种映射文件(个体stem<TAB>物种名)｜individual-to-species map |
| `--waster-path` | — |  | waster 路径(默认~/software/ASTER/bin/waster)｜waster path |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- waster / waster_branchlength(ASTER 包 v1.25.2.5+,`~/software/ASTER/bin/`)
- Python 3 + click + pyyaml
- 内存 ≥64GB(waster 固定 K=9),必须提交计算节点运行

## 常见问题 | FAQ

**Q1: 能不能降内存(比如 16G 服务器)?**
不能。当前 waster 硬编码 K=9,官方文档里的 `-k 7` 是旧版(waster-old)参数,新版已移除;只能申请 ≥64GB 的计算节点。

**Q2: 能传测序 reads(fastq)吗?**
不能。模块只认 fasta 系后缀(`.fa/.fasta/.fna/.fas`),fastq 等文件一律忽略。reads 需先组装成 fasta(低覆盖可用拼接或组装流程)再传入;若两端重叠,官方建议先 BBMerge 合并。

**Q3: 什么时候不该用 genome2tree?**
已有全基因组比对(MAF/多序列比对)→ 用 caster 更准;已有基因家族树 → 用 astral-pro3;物种分歧极远(ANI<80%)→ k-mer 信号衰减,建议传统基因树流程。

**Q4: 树文件用什么看?**
iTOL(在线拖入即可)或 FigTree。

**Q5: 重跑会重复算吗?**
不会。02/03 步产物存在即跳过;换参数重跑请换输出目录或删对应步骤目录。
