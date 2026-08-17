# PGGB 泛基因组图构建 | PGGB Pangenome Graph Builder

一句话理解：把多个基因组的 FASTA 序列「揉成一张图」——相同的序列只存一份、有差异的地方开叉，得到一张泛基因组变异图，供后续比较、可视化与变异调用。

## 功能概述 | Overview { #overview }

- 用业界标准流程 wfmash（比对）+ seqwish（建图）+ smoothxg（平滑/归一化）一步构建泛基因组变异图
- 输入一个多序列 FASTA（所有基因组合并到一个文件），输出 GFA / OG / GBZ / PAF 等多种格式
- 可选输出 VCF（需 `--vcf-spec` 指定参考）与多序列比对 MAF，可选生成统计（`--stats`）
- 自动为输入 FASTA 建索引（samtools faidx）；支持断点续传（`--resume`）与输出压缩（`--compress`）

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools pggb -i genomes.fa -o output/
```

最小输入：一个包含多基因组序列的 FASTA 文件。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语 | 通俗理解 |
|------|----------|
| 泛基因组变异图 | 用「图」表示多条基因组：节点是序列片段，边是连接关系，差异处形成分叉路径 |
| wfmash | 负责「找相同」，把不同基因组间相似的片段对齐定位 |
| seqwish | 负责「搭骨架」，根据比对把片段连成图结构 |
| smoothxg | 负责「打磨」，给图去重、归一化、平滑，让图更规整 |
| 单倍型（haplotype） | 一条完整基因组的「一份拷贝」；二倍体有两个单倍型 |
| GFA / GBZ / OG | 图的几种存储格式：GFA 是文本图格式，GBZ 是压缩图格式，OG 是 odgi 工具的图格式 |
| PAF | 记录「哪段和哪段相似」的比对结果表格 |
| 比对一致度（map-pct-id） | 认为两段「相同」的最低相似百分比，越高要求越严 |

## 输入 | Input { #input }

### FASTA 文件

- 一个多序列 FASTA，把要纳入泛基因组的全部基因组序列合并写入同一个文件（不同基因组用不同序列名区分）
- 程序会自动用 samtools faidx 建 `.fai` 索引（已存在则跳过）
- 可选：用外部 PAF 跳过 wfmash 比对（直调参数 `--input-paf`，click 包装器未暴露）

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** 告诉程序「读哪个 FASTA、结果放哪」。没有别的讲究。

### 运行资源 | Runtime

**通俗理解|In plain words:** 线程数管并行度，PGGB 是大计算，`-t` 调大能明显提速（默认 24），但也要看机器有多少核。`--conda-env` 是 PGGB 所在 conda 环境名，只在你的环境不叫 `pggb_v.0.7.4` 时才需要改，一般不用动。

### 比对参数（wfmash）| Mapping parameters

**通俗理解|In plain words:** `-s` 决定把序列切成多长的小段去比对，段越短越「细」、越慢，默认 5000 是通用折中。`-p` 是最低相似度，越高要求越严（漏配），越低越宽松（错配），默认 90 适合亲缘较近的物种。`-n` 是单倍型数，默认 0 让程序从序列名前缀自动推断，只有推断错了才需要手动给。

### 输出选项 | Output options

**通俗理解|In plain words:** `--vcf-spec` 指定以哪条参考序列为标准导出 VCF（不填就不出 VCF）。`--stats` 额外生成统计信息，会多花一点时间。`--resume` 是断点续传，中断后重跑同一个输出目录可以接着算。`--compress` 把大文件压缩以省磁盘。

> 更多精细参数（block 长度、mash kmer、稀疏映射、POA 等）未在 `biopytools pggb` 暴露，需直接调用 `python -m biopytools.pggb` 使用，一般用户用不到。

## 分析流程 | Pipeline { #pipeline }

```text
多序列 FASTA
  -> 检查依赖（pggb 可用性，VCF 输出时另查 vg/bcftools）
  -> 确保 .fai 索引存在（samtools faidx）
  -> wfmash 全对全比对（或用外部 PAF 跳过）
  -> seqwish 建图
  -> smoothxg 平滑 / 归一化
  -> 输出 GFA / OG / GBZ / PAF（可选 VCF / MAF / 统计）
```

## 输出 | Output { #output }

```text
output_dir/
├── *.gfa / *.gfa.gz    # 变异图（GFA 文本格式，可压缩）
├── *.og / *.og.gz      # odgi 图格式
├── *.gbz               # GBZ 压缩图格式
├── *.paf / *.paf.gz    # wfmash 比对结果
├── *.vcf.gz / *.vcf    # 变异（仅指定 --vcf-spec 时）
├── *.maf / *.maf.gz    # 多序列比对（可选）
├── *.fai               # 输入 FASTA 索引
├── *.log / *.yml       # 流程日志与参数记录
└── 99_logs/
    └── pggb.log        # biopytools 运行日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

- **GFA / GBZ**：图本体，是后续分析（比对、可视化、变异调用）的入口；GFA 是文本可读，GBZ 更省空间
- **OG**：odgi 工具的图，可用 odgi 直接做统计和可视化
- **PAF**：原始比对表，看不同基因组片段间的对应关系，排查问题时很有用
- **VCF**（可选）：以 `--vcf-spec` 指定的参考序列为基准导出的变异，可与传统变异分析衔接
- **MAF**（可选）：多序列比对，适合做进化分析
- 运行日志 `99_logs/pggb.log` 末尾会按类别列出所有输出文件及其大小

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 亲缘近的物种可把 `-p` 提到 95 减少误配；亲缘远/结构差异大降到 80 左右避免漏配
- 基因组很大、内存吃紧时，把 `-s` 调小（如 3000）可降低单次比对内存，但整体更慢
- 已跑出 PAF、只想换建图参数时，可用直调参数 `--input-paf` 跳过 wfmash 省时间
- 中断重跑务必加 `--resume`；磁盘紧张加 `--compress`
- `-n`、`--conda-env` 等一般不用动

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入FASTA文件｜Input FASTA file |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-t, --threads` | `24` | int | 线程数｜Number of threads (default: 24) |
| `--conda-env` | `pggb_v.0.7.4` |  | conda环境名｜Conda environment name (default: pggb_v.0.7.4) |
| `-s, --segment-length` | `5000` | int | 比对分段长度｜Segment length (default: 5000) |
| `-p, --map-pct-id` | `90` | int | 比对一致度｜Map percent identity (default: 90) |
| `-n, --n-haplotypes` | `0` | int | 单倍型数(0=自动)｜N haplotypes (0=auto) |
| `--vcf-spec` | `` |  | VCF输出参考规范｜VCF output reference spec |
| `--resume` | — |  | 断点续传｜Resume from existing outputs |
| `--compress` | — |  | 压缩输出｜Compress output files |
| `--stats` | — |  | 生成统计信息｜Generate statistics |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | [FILE] 输入FASTA文件｜Input FASTA file |
| `-o, --output` | 必填 |  | [DIR] 输出目录｜Output directory |
| `-t, --threads` | `24` | int | [INT] 线程数 (default: 24) |
| `--conda-env` | `pggb_v.0.7.4` |  | [STR] conda环境名 (default: pggb_v.0.7.4) |
| `-s, --segment-length` | `5000` | int | [INT] 比对分段长度 (default: 5000) |
| `-l, --block-length` | `0` | int | [INT] 最小block长度 (default: 0, auto=5*segment-length) |
| `-p, --map-pct-id` | `90` | int | [INT] 比对一致度 (default: 90) |
| `-c, --n-mappings` | `1` | int | [INT] 每segment的mapping数 (default: 1) |
| `-K, --mash-kmer` | `19` | int | [INT] mash kmer大小 (default: 19) |
| `--no-splits` | — | store_true | 禁用序列拆分｜Disable sequence splitting |
| `--sparse-map` | `` |  | [STR] 稀疏映射比例｜Sparse mapping fraction |
| `--input-paf` | `` |  | [FILE] 外部PAF文件(跳过wfmash)｜External PAF file |
| `-n, --n-haplotypes` | `0` | int | [INT] 单倍型数 (default: 0, auto-detect) |
| `--skip-normalization` | — | store_true | 跳过图归一化｜Skip graph normalization |
| `--vcf-spec` | `` |  | [STR] VCF输出参考规范｜VCF output reference spec |
| `--stats` | — | store_true | 生成统计信息｜Generate statistics |
| `--resume` | — | store_true | 断点续传｜Resume from existing outputs |
| `--compress` | — | store_true | 压缩输出｜Compress output files |
| `--keep-temp` | — | store_true | 保留临时文件｜Keep intermediate files |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- PGGB（conda 环境 `pggb_v.0.7.4`，通过 `--conda-env` 可改）
- samtools（用于 `faidx` 建索引）
- vg + bcftools（仅当指定 `--vcf-spec` 输出 VCF 时额外需要）

## 常见问题 | FAQ { #faq }

**Q1：跑得很慢正常吗？**
正常。PGGB 是全对全比对 + 建图的大计算，程序不设超时，几小时到几十小时都可能，建议用多线程并放后台/提交作业运行。

**Q2：断点续传怎么用？**
重跑时给 `--resume`（对应 pggb 的 `-r`），已完成的步骤会被复用。没给 `--resume` 则从头跑。

**Q3：为什么输入 FASTA 旁边多了个 .fai 文件？**
程序自动用 samtools faidx 建的索引，属正常产物，可放心保留。

**Q4：想输出 VCF 但不知道 --vcf-spec 填什么？**
`--vcf-spec` 填你想作为基准的那条参考序列名（需与输入 FASTA 中某条序列名一致），程序据此导出变异。

**Q5：click 命令里没有的 PGGB 参数怎么用？**
`biopytools pggb` 只暴露常用参数；其余参数需直接调用 `python -m biopytools.pggb`（argparse 入口）传入。
