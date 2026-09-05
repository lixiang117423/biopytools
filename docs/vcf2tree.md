# vcf2tree - VCF 直接建树 | vcf2tree Phylogenetic Tree from VCF

一句话理解：**从 VCF(SNP) 一键建进化树：先转成 IUPAC FASTA 比对，再交给 IQ-TREE(默认,自动选模型+UFBoot)或 FastTree 建树。**

## 功能概述 | Overview

- 两种建树后端：`iqtree`(默认,精度高)和 `fasttree`(速度快)
- IQ-TREE 默认 ModelFinder 自动选模型 + UFBoot 1000 次，并默认开启 ASC 校正
- 内置 VCF -> FASTA(IUPAC 编码)转换，只保留二等位单碱基 SNP
- 两步均支持断点续传(已有结果则跳过)

## 快速开始 | Quick Start

```bash
biopytools vcf2tree -i variants.vcf -o tree_output
```

最小输入：一个 VCF 文件(`-i`)。默认用 IQ-TREE 建树，结果写在 `-o` 目录。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| VCF | 记录「每个样本在每个位点是什么基因型」的标准变异文件 |
| IUPAC 码 | 用单个字母表示「可能是哪几个碱基」的模糊编码(杂合用 R/Y/S 等) |
| 最大似然(ML) | 在众多候选树里挑「最可能产生眼前数据」的那棵 |
| ModelFinder | IQ-TREE 内置的「自动挑最合适进化模型」的模块 |
| UFBoot | IQ-TREE 的「快速 bootstrap」，给每个分叉打支持值 |
| ASC 校正 | 针对「只含变异位点、缺恒定位点」的 SNP 数据的校正，防止分支长度被低估 |
| 恒定/可变位点 | 恒定=所有样本都一样(没信息)；可变=有差异(有信息)。SNP 数据天然只保留可变位点 |

## 输入 | Input

### VCF 文件

支持 `.vcf` 和 `.vcf.gz`。需含 `#CHROM` 样本头行。示例：

```text
##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1  sample2  sample3
Chr1  100  .  A  G  .  PASS  .  GT  0/0  0/1  1/1
```

只有二等位单碱基 SNP 会被保留；多等位、InDel 会被跳过。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 只有 `-i` 是必填(输入 VCF)。`-o/--output-dir` 默认 `./vcf2tree_output`。

### 建树方法 | Tree method

**通俗理解|In plain words:** `-m/--method` 选后端：`iqtree`(默认)精度高、自动选模型；`fasttree` 快、适合超大 SNP 数据集但精度略低。**追求准确用默认 iqtree，数据特别大赶时间才换 fasttree。**

### 运行参数 | Run parameters

**通俗理解|In plain words:** `-t/--threads` 是线程数(默认 12)。`-g/--outgroup` 指定外群样本名给树定方向。`--min-samples-locus` 是「位点最少样本数」(默认 4)——缺失太严重的位点直接丢弃，**一般不用动**。

### FastTree 参数 | FastTree options

**通俗理解|In plain words:** 只在 `--method fasttree` 时生效。`--fasttree-path` 是软件路径(默认 `~/.local/bin/FastTree`)。`--fasttree-params` 可追加额外 FastTree 参数(如 `-gtr -nosupport`)。**注意 FastTree 没有 ASC 校正，SNP 数据的分支长度会被低估**，需要 ASC 时请用 IQ-TREE 后端。

### IQ-TREE 参数 | IQ-TREE options

**通俗理解|In plain words:** 只在 `--method iqtree` 时生效。`--iqtree-bootstrap` 是 UFBoot 重复次数(默认 1000，0 表示跳过 bootstrap)。`--iqtree-model` 手动指定模型(留空则 ModelFinder 自动选)。`--iqtree-path` 是软件路径(默认已配好)。`--no-asc` **关闭 ASC 校正**——ASC 对 SNP 数据很重要，**一般不要关**。

## 分析流程 | Pipeline

```text
VCF 文件
    │
    ▼
步骤1: VCF -> FASTA(IUPAC 编码比对)，只留二等位单碱基 SNP
    │
    ▼
步骤2: 建树
    ├─ IQ-TREE(默认): ModelFinder 自动选模型 + UFBoot + ASC 校正
    └─ FastTree(可选): -nt 核酸模式快速建树(无 ASC)
    │
    ▼
步骤3: 记录软件版本(00_pipeline_info/software_versions.yml)
```

## 输出 | Output

```text
tree_output/
├── 01_vcf2fasta/
│   └── <base>_snps.fa            # SNP 的 IUPAC FASTA 比对
├── 02_tree/
│   ├── <base>_<method>.nwk       # 最终进化树(Newick)
│   └── <base>.treefile 等        # IQ-TREE 中间文件(仅 iqtree 后端)
├── 00_pipeline_info/
│   └── software_versions.yml     # 软件版本与参数记录
└── 99_logs/
    └── vcf2tree.log              # 运行日志
```

其中 `<base>` 是输入 VCF 的 basename，`<method>` 是 `iqtree` 或 `fasttree`。

## 结果解读 | Interpreting Results

### 1. 进化树(<base>_<method>.nwk)

**通俗理解|In plain words:** 最终结果，Newick 格式，用 FigTree、iTOL 打开即见树形图。IQ-TREE 后端的分叉上带 UFBoot 支持值(越高越可信)。

### 2. SNP FASTA 比对(<base>_snps.fa)

**通俗理解|In plain words:** 中间产物，VCF 转成的 IUPAC 比对。杂合位点用模糊码(如 R、Y)表示。可以单独拿去给其它建树软件用。

### 3. 版本信息(software_versions.yml)

**通俗理解|In plain words:** 记录本次用到的建树软件版本和关键参数，写论文 Methods 时直接抄。

## 参数选择建议 | Parameter Guidance

- **默认方案**：`-m iqtree`(默认)，自动选模型 + UFBoot 1000，精度最好
- **超大 SNP 数据集、赶时间**：`-m fasttree`(快，但无 ASC、分支长度低估)
- **不要 bootstrap、只要一棵树**：IQ-TREE 下 `--iqtree-bootstrap 0`
- **已有明确模型**：`--iqtree-model GTR+G`(跳过 ModelFinder，省时)
- **给树定方向**：`-g <外群样本名>`

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `--method, -m` | `iqtree` | fasttree/iqtree | 建树方法｜Tree method (默认: iqtree｜default: iqtree) |
| `--output-dir, -o` | `./vcf2tree_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--outgroup, -g` | `` |  | 外群样本名称｜Outgroup sample name |
| `--min-samples-locus` | `4` | int | 位点最少样本数｜Minimum samples per locus |
| `--fasttree-path` | `~/.local/bin/FastTree` |  | FastTree软件路径｜FastTree software path |
| `--fasttree-params` | `` |  | FastTree额外参数｜Additional FastTree parameters |
| `--iqtree-path` | — |  | IQ-TREE软件路径｜IQ-TREE software path |
| `--iqtree-bootstrap` | `1000` | int | IQ-TREE UFBoot重复次数｜IQ-TREE UFBoot replicates |
| `--iqtree-model` | — |  | IQ-TREE模型(默认ModelFinder)｜IQ-TREE model (default: ModelFinder) |
| `--no-asc` | `False` |  | 关闭SNP数据的ASC校正(默认开启)｜Disable ASC correction for SNP data (on by default) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `--method` | `iqtree` | fasttree/iqtree | 建树方法｜Tree method: fasttree or iqtree (default: iqtree) |
| `-o, --output-dir` | `./vcf2tree_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Number of threads |
| `-g, --outgroup` | `` |  | 外群样本名称｜Outgroup sample name |
| `--min-samples-locus` | `4` | int | 位点最少样本数｜Minimum samples per locus |
| `--fasttree-path` | `~/.local/bin/FastTree` |  | FastTree软件路径｜FastTree software path |
| `--fasttree-params` | `` |  | FastTree额外参数｜Additional FastTree parameters |
| `--iqtree-path` | — |  | IQ-TREE软件路径｜IQ-TREE software path |
| `--iqtree-bootstrap` | `1000` | int | IQ-TREE UFBoot重复次数｜IQ-TREE UFBoot replicates |
| `--iqtree-model` | — |  | IQ-TREE进化模型(默认ModelFinder自动)｜IQ-TREE model (default: ModelFinder) |
| `--no-asc` | — | store_true | 关闭SNP数据的ASC校正(默认开启)｜Disable ASC correction for SNP data (on by default) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

| 软件<br>Software | 说明<br>Description |
|------|------|
| IQ-TREE | 默认后端，默认路径 `~/miniforge3/envs/phylo/bin/iqtree`(conda 环境 phylo)，可用 `--iqtree-path` 覆盖 |
| FastTree | 可选后端，默认路径 `~/.local/bin/FastTree`，可用 `--fasttree-path` 覆盖 |

两者均为独立二进制：在 conda 环境内自动用 `conda run` 包装，否则直接调用。

## 常见问题 | FAQ

### 1. 会断点续传吗？

会。**步骤1(VCF->FASTA)和步骤2(建树)各自独立续传**：对应输出文件已存在且非空时自动跳过。想强制重跑，删除 `01_vcf2fasta/<base>_snps.fa` 或 `02_tree/<base>_<method>.nwk` 即可。

### 2. 为什么用 SNP 数据建树要开 ASC 校正？

SNP 数据通常只含「可变位点」(恒定不变的位置在叫变异时被去掉了)。若不做 ASC 校正，模型会误以为「没有恒定位点」，导致**分支长度被系统性低估**。本工具默认开启 ASC(`+ASC`)，**一般不要用 `--no-asc` 关掉**。ASC 要求比对中不含恒定位点，本工具会在转 FASTA 时自动移除「仅杂合变异」的恒定位点(见 FAQ 7)。

### 3. FastTree 和 IQ-TREE 结果差在哪？

FastTree 快但**没有 ASC 校正**，SNP 数据下分支长度会低估；IQ-TREE 支持 ASC 且能自动选模型，精度更高。需要严谨的分支长度和支持值，用 IQ-TREE(默认)。

### 4. 输入 VCF 需要先做 LD 剪枝吗？

本工具**不做 LD 剪枝和位点 QC**，假设输入 VCF 已质控。紧密连锁的 SNP 彼此不独立，会影响「假设位点独立」的建树结果——建议先做 LD 剪枝再传入。

### 5. 为什么我的位点被过滤了？

VCF 转 FASTA 阶段只保留二等位单碱基 SNP，多等位、InDel、样本数不足(`--min-samples-locus`)的位点会被跳过，日志里会打印各类过滤统计。

### 6. IQ-TREE 3 报 "Unknown sequence type" 是怎么回事？

这是 IQ-TREE 3.x 的行为变化：IQ-TREE 2 能自动把含大量简并碱基(R/Y/S/W 等)和缺失(N)的比对识别为 DNA；IQ-TREE 3.1.x 的自动识别严格得多——**全比对里非 A/C/G/T 字符占比约 10% 以上时，它放弃识别**并报 `Unknown sequence type`。VCF 转来的比对天然富含简并码和缺失，很容易触发(尤其含许多公共数据样本时)。

本工具自 v1.0.1 起在 IQ-TREE 后端**显式传 `-st DNA`**，正常使用不会再遇到此问题。如果你手动跑 iqtree 撞上这个报错，命令里加 `-st DNA` 即可。

### 7. 日志里「移除仅杂合变异位点」是什么意思？

有些位点在 VCF 里标记为变异，但次要等位基因**只以杂合形式出现**(基因型 0/1)。IQ-TREE 把简并碱基当缺失看待，这类位点在它眼里全部样本都是同一个碱基——属于「恒定位点」：对建树**零信息贡献**，还会让 IQ-TREE 3 的 +ASC 校正直接报错(`Invalid use of +ASC because of ... invariant sites`)。

因此本工具在写 FASTA 前自动移除这类位点(日志打印移除数量，如 417 样本的实例数据移除了 26,639/147,693 个)。**移除无科学损失**——被移掉的位点本来就提供不了任何建树信息，移除后对齐长度变短、建树更快。
