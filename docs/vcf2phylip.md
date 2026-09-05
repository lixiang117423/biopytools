# vcf2phylip - VCF 转建树格式 | vcf2phylip VCF to PHYLIP/FASTA/NEXUS

一句话理解：**把 VCF 里的 SNP 数据转成下游建树软件(PHYLIP / FASTA / NEXUS)能直接读的比对矩阵文件。**

## 功能概述 | Overview

- 输出 4 种格式：PHYLIP(默认)、FASTA、NEXUS、二进制 NEXUS
- 只保留「二等位单碱基 SNP」，多等位和 InDel 自动剔除
- 可按「位点最少样本数」过滤缺失过多的位点
- 支持随机解析杂合基因型、记录通过筛选的位点坐标
- 纯 Python 实现，无任何外部软件依赖

## 快速开始 | Quick Start

```bash
biopytools vcf2phylip -i variants.vcf -o converted_results
```

最小输入：一个 VCF 文件(`-i`)。默认输出 PHYLIP 格式(`.phy`)到 `-o` 目录。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| VCF | 记录「每个样本在每个位点是什么基因型」的标准变异文件 |
| SNP | 单碱基变异，即某个位置 A 变成了 G |
| 二等位位点 | 一个位置只有「参考/替代」两种可能的碱基(最常见) |
| 杂合基因型 | 样本在该位点两条染色体不一样(如 0/1) |
| IUPAC 码 | 用单个字母表示「可能是哪几个碱基」的模糊编码(如 R=A 或 G) |
| PHYLIP / NEXUS / FASTA | 三种建树软件常用的比对存储格式 |

## 输入 | Input

### VCF 文件

支持 `.vcf` 和 `.vcf.gz`(gzip 压缩)。需含 `#CHROM` 样本头行。示例：

```text
##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1  sample2  sample3
Chr1  100  .  A  G  .  PASS  .  GT  0/0  0/1  1/1
```

只有二等位单碱基 SNP 会被保留；多等位位点、InDel(长度不是 1)会被跳过。

## 参数说明 | Parameters

### 必需参数 | Required

**通俗理解|In plain words:** 只有 `-i` 是必填(输入 VCF)。`-o/--output` 默认 `./converted_output`，`--output-prefix` 是输出文件名前缀(默认取输入文件名去掉 `.vcf`)。

### 转换参数 | Conversion

**通俗理解|In plain words:** `-m/--min-samples-locus` 是「一个位点至少要有多少个样本有数据才保留」(默认 4)——缺失太严重的位点信息量不足，直接扔掉；**调大=更严格(留更少位点)，调小=更宽松**。`-g/--outgroup` 指定外群样本名(可选)，指定后外群序列会排在最前。

### 输出格式 | Output format

**通俗理解|In plain words:** 默认只输出 PHYLIP(`.phy`)。要别的格式就加对应开关：`-f`(FASTA)、`-n`(NEXUS)、`-b`(二进制 NEXUS)。`-p/--phylip-disable` 可关掉默认的 PHYLIP 输出——但**至少得保留一种输出格式**，否则会报错。

### 处理选项 | Processing options

**通俗理解|In plain words:** `-r/--resolve-IUPAC` 把杂合位点「随机选一个碱基」变成纯合(注意**结果不可复现**，因为用了随机选择)。`-w/--write-used-sites` 额外输出一张「哪些位点通过了筛选」的坐标表。`-t/--threads` 是线程数——**本工具当前实际未用它做并行**，仅作日志记录。

## 输出 | Output

```text
converted_results/
├── <prefix>_min4.phy          # PHYLIP 比对矩阵(默认)
├── <prefix>_min4.fasta        # FASTA 矩阵(仅 -f)
├── <prefix>_min4.nexus        # NEXUS 矩阵(仅 -n)
├── <prefix>_min4.bin.nexus    # 二进制 NEXUS 矩阵(仅 -b)
├── <prefix>_min4_used_sites.tsv  # 通过筛选的位点坐标(仅 -w)
└── vcf_conversion.log         # 运行日志
```

其中 `<prefix>` 是 `--output-prefix`(默认输入 basename)，`min4` 表示 `--min-samples-locus` 取 4。

## 结果解读 | Interpreting Results

### 1. PHYLIP 矩阵(.phy)

**通俗理解|In plain words:** 首行是 `样本数 位点数`，之后每行一个样本的序列。这是给 RAxML 等经典软件的标准输入，可直接喂给建树工具。

```text
3 1000
sample1  ACGTACGT...
sample2  ACGTACGT...
sample3  ACGTACGT...
```

### 2. FASTA / NEXUS 矩阵

**通俗理解|In plain words:** 同样的数据换种格式存，FASTA 给 IQ-TREE 等、NEXUS 给 MrBayes 等。二进制 NEXUS(`.bin.nexus`)用 0/1/2 编码纯合/杂合，只含二等位 SNP。

### 3. 位点坐标表(_used_sites.tsv)

**通俗理解|In plain words:** 记录「哪些位点通过了筛选」(染色体、位置、样本数)。需要把结果里的位点映射回原始基因组坐标时用它。

## 参数选择建议 | Parameter Guidance

- **只给 RAxML 用**：默认输出 PHYLIP 即可
- **要多种格式**：默认 PHYLIP 之外，按需加 `-f`(FASTA)、`-n`(NEXUS)
- **位点缺失太多想收紧**：调大 `--min-samples-locus`
- **给树定方向**：用 `-g <外群样本名>` 让外群排最前
- **不要 PHYLIP 只留 FASTA**：`-p -f`(关 PHYLIP + 开 FASTA)

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 | Path | 输入VCF文件路径｜Input VCF file path |
| `--output, -o` | `./converted_output` | Path | 输出目录｜Output directory |
| `--output-prefix` | — | str | 输出文件名前缀｜Output filename prefix |
| `--min-samples-locus, -m` | `4` | int | 位点最少样本数｜Minimum samples per locus |
| `--outgroup, -g` | `` | str | 外群样本名称｜Outgroup sample name |
| `--phylip-disable, -p` | — |  | 禁用PHYLIP输出｜Disable PHYLIP output |
| `--fasta, -f` | — |  | 启用FASTA输出｜Enable FASTA output |
| `--nexus, -n` | — |  | 启用NEXUS输出｜Enable NEXUS output |
| `--nexus-binary, -b` | — |  | 启用二进制NEXUS输出｜Enable binary NEXUS output |
| `--resolve-IUPAC, -r` | — |  | 随机解析杂合子基因型｜Resolve heterozygous genotypes |
| `--write-used-sites, -w` | — |  | 保存筛选通过的位点坐标｜Save used sites coordinates |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件路径｜Input VCF file path |
| `-o, --output` | `./converted_output` |  | 输出目录｜Output directory |
| `--output-prefix` | — |  | 输出文件名前缀｜Output filename prefix |
| `-m, --min-samples-locus` | `4` | int | 位点最少样本数｜Minimum samples per locus |
| `-g, --outgroup` | `` |  | 外群样本名称｜Outgroup sample name |
| `-p, --phylip-disable` | — | store_true | 禁用PHYLIP输出｜Disable PHYLIP output |
| `-f, --fasta` | — | store_true | 启用FASTA输出｜Enable FASTA output |
| `-n, --nexus` | — | store_true | 启用NEXUS输出｜Enable NEXUS output |
| `-b, --nexus-binary` | — | store_true | 启用二进制NEXUS输出｜Enable binary NEXUS output |
| `-r, --resolve-IUPAC` | — | store_true | 随机解析杂合子基因型｜Resolve heterozygous genotypes |
| `-w, --write-used-sites` | — | store_true | 保存筛选通过的位点坐标｜Save used sites coordinates |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `-v, --version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

无外部软件依赖(纯 Python 标准库实现：gzip 读压缩、random 解析杂合)。无需 conda 环境。

## 常见问题 | FAQ

### 1. 会断点续传吗？

**不会。** 每次运行都会重新解析 VCF 并重写输出文件。

### 2. 为什么我的位点少了很多？

本工具**只保留二等位单碱基 SNP**：多等位位点、InDel(插入/缺失)、缺失样本数超过 `--min-samples-locus` 的位点都会被跳过。这是正常过滤，日志里会打印各步骤的过滤统计。

### 3. `-r` 随机解析杂合的结果能复现吗？

不能。`--resolve-IUPAC` 用随机选择把杂合变成纯合，每次运行结果可能不同。需要可复现时不要用 `-r`(保留 IUPAC 模糊码即可)。

### 4. `-t` 线程数有用吗？

当前版本 `-t/--threads` 仅记录到日志，**并未真正用于并行处理**。想加速可后续版本关注，或直接不改它。

### 5. 二进制 NEXUS 和普通 NEXUS 有什么区别？

普通 NEXUS 用 IUPAC 碱基(A/C/G/T/R/Y...)表示；二进制 NEXUS(`.bin.nexus`)用 0/1/2 表示纯合参考/杂合/纯合替代，只含二等位 SNP，适合某些专门软件(如 MrBayes 的 SNP 模型)。
