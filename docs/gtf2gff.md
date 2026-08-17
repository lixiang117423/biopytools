# gtf2gff - GTF 转 GFF | GTF to GFF Converter

一句话理解：**把 GTF 格式的注释文件转换成 GFF 格式，并按需清理属性、移除内含子、给特征加统一 ID 前缀**，方便用只认 GFF 的工具继续处理。

## 功能概述 | Overview

- 将 GTF 注释转换为 GFF 文件（-i 输入 GTF、-o 输出 GFF）
- 可选移除 intron 特征（--remove-introns）
- 可选保留全部原始属性（--keep-all-attributes）、跳过属性清理（--no-clean）
- 可选指定 ID 前缀（-p）与物种缩写（-s）
- 注意：当前代码库中该命令的 CLI 包装器已注册，但其底层实现模块（biopytools.gtf_to_gff）尚未存在，命令实际运行会报导入错误，见 FAQ

## 快速开始 | Quick Start

```bash
biopytools gtf2gff -i input.gtf -o output.gff
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| GTF | 一种基因组注释格式，与 GFF 类似但属性用键值对加引号表示，常用于转录组注释 |
| GFF | 另一种更通用的注释格式，部分工具只认 GFF，所以需要转换 |
| intron | 内含子，转录本里被剪掉、不出现在成熟 mRNA 里的部分 |
| 属性 attribute | 每个特征第 9 列携带的信息（ID、Parent、gene_id 等） |
| ID 前缀 prefix | 给特征 ID 统一加的前缀，如 CDRT、AGIS |
| 物种缩写 species | 物种短名，如 Ov、Os，拼进 ID 里 |

## 输入 | Input

### GTF 文件

标准 GTF 格式，由 -i 指定；-o 指定输出的 GFF 文件路径。

```text
chr1    source  exon    1000    1200    .    +    .    gene_id "gene1"; transcript_id "gene1.1";
chr1    source  exon    1300    1500    .    +    .    gene_id "gene1"; transcript_id "gene1.1";
```

## 参数说明 | Parameters

### 输入输出 | Input & output

**通俗理解|In plain words:** -i 输入 GTF、-o 输出 GFF，两个必填。这是转换的「从哪来、到哪去」。

### 转换行为 | Conversion behavior

**通俗理解|In plain words:** --remove-introns 把 intron 特征从输出里去掉；--keep-all-attributes 保留输入里的全部原始属性（否则默认会做清理）；--no-clean 跳过属性清理。一般默认清理即可，除非你有特殊属性要保留。

### ID 前缀与物种 | ID prefix & species

**通俗理解|In plain words:** -p/--prefix 给特征 ID 加统一前缀，-s/--species 指定物种缩写。二者可选，仅在需要给转换后的 ID 加前缀时才用。

### 线程 | Threads

**通俗理解|In plain words:** -t/--threads 线程数，默认 12。一般不用动。

## 输出 | Output

```text
<输出路径>/
└── output.gff                    # 转换后的 GFF 文件 (文件名由 -o 指定)
```

输出为单个 GFF 文件，内容为转换后的注释特征。

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 转换后的 GFF 应包含与输入 GTF 对应的特征，结构一致、格式变为 GFF。

- 检查特征类型与坐标是否与输入 GTF 一致（除非用了 --remove-introns）
- 若指定了 -p/-s，确认输出 ID 带上了相应前缀

## 参数选择建议 | Parameter Guidance

- 只想简单 GTF 转 GFF：biopytools gtf2gff -i in.gtf -o out.gff 即可
- 不需要内含子：加 --remove-introns
- 需要保留全部原始属性：加 --keep-all-attributes
- 需要给 ID 加前缀：加 -p <前缀> -s <物种缩写>

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入GTF文件｜Input GTF file |
| `-o, --output` | 必填 |  | 输出GFF文件｜Output GFF file |
| `--remove-introns` | — |  | 移除intron特征｜Remove intron features |
| `--keep-all-attributes` | — |  | 保留所有原始属性｜Keep all original attributes |
| `--no-clean` | — |  | 不清理属性字段｜Do not clean attributes |
| `-p, --prefix` | — |  | ID前缀｜ID prefix (e.g., CDRT, AGIS) |
| `-s, --species` | — |  | 物种缩写｜Species abbreviation (e.g., Ov, Os) |
| `-t, --threads` | `12` |  | 线程数｜Number of threads [default: 12] |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3
- 预期底层实现依赖 biopytools.gtf_to_gff 模块（当前缺失，见 FAQ）

## 常见问题 | FAQ

**Q1：运行 biopytools gtf2gff 报 ImportError（No module named biopytools.gtf_to_gff）怎么办？**
当前代码库中，gtf2gff 命令的 CLI 包装器已注册，但它延迟导入的底层实现模块 biopytools.gtf_to_gff 尚未实现，因此命令实际运行会报导入错误。请等待该模块补齐后再使用，或改用其他 GFF/GTF 处理工具。

**Q2：这个命令有哪些参数？**
-i（输入 GTF，必填）、-o（输出 GFF，必填）、--remove-introns、--keep-all-attributes、--no-clean、-p/--prefix、-s/--species、-t/--threads（默认 12）。

**Q3：能断点续传吗？**
该命令当前未实现底层逻辑，尚无断点续传能力。