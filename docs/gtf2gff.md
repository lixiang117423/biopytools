# GTF 转 GFF | GTF to GFF Converter

**将 GTF 注释文件转换为规范的 GFF3 格式，支持属性清理、intron 移除与 ID 重命名 | Convert GTF annotation to well-formed GFF3 with attribute cleaning, intron removal and ID renaming**

## 功能概述 | Overview

`gtf2gff` 用于将常见基因组浏览器/注释工具产出的 GTF 文件（如 StringTie、Scallop、Cufflinks、MAKER 等）转换为 GFF3 规范格式。GFF3 采用 `key=value;key=value` 的属性语法并要求显式的 `Parent` 层级关系（gene → mRNA → exon/CDS），而 GTF 使用 `key "value";` 语法、靠 `gene_id`/`transcript_id` 关联特征，本工具完成这一结构与语法层面的完整转换。

转换过程中可同时进行属性清理（去除冗余字段、统一命名）、移除 `intron` 特征行（部分注释流程会显式写出 intron 行，GFF3 中通常不需要）、以及为 gene/mRNA 重新生成统一前缀的 ID。`--prefix` 与 `--species` 需同时提供或同时省略，用于按物种缩写批量重命名 ID（例如前缀 `Ov`、物种 `Os`），便于多个注释文件合并前的统一。

## 快速开始 | Quick Start

```bash
# 基本转换
biopytools gtf2gff -i input.gtf -o output.gff

# 带 ID 重命名与多线程
biopytools gtf2gff -i stringtie.gtf -o cleaned.gff3 \
    -p CDRT -s Os -t 16 --remove-introns
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入 GTF 文件 |
| `-o, --output` | 输出 GFF 文件 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-t, --threads` | `12` | 线程数 |
| `--remove-introns` | 关 | 移除 intron 特征行 |
| `--keep-all-attributes` | 关 | 保留所有原始属性（不做精简）|
| `--no-clean` | 关 | 不清理属性字段 |
| `-p, --prefix` | 无 | ID 前缀（如 `CDRT`、`AGIS`），需与 `-s` 同时使用 |
| `-s, --species` | 无 | 物种缩写（如 `Ov`、`Os`），需与 `-p` 同时使用 |

（运行 `biopytools gtf2gff -h` 查看完整参数列表）

约束：`--prefix` 和 `--species` 必须同时提供，或同时不提供。

## 输出 | Output

```
./
├── output.gff3         # 标准 GFF3 文件，含 gene/mRNA/exon/CDS 层级
└── gtf_to_gff.log      # 转换日志（位于输出文件同目录）
```

日志中会输出基因数、转录本数、外显子数、CDS 数等统计信息。

## 依赖 | Dependencies

- Python 3.7+（纯 Python 实现，无第三方依赖）

## 引用 | Citation

- GFF3 规范：Sequence Ontology Project, Generic Feature Format Version 3 (https://github.com/The-Sequence-Ontology/Specifications)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
