# 基因 DNA 序列提取 | Gene DNA Sequence Extraction

**根据 GFF 注释从基因组 FASTA 中提取指定类型的特征序列（基因/CDS/mRNA 等） | Extract feature sequences (gene/CDS/mRNA, etc.) from a genome FASTA based on GFF annotation.**

## 功能概述 | Overview

- 加载基因组 FASTA 与 GFF 注释文件
- 按特征类型（`gene`、`mRNA`、`CDS`、`exon` 等）过滤 GFF
- 根据正负链自动反向互补提取序列
- 支持最小长度过滤，跳过坐标越界或染色体不存在的记录
- 可控制 FASTA 输出行宽，并在 verbose 模式输出处理进度

## 快速开始 | Quick Start

```bash
# 默认提取gene特征
biopytools parse-gene-dna -g genome.fasta -f annotation.gff -o genes.fasta

# 提取CDS并过滤短序列
biopytools parse-gene-dna -g genome.fasta -f annotation.gff -o cds.fasta \
    --feature-type CDS --min-length 100 -v
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-g, --genome` | 基因组 FASTA 文件路径 |
| `-f, --gff` | GFF 注释文件路径 |
| `-o, --output` | 输出 FASTA 文件路径 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--feature-type` | `gene` | 要提取的 GFF 第三列特征类型 |
| `--min-length` | `0` | 最小基因长度过滤（bp） |
| `-t, --threads` | `12` | 线程数 |
| `--line-width` | `60` | 输出 FASTA 序列行宽度（0 表示单行） |
| `-v, --verbose` | off | 显示详细信息，包括跳过的记录 |

（运行 `biopytools parse-gene-dna -h` 查看完整参数列表）

## 输出 | Output

- 单个 FASTA 文件，每条序列的 header 形如：
  `>{ID} {Name} {seqid}:{start}-{end}({strand}) length={len}`
- ID 取自 GFF 属性 `ID`，Name 取自 `Name`（若存在且与 ID 不同）
- 若 GFF 中不存在指定特征类型，程序会报错退出
- 日志中输出成功提取数量与跳过数量

## 依赖 | Dependencies

- Python 标准库，无第三方依赖

## 引用 | Citation

- GFF3 规范：Sequence Ontology Project GFF3 specification.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
