# 基因分组核苷酸多样性 | Nucleotide Diversity per Gene Group (pi4gene)

**按分组提取序列、MAFFT比对、计算核苷酸多样性pi | Extract sequences by group, align with MAFFT, calculate nucleotide diversity (pi)**

## 功能概述 | Overview

pi4gene 用于计算基因分组的核苷酸多样性(nucleotide diversity, pi)。用户提供一个序列FASTA文件和一个分组ID文件(两列:分组名 + 序列ID), 工具会按分组提取序列, 使用 MAFFT 进行多重序列比对, 然后基于比对结果计算每个分组的 pi 值。适用于基因家族、同源基因、等位基因多样性等分析场景。

## 快速开始 | Quick Start

```bash
# 标准用法
biopytools pi4gene -i genes.fasta -d groups.txt -o pi4gene_output

# 多线程加速
biopytools pi4gene -i genes.fasta -d groups.txt -o pi4gene_output -t 24
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入序列FASTA文件路径 |
| `-d, --id-file` | 分组ID文件(第一列分组名, 第二列序列ID) |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./pi4gene_output` | 输出目录 |
| `-t, --threads` | `12` | 线程数 |
| `--mafft-path` | `None` | MAFFT可执行文件路径(默认从PATH查找) |

(运行 `biopytools pi4gene -h` 查看完整参数列表)

## 输出 | Output

输出目录包含每个分组的比对结果和pi统计表:

```
pi4gene_output/
├── alignments/             # 每个分组的MAFFT比对结果
│   ├── group1.fasta
│   └── group2.fasta
├── pi_results.csv          # pi统计结果表
└── pi4gene.log             # 运行日志
```

## 依赖 | Dependencies

- **MAFFT**: 多重序列比对 (https://mafft.cbrc.jp/alignment/software/)
- **Python库**: biopython (用于序列处理和pi计算)
- **标准库**: multiprocessing (并行化)

## 引用 | Citation

- Katoh K., Standley D.M. (2013) MAFFT multiple sequence alignment software version 7. Molecular Biology and Evolution. 30(4):772-780.
- Nei M., Li W.H. (1979) Mathematical model for studying genetic variation in terms of restriction endonucleases. PNAS. 76(10):5269-5273.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
