# Phobius跨膜拓扑+信号肽预测 | Phobius TM Topology & Signal Peptide

**用Phobius预测蛋白质跨膜螺旋拓扑与信号肽, 输出清洗后的TSV | Predict transmembrane helix topology and signal peptides with Phobius, output a clean TSV**

## 功能概述 | Overview

phobius 模块封装了 [Phobius](https://phobius.sbc.su.se/), 对蛋白质集预测跨膜(TM)拓扑与信号肽, 与 DeepTMHMM 互补的经典算法。模块分别跑 `-short` 与 `-long` 两种输出, 再解析合并为一份清洗 TSV(每蛋白: 是否信号肽、TM 区数等), 并统计信号肽/TM 蛋白数量。

## 快速开始 | Quick Start

```bash
# 基础用法
biopytools phobius -i proteins.fa -o output_dir/

# 指定输出前缀
biopytools phobius -i proteins.fa -o output_dir/ --prefix SAMPLE1
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入蛋白质 FASTA |
| `-o, --output-dir` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--prefix` | 输入文件名 | 输出前缀 |
| `--phobius-path` | `~/miniforge3/envs/phobius_v.1.0.1/bin/phobius.pl` | phobius.pl 路径 |

(运行 `biopytools phobius -h` 查看完整参数列表)

## 输出 | Output

- `{prefix}.phobius.short.txt`: Phobius short 原始输出
- `{prefix}.phobius.long.txt`: Phobius long 原始输出
- `{prefix}.phobius.tsv`: 清洗后的合并 TSV(每蛋白一行)

## 断点续传 | Checkpoint Resume

short / long 输出文件已存在时分别跳过对应步骤。

## 依赖 | Dependencies

- **Phobius**: 跨膜拓扑 + 信号肽预测 (HMM)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
