# 系统发育树构建 | Phylogenetic Tree Builder (MAFFT + FastTree)

**使用 MAFFT 进行多序列比对，FastTree 构建最大似然系统发育树 | Perform multiple sequence alignment with MAFFT and build a maximum-likelihood phylogenetic tree with FastTree.**

## 功能概述 | Overview

- 一键完成"清理序列 → MAFFT 比对 → FastTree 建树"完整流程
- 自动检测蛋白/核酸序列类型，也可手动指定
- 自动清洗序列 ID 中的特殊字符，并生成 ID 映射文件
- 支持自定义 MAFFT 与 FastTree 额外参数及可执行文件路径
- 输出 Newick 格式系统发育树

## 快速开始 | Quick Start

```bash
# 自动检测序列类型并建树
biopytools mafft-fasttree -i sequences.fa -o phylo_results

# 指定为蛋白序列并自定义MAFFT参数
biopytools mafft-fasttree -i proteins.fa -o results --seq-type protein \
    --mafft-params "--auto --maxiterate 1000"
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入序列文件 (FASTA 格式) |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `./phylo_output` | 输出目录 |
| `-t, --threads` | `12` | 线程数 |
| `--seq-type` | 自动检测 | 序列类型：`protein` 或 `nucleotide` |
| `--mafft-params` | `--auto` | MAFFT 额外参数 |
| `--fasttree-params` | 空 | FastTree 额外参数（蛋白建树可加 `-lg`） |
| `--mafft-path` | `mafft` | MAFFT 可执行文件路径 |
| `--fasttree-path` | `fasttree` | FastTree 可执行文件路径（部分系统为 `FastTree` 或 `FastTreeMP`） |

（运行 `biopytools mafft-fasttree -h` 查看完整参数列表）

## 输出 | Output

输出目录下生成（`{base}` 为输入文件主名）：

- `{base}.cleaned.fa`：清理 ID 后的序列
- `{base}.mafft.fa`：MAFFT 比对结果
- `{base}.fasttree.nwk`：FastTree 构建的 Newick 树文件
- `{base}.id_mapping.txt`：原始 ID 与清理后 ID 的映射表（如有重命名）
- 运行日志文件

## 依赖 | Dependencies

- MAFFT（>= 7）
- FastTree（建议 FastTreeMP 以支持多线程）
- Python 标准库

## 引用 | Citation

- Katoh K. & Standley D.M. MAFFT multiple sequence alignment software version 7. *Molecular Biology and Evolution*, 2013.
- Price M.N. et al. FastTree 2 — Approximately Maximum-Likelihood Trees for Large Alignments. *PLoS ONE*, 2010.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
