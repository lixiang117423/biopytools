# IQ-TREE系统发育树 | IQ-TREE Phylogenetic Tree

**基于IQ-TREE构建最大似然系统发育树, 支持模型选择、分区分析、UFBoot和一致性因子 | Maximum likelihood phylogenetic tree with model selection and ultrafast bootstrap**

## 功能概述 | Overview

iqtree 模块封装了 IQ-TREE, 用于基于多重序列比对构建最大似然(ML)系统发育树。IQ-TREE支持DNA、蛋白质、密码子模型, 内置ModelFinder进行自动模型选择, 超快速bootstrap(UFBoot2)和SH-aLRT分支支持度评估, 还支持分区分析(多基因联合建树)、一致性因子计算和祖先状态重建。

## 快速开始 | Quick Start

```bash
# 基础用法(自动模型选择)
biopytools iqtree -i alignment.fasta -o tree_results -p my_tree

# 指定模型和bootstrap
biopytools iqtree -i alignment.fasta -o tree_results -p my_tree -m GTR+I+G -b 1000 --outgroup outgroup1,outgroup2

# 分区分析(多基因)
biopytools iqtree -i concat.fasta -o tree_results -p my_tree --partition partitions.nex --partition-mode edge-linked
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入比对文件 |
| `-o, --output` | 输出目录 |
| `-p, --prefix` | 输出文件前缀 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --model` | `None` | 进化模型(不指定则用ModelFinder自动选择) |
| `-t, --threads` | `12` | 线程数 |
| `-b, --bootstrap` | `1000` | Bootstrap重复次数 |
| `--boot-type` | `ufboot` | Bootstrap类型(ufboot/standard) |
| `--save-boot-trees` | `False` | 保存所有bootstrap树 |
| `--outgroup` | `None` | 外群名称(逗号分隔) |
| `--constraint` | `None` | 约束树文件 |
| `--partition` | `None` | 分区文件 |
| `--partition-mode` | `edge-linked` | 分区模式(edge-linked/edge-equal/edge-unlinked) |
| `--concordance` | `None` | 一致性因子基因树文件 |
| `--ancestral` | `False` | 启用祖先状态重建 |
| `--seed` | `None` | 随机种子 |
| `--runs` | `1` | 独立运行次数 |
| `--redo` | `False` | 重新运行分析 |
| `--iqtree-path` | `iqtree` | IQ-TREE软件路径 |

(运行 `biopytools iqtree -h` 查看完整参数列表)

## 输出 | Output

- `{prefix}.treefile`: 最佳树(Newick格式)
- `{prefix}.best_scheme.nex`: 最优模型和分区方案
- `{prefix}.iqtree`: 完整文本报告
- `{prefix}.contree`: 一致性树(若使用bootstrap)
- `{prefix}.cf.tree`: 一致性因子树(若启用--concordance)

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | 输入比对文件｜Input alignment file |
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--prefix, -p` | 必填 |  | 输出文件前缀｜Output file prefix |
| `--model, -m` | — |  | 进化模型｜Evolutionary model |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--bootstrap, -b` | `1000` | int | Bootstrap重复次数｜Bootstrap replicates |
| `--boot-type` | `ufboot` | ufboot/standard | Bootstrap类型｜Bootstrap type |
| `--save-boot-trees` | — |  | 保存所有bootstrap树｜Save all bootstrap trees |
| `--outgroup` | — |  | 外群名称(逗号分隔)｜Outgroup taxon names (comma-separated) |
| `--constraint` | — | Path | 约束树文件｜Constraint tree file |
| `--partition` | — | Path | 分区文件｜Partition file |
| `--partition-mode` | `edge-linked` | edge-linked/edge-equal/edge-unlinked | 分区模式｜Partition mode |
| `--concordance` | — | Path | 一致性因子基因树文件｜Concordance factor gene tree file |
| `--ancestral` | — |  | 启用祖先状态重建｜Enable ancestral state reconstruction |
| `--seed` | — | int | 随机种子｜Random seed |
| `--runs` | `1` | int | 独立运行次数｜Number of independent runs |
| `--redo` | — |  | 重新运行分析｜Redo analysis |
| `--iqtree-path` | `iqtree` |  | IQ-TREE软件路径｜IQ-TREE program path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入比对文件 (FASTA/PHYLIP/NEXUS)｜Input alignment file (FASTA/PHYLIP/NEXUS) |
| `-o, --output-dir` | 必填 |  | 输出目录｜Output directory |
| `-p, --prefix` | 必填 |  | 输出文件前缀｜Output file prefix |
| `-m, --model` | — |  | 进化模型 (不指定则自动选择)｜Evolutionary model (auto-select if not specified) |
| `-t, --threads` | `12` | int | 线程数｜Number of threads (default: 12) |
| `-b, --bootstrap` | `1000` | int | Bootstrap重复次数｜Bootstrap replicates (default: 1000) |
| `--boot-type` | `ufboot` | ufboot/standard | Bootstrap类型｜Bootstrap type (default: ufboot) |
| `--save-boot-trees` | — | store_true | 保存所有bootstrap树到文件｜Save all bootstrap trees to file |
| `--outgroup` | — |  | 外群名称 (多个用逗号分隔)｜Outgroup taxon names (comma-separated) |
| `--constraint` | — |  | 约束树文件｜Constraint tree file |
| `--partition` | — |  | 分区文件｜Partition file |
| `--partition-mode` | `edge-linked` | edge-linked/edge-equal/edge-unlinked | 分区模式｜Partition mode (default: edge-linked) |
| `--concordance` | — |  | 一致性因子分析：基因树文件｜Concordance factor: gene tree file |
| `--ancestral` | — | store_true | 启用祖先状态重建｜Enable ancestral state reconstruction |
| `--seed` | — | int | 随机种子｜Random seed |
| `--runs` | `1` | int | 独立运行次数｜Number of independent runs (default: 1) |
| `--redo` | — | store_true | 重新运行分析｜Redo analysis |
| `--iqtree-path` | — |  | IQ-TREE程序路径｜IQ-TREE program path |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **IQ-TREE 2**: 系统发育推断 (http://www.iqtree.org/)

## 引用 | Citation

- Minh B.Q., Schmidt H.A., Chernomor O., et al. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution. 37(5):1530-1534.
- Hoang D.T., Chernomor O., von Haeseler A., et al. (2018) UFBoot2: improving the ultrafast bootstrap approximation. Molecular Biology and Evolution. 35(2):518-522.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
