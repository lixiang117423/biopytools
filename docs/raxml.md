# RAxML系统发育树 | RAxML Phylogenetic Analysis

**基于RAxML构建最大似然系统发育树, 支持快速bootstrap和收敛标准 | Maximum likelihood phylogeny with rapid bootstrap and convergence criteria**

## 功能概述 | Overview

raxml 模块封装了 RAxML(随机加速最大似然法), 用于构建大尺度数据集的系统发育树。RAxML在超大规模数据(数百万位点)上有出色的计算效率, 支持DNA、蛋白质、二进制模型, 提供快速bootstrap算法、bootstrap收敛自动停止、ML搜索收敛标准等高级功能。输入需为PHYLIP格式比对文件。

## 快速开始 | Quick Start

```bash
# 基础ML树
biopytools raxml -s alignment.phy -n my_tree

# 快速bootstrap + ML搜索
biopytools raxml -s alignment.phy -n my_tree -f a -x 12345 -p 12345 -N 1000 -m GTRGAMMA -o outgroup

# 自动bootstrap收敛
biopytools raxml -s alignment.phy -n my_tree -f a -x 12345 -p 12345 -I autoMRE -m GTRGAMMA
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-s, --sequence-file` | 输入序列文件(PHYLIP格式) |
| `-n, --output-name` | 输出文件名称 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --model` | `GTRGAMMA` | 替换模型(GTRGAMMA/PROTGAMMAWAG等) |
| `-c, --categories` | `25` | 速率异质性类别数 |
| `-e, --likelihood-epsilon` | `0.1` | 似然优化精度 |
| `-f, --algorithm` | `d` | 算法类型(d=hill-climbing,a=rapid bootstrap+ML) |
| `-p, --parsimony-seed` | `None` | 简约法随机种子 |
| `-N, --runs` | `1` | 运行次数或bootstrap次数 |
| `-b, --bootstrap-seed` | `None` | Bootstrap随机种子 |
| `-x, --rapid-bootstrap-seed` | `None` | 快速bootstrap随机种子 |
| `-I, --bootstrap-convergence` | `None` | Bootstrap收敛标准(autoFC/autoMR/autoMRE/autoMRE_IGN) |
| `-B, --bootstop-threshold` | `0.03` | Bootstrap停止阈值 |
| `-k, --print-bootstrap-trees` | `False` | 输出带分支长度的bootstrap树 |
| `-t, --starting-tree` | `None` | 起始树文件 |
| `-g, --constraint-tree` | `None` | 约束树文件 |
| `-o, --outgroup` | `None` | 外群名称(逗号分隔) |
| `-T, --threads` | `12` | 线程数 |
| `-w, --output-dir` | `./raxml_output` | 输出目录 |
| `--raxml-path` | `raxmlHPC-PTHREADS` | RAxML程序路径 |

(运行 `biopytools raxml -h` 查看完整参数列表)

## 输出 | Output

输出到`{output_dir}/`, 文件名以`RAxML_{type}.{output_name}`格式命名:

- `RAxML_bestTree.{name}`: 最佳ML树
- `RAxML_bipartitions.{name}`: 带分支支持度的树
- `RAxML_bootstrap.{name}`: 所有bootstrap树
- `RAxML_info.{name}`: 运行信息和参数

## 依赖 | Dependencies

- **RAxML**: 系统发育推断 (https://github.com/stamatak/standard-RAxML)
- 推荐 `raxmlHPC-PTHREADS` 多线程版本

## 引用 | Citation

- Stamatakis A. (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics. 30(9):1312-1313.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
