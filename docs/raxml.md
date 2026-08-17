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

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--sequence-file, -s` | 必填 |  | 输入序列文件(PHYLIP格式)｜Input sequence file (PHYLIP format) |
| `--output-name, -n` | 必填 | str | 输出文件名称｜Output file name |
| `--model, -m` | `GTRGAMMA` | str | 替换模型｜Substitution model (GTRGAMMA, PROTGAMMAWAG, etc.) |
| `--categories, -c` | `25` | int | 速率异质性类别数｜Number of rate heterogeneity categories |
| `--likelihood-epsilon, -e` | `0.1` | float | 似然优化精度｜Likelihood optimization precision |
| `--algorithm, -f` | `d` | str | 算法类型｜Algorithm type (d=rapid hill-climbing, a=rapid bootstrap, etc.) |
| `--parsimony-seed, -p` | — | int | 简约法随机种子｜Parsimony random seed |
| `--runs, -N` | `1` | str | 运行次数或bootstrap次数｜Number of runs or bootstrap replicates |
| `--bootstrap-seed, -b` | — | int | Bootstrap随机种子｜Bootstrap random seed |
| `--rapid-bootstrap-seed, -x` | — | int | 快速bootstrap随机种子｜Rapid bootstrap random seed |
| `--bootstrap-convergence, -I` | — | autoFC/autoMR/autoMRE/autoMRE_IGN | Bootstrap收敛标准｜Bootstrap convergence criterion |
| `--bootstop-threshold, -B` | `0.03` | float | Bootstrap停止阈值｜Bootstrap stop threshold |
| `--bootstop-perms` | `100` | int | Bootstrap停止检验次数｜Bootstrap stop test permutations |
| `--print-bootstrap-trees, -k` | — |  | 输出带分支长度的bootstrap树｜Print bootstrap trees with branch lengths |
| `--starting-tree, -t` | — | Path | 起始树文件｜Starting tree file |
| `--constraint-tree, -g` | — | Path | 约束树文件｜Constraint tree file |
| `--outgroup, -o` | — | str | 外群名称(逗号分隔多个)｜Outgroup name(s) (comma-separated) |
| `--threads, -T` | `12` | int | 线程数｜Number of threads |
| `--memory-saving, -U` | — |  | 启用内存节省模式｜Enable memory saving mode |
| `--ml-search-convergence, -D` | — |  | 启用ML搜索收敛标准｜Enable ML search convergence criterion |
| `--random-starting-tree, -d` | — |  | 使用随机起始树｜Use random starting tree |
| `--disable-rate-heterogeneity, -V` | — |  | 禁用速率异质性模型｜Disable rate heterogeneity model |
| `--gamma-median, -u` | — |  | 使用GAMMA模型中位数｜Use median for GAMMA model |
| `--disable-pattern-compression, -H` | — |  | 禁用模式压缩｜Disable pattern compression |
| `--output-dir, -w` | `./raxml_output` | Path | 输出目录｜Output directory |
| `--raxml-path` | `raxmlHPC-PTHREADS` | str | RAxML程序路径｜RAxML program path |
| `--no-seq-check` | — |  | 跳过序列检查｜Skip sequence checking |
| `--silent` | — |  | 静默模式｜Silent mode |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-s, --sequence-file` | 必填 |  | 输入序列文件(PHYLIP格式)｜Input sequence file (PHYLIP format) |
| `-n, --output-name` | 必填 |  | 输出文件名称｜Output file name |
| `-m, --model` | `GTRGAMMA` |  | 替换模型｜Substitution model (GTRGAMMA, PROTGAMMAWAG, etc.) |
| `-c, --categories` | `25` | int | 速率异质性类别数｜Number of rate heterogeneity categories |
| `-e, --likelihood-epsilon` | `0.1` | float | 似然优化精度｜Likelihood optimization precision |
| `-f, --algorithm` | `d` |  | 算法类型｜Algorithm type (d=rapid hill-climbing, a=rapid bootstrap, etc.) |
| `-p, --parsimony-seed` | — | int | 简约法随机种子｜Parsimony random seed |
| `-#, -N, --runs` | `1` |  | 运行次数或bootstrap次数｜Number of runs or bootstrap replicates |
| `-b, --bootstrap-seed` | — | int | Bootstrap随机种子｜Bootstrap random seed |
| `-x, --rapid-bootstrap-seed` | — | int | 快速bootstrap随机种子｜Rapid bootstrap random seed |
| `-I, --bootstrap-convergence` | — | autoFC/autoMR/autoMRE/autoMRE_IGN | Bootstrap收敛标准｜Bootstrap convergence criterion |
| `-B, --bootstop-threshold` | `0.03` | float | Bootstrap停止阈值｜Bootstrap stop threshold |
| `--bootstop-perms` | `100` | int | Bootstrap停止检验次数｜Bootstrap stop test permutations |
| `-k, --print-bootstrap-trees` | — | store_true | 输出带分支长度的bootstrap树｜Print bootstrap trees with branch lengths |
| `-t, --starting-tree` | — |  | 起始树文件｜Starting tree file |
| `-g, --constraint-tree` | — |  | 约束树文件｜Constraint tree file |
| `-o, --outgroup` | — |  | 外群名称 (逗号分隔多个)｜Outgroup name(s) (comma-separated) |
| `-T, --threads` | `88` | int | 线程数｜Number of threads |
| `-U, --memory-saving` | — | store_true | 启用内存节省模式｜Enable memory saving mode |
| `-D, --ml-search-convergence` | — | store_true | 启用ML搜索收敛标准｜Enable ML search convergence criterion |
| `-d, --random-starting-tree` | — | store_true | 使用随机起始树｜Use random starting tree |
| `-V, --disable-rate-heterogeneity` | — | store_true | 禁用速率异质性模型｜Disable rate heterogeneity model |
| `-u, --gamma-median` | — | store_true | 使用GAMMA模型中位数｜Use median for GAMMA model |
| `-H, --disable-pattern-compression` | — | store_true | 禁用模式压缩｜Disable pattern compression |
| `-w, --output-dir` | `./raxml_output` |  | 输出目录｜Output directory |
| `--raxml-path` | `raxmlHPC-PTHREADS` |  | RAxML程序路径｜RAxML program path |
| `--no-seq-check` | — | store_true | 跳过序列检查｜Skip sequence checking |
| `--silent` | — | store_true | 静默模式｜Silent mode |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **RAxML**: 系统发育推断 (https://github.com/stamatak/standard-RAxML)
- 推荐 `raxmlHPC-PTHREADS` 多线程版本

## 引用 | Citation

- Stamatakis A. (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics. 30(9):1312-1313.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
