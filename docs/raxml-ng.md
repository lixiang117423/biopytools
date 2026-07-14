# RAxML-NG系统发育树 | RAxML-NG Phylogenetic Tree

**基于RAxML-NG用最大似然法构建系统发育树, 支持自适应模型选择、bootstrap分支支持度与断点续传 | Maximum likelihood phylogenetic tree with adaptive model selection, bootstrap support and checkpointing**

## 功能概述 | Overview

raxml-ng 模块封装了 RAxML-NG, 用于基于多重序列比对用最大似然(ML)法构建系统发育树。RAxML-NG是经典RAxML的继任者, 采用SPR启发式搜索, 内置自适应搜索策略、自动模型选择和快速分支支持度计算(FBP/TBE/SH/EBG等), 支持DNA、蛋白质和二进制数据, 并原生支持断点续传。

支持三种核心工作流:
- **all(默认)**: 一站式 ML搜索 + bootstrap + 支持值映射, 适合发表级分析
- **search**: 仅ML搜索, 快速建树
- **support**: 给已有参考树映射bootstrap支持值

## 快速开始 | Quick Start

```bash
# 基础用法(一站式: ML + 1000次bootstrap + 支持值, 自适应模型)
biopytools raxml-ng -i alignment.fasta -o tree_results -p my_tree

# 指定模型与bootstrap数, 设置外群
biopytools raxml-ng -i alignment.fasta -o tree_results -p my_tree -m GTR+G -b 1000 --outgroup outgroup1,outgroup2

# 仅快速ML搜索(无bootstrap)
biopytools raxml-ng -i alignment.fasta -o tree_results -p my_tree --mode search
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入比对文件(FASTA/PHYLIP) |
| `-o, --output` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-p, --prefix` | 输入文件名 | 输出文件前缀(决定输出名与续传键) |
| `--mode` | `all` | 分析模式(all/search/support) |
| `-m, --model` | `None` | 进化模型(不指定则RAxML-NG自适应选择) |
| `-t, --threads` | `12` | 线程数 |
| `-b, --bs-trees` | `1000` | Bootstrap重复数,可传`autoMRE{1000}`走自动收敛 |
| `--bs-metric` | `fbp` | 支持值类型(fbp/tbe/sh/ebg/rbs) |
| `--tree` | `None` | 起始树(search)/参考树(support) |
| `--outgroup` | `None` | 外群名称(逗号分隔,仅画图) |
| `--seed` | `None` | 随机种子(可复现) |
| `--redo` | `False` | 覆盖已有结果,忽略checkpoint |
| `--raxml-ng-path` | 默认二进制 | RAxML-NG软件路径 |

(运行 `biopytools raxml-ng -h` 查看完整参数列表)

## 模式说明 | Modes

| 模式 | 底层命令 | 用途 |
|------|----------|------|
| `all`(默认) | `raxml-ng --all` | ML搜索 + bootstrap + 支持值映射,一站式中发表级 |
| `search` | `raxml-ng --search` | 仅ML搜索,快速建树 |
| `support` | `raxml-ng --support` | 给参考树映射bootstrap支持值(需`--tree`参考树与`-b`bootstrap树文件) |

## 输出 | Output

RAxML-NG 以 `{prefix}.raxml.*` 形式输出到输出目录(扁平结构,保留原生断点续传):

- `{prefix}.raxml.bestTree`: 最佳ML树(Newick)
- `{prefix}.raxml.supportTree`: 带支持值的树(all模式)
- `{prefix}.raxml.bootstraps`: bootstrap树
- `{prefix}.raxml.log`: RAxML-NG运行日志
- `{prefix}.raxml.ckp`: checkpoint(断点续传用)
- `99_logs/raxml_ng_analysis.log`: biopytools运行日志
- `00_pipeline_info/software_versions.yml`: 软件版本与参数记录

## 断点续传 | Checkpoint Resume

默认不传`--redo`, RAxML-NG自动检测checkpoint并续传(中断后重跑同一命令即可继续)。加`--redo`则覆盖已有结果、忽略checkpoint。

## 路径配置 | Path Configuration

二进制路径按优先级: 环境变量`RAXML_NG_PATH` > `~/.config/biopytools/config.yml`的`tools.raxml_ng` > 代码默认值。默认指向`~/software/RAxML/raxml-ng_v2.0.2_linux_x86_64/raxml-ng`。

## 依赖 | Dependencies

- **RAxML-NG 2.0.2**: 系统发育推断 (https://github.com/amkozlov/raxml-ng)

## 引用 | Citation

- Kozlov A.M., Darriba D., Flouri T., Morel B., Stamatakis A. (2019) RAxML-NG: A fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35(21):4453-4455.
- Togkousidis A. et al. (2023) Adaptive strategies for RAxML-NG. *Mol. Biol. Evol.* 40(10):msad227.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
- [RAxML-NG wiki](https://codeberg.org/amkozlov/raxml-ng/wiki)
