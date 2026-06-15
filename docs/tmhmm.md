# TMHMM 跨膜螺旋预测模块

**预测蛋白质跨膜螺旋（TMH）结构域 | Predict transmembrane helices in proteins**

## 功能概述 | Overview

TMHMM模块封装了 TMHMM 2.0c 工具，基于隐马尔可夫模型（HMM）预测蛋白质序列中的跨膜螺旋（Transmembrane Helix, TMH）结构域。模块自动运行tmhmm命令，解析带 `#` 注释的原始输出，生成结构化的TSV表格（带英文列名，便于R/Excel处理），并统计跨膜蛋白数量分布。支持断点续传：若输出已存在则自动跳过。

## 快速开始 | Quick Start

```bash
# 基本预测（不生成图形）
biopytools tmhmm -i proteins.fa -o output_dir/

# 生成图形
biopytools tmhmm -i proteins.fa -o output_dir/ --plot

# 指定输出前缀
biopytools tmhmm -i proteins.fa -o output_dir/ --prefix sample1
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入蛋白质FASTA文件 |
| `-o, --output-dir` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--plot` | `False` | 生成图形（默认 `-noplot`，不生成） |
| `--prefix` | 输入文件名（去扩展名） | 输出文件前缀 |

（运行 `biopytools tmhmm -h` 查看完整参数列表）

### 工具路径配置

tmhmm可执行文件路径按以下优先级查找：
1. `TmhmmConfig` 构造时传入的 `tmhmm_path`
2. 环境变量 `TMHMM_PATH`
3. 默认路径：`~/miniforge3/envs/tmmhmm_v.2.0c/bin/tmhmm`

## 输出 | Output

```
output_dir/
├── {prefix}_tmhmm_raw.txt   # tmhmm原始输出（带#注释行）
├── {prefix}_tmhmm.tsv       # 结构化TSV表格 ⭐
└── tmhmm.log                # 运行日志
```

`{prefix}_tmhmm.tsv` 列说明：

| 列名 | 描述 |
|------|------|
| `ID` | 蛋白质ID |
| `Length` | 蛋白长度（氨基酸数） |
| `Pred_TMHs` | 预测的跨膜螺旋数量 |
| `Exp_AAs_in_TMHs` | 跨膜区期望氨基酸数 |
| `Exp_AAs_first60` | 前60个氨基酸中跨膜区期望氨基酸数 |
| `Prob_N_in` | N端在膜内的总概率 |
| `Signal_Seq` | 是否可能存在N端信号序列（`yes` 或空） |

运行结束后日志会输出统计：`0个TMH`、`1个TMH`、`>=2个TMH` 的蛋白分布。

## 依赖 | Dependencies

- TMHMM 2.0c（需单独安装，学术用户可从DTU申请）
- 默认程序路径：`~/miniforge3/envs/tmmhmm_v.2.0c/bin/tmhmm`
- 推荐通过Conda环境隔离安装

## 引用 | Citation

- Krogh, A., Larsson, B., von Heijne, G., & Sonnhammer, E. L. L. Predicting transmembrane protein topology with a hidden Markov model: application to complete genomes. *Journal of Molecular Biology*, 2001, 305(3):567-580.

## 相关链接 | References

- [TMHMM (DTU)](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0)
- [项目主页](https://github.com/lixiang117423/biopytools)
