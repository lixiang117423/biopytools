# TreeMix 群体历史与基因流 | TreeMix Population History & Gene Flow

一句话理解：**用群体间的等位基因频率推断「谁先谁后分家」的群体分化树，并检测群体之间发生过多少次基因流(迁移)事件**。

## 功能概述 | Overview

- 封装 TreeMix，从等位基因频率协方差构建群体迁移图(migration graph)
- 三个子命令：`prepare`(VCF 转 TreeMix 格式)、`run`(m 值扫描 + OptM + 绘图)、`all`(一键完整流程)
- 自动完成 LD 过滤、频率计算、TreeMix 格式转换
- 扫描 m=0..m_max，用 R 包 OptM 确定最优迁移边数
- 绘制群体树图与残差图

## 快速开始 | Quick Start

```bash
biopytools treemix all -i input.vcf.gz -o treemix_output --cluster cluster.txt
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| 群体分化树 | 群体像「家族树」一样逐级分家的先后关系 |
| 基因流/迁移边(m) | 分家后仍发生的「基因互通」事件，TreeMix 用一条边表示一次迁移 |
| m 值 | 假设发生过多少次迁移；m=0 是纯树，m 越大允许的迁移越多 |
| 外群(root) | 用来「定根」的群体，确定整棵树的方向 |
| OptM | 用统计方法挑出「迁移边加几条最划算」的 R 工具包 |
| bootstrap | 反复重抽样看结果稳不稳，衡量支持度 |
| 残差 | 模型「预测值」和「实际值」的差距，越小拟合越好 |

## 输入 | Input

### VCF 文件(prepare/all 用)

标准 VCF(支持 .vcf.gz)，含基因型。

### 分组文件(可选，`--cluster`)

两列：`样本ID` `群体名`；不提供时按样本名分隔符(默认 `_`)自动推断群体：

```text
sample1	popA
sample2	popA
sample3	popB
```

### TreeMix 输入(run 用)

`prepare` 产出的 `.frq.gz` 格式文件(也可用其他工具预先生成)。

## 子命令 | Subcommands

- `prepare`：VCF → 筛选样品 → LD 过滤 → plink 频率计算 → TreeMix 格式(`input.treemix.frq.gz`)
- `run`：对已有 `.frq.gz` 输入做 m 值扫描 + OptM + 绘图
- `all`：prepare + run 一步完成

## 分析流程 | Pipeline

```text
prepare(输入准备):
  VCF → bcftools 筛选样品 → LD 过滤 → plink 频率(--freq) → plink2treemix → input.treemix.frq.gz

run(核心分析):
  阶段1: 扫描 m=0..m_max，每个 m 跑 N 次重复(含 bootstrap)
  阶段2: OptM 确定最优 m 值
  阶段3: 绘制 m=0 与最优 m 的树图和残差图
```

## 输出 | Output

```text
treemix_output/
├── 01_prepare/
│   ├── selected_samples.vcf.gz        # 按分组筛选后的 VCF
│   ├── pop.cov                        # 群体分组文件
│   ├── input.frq.strat.gz             # plink 频率输出
│   ├── input.treemix.frq.gz           # TreeMix 输入格式(核心中间产物)
│   └── pop.order.txt                  # 群体排序
├── 02_treemix/m{m}/                   # 每个 m 值一个目录
│   └── rep_XXXX.{treeout,vertices,edges,cov,llik}.gz  # 各重复的树/边/协方差/似然
├── 03_optm/
│   ├── run_optm.R                     # OptM R 脚本
│   ├── optimal_m.txt                  # 最优 m 值
│   └── optm_selection.pdf             # m 选择图
├── 04_plot/                           # 最终图
│   ├── m0_tree.pdf                    # m=0 纯树
│   ├── m{m}_tree.pdf                  # 最优 m 的树(含迁移边)
│   └── m{m}_residual.pdf              # 残差图
├── treemix_llik_summary.txt           # 各 m 值的似然汇总
├── optimal_m.txt                      # 最优 m 值(根目录副本)
└── 99_logs/treemix.log                # 运行日志
```

## 结果解读 | Interpreting Results

**通俗理解|In plain words:** 看最优 m 的树图：主干是「分家顺序」，树枝上画的带颜色箭头是「迁移事件」，箭头越粗/颜色越深说明迁移权重越大；再看残差图，残差越小说明模型越能把数据解释干净。

- **树图(`m{m}_tree.pdf`)**：节点为群体，长度表示分化程度；迁移边显示基因流方向与权重
- **残差图(`m{m}_residual.pdf`)**：正值/负值残差表示模型在该群体对上的「偏差」；加迁移边后残差普遍变小，说明该边解释了大量信号
- **`optimal_m.txt` / `03_optm/optimal_m.txt`**：OptM 选出的最优迁移边数
- **`treemix_llik_summary.txt`**：各 m 值的对数似然，用于判断「加边收益是否递减」

## 参数选择建议 | Parameter Guidance

- **`--root`**：强烈建议指定外群群体名，否则树没有方向、迁移方向难解读
- **`--m-max`**：默认 10，测试 m=0..10；群体数多或怀疑多次迁移可加大
- **`-k`**：SNP block 大小(默认 500)，控制相邻 SNP 的相关性，**一般不用动**
- **`--bootstrap 1000`/`--replicates 10`**：默认即可；样本少可提高 replicate 数增加稳定性
- **`--noss`**：单样本群体时关闭样本量校正用
- **`--cluster`**：有明确分组就用它；没有则靠样本名分隔符(默认 `_`)自动推断


<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

_未找到 CLI 参数定义|No CLI definitions found_

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- TreeMix(默认 `~/miniforge3/envs/pop/bin/treemix`)
- PLINK(默认 `~/miniforge3/envs/pop/bin/plink`，用于频率计算与 LD 过滤)
- bcftools(默认 `~/miniforge3/envs/align/bin/bcftools`)
- R(默认 `~/miniforge3/envs/treemix_v.1.13/bin/R`，含 OptM 包)
- plotting_funcs.R(TreeMix 自带绘图脚本，可 `--plotting-funcs` 指定)

## 常见问题 | FAQ

**Q1：换参数重跑，为什么部分步骤被跳过？**
run 阶段有断点续传：某个 m 值下已完成的重复(`rep_XXXX.treeout.gz`)会自动跳过。换 `--bootstrap`、`-k` 等参数重跑旧目录前，先删除 `02_treemix/` 下对应产物，否则复用旧结果。

**Q2：必须指定外群吗？**
不是硬性必须，但强烈建议。不指定 `--root`，TreeMix 也能跑，但树的方向和迁移方向难以确定，解读时容易出错。

**Q3：OptM 失败会怎样？**
程序会跳过绘图但保留扫描结果(`02_treemix/` 和 `treemix_llik_summary.txt`)，并在日志里提示。常见原因是 R 环境缺 OptM 包，检查 `--r-path` 指向的环境。

**Q4：prepare 和 run 能分开跑吗？**
能。先 `prepare` 生成 `input.treemix.frq.gz`，之后随时 `run -i input.treemix.frq.gz` 复用它，不必每次从 VCF 重做。