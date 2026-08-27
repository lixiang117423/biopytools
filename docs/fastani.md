# fastANI 全基因组 ANI 计算

给一批基因组两两算"拼写相似度平均分"(ANI):分数 95% 以上通常就是同一个种,80% 以下工具不报告。解决"我手里这些分离物到底是不是同一个种"的问题。

## 功能概述 | Overview

- 一条命令算出全部基因组两两 ANI,输出对称矩阵 + 每个基因组"最像谁"表
- 支持两种模式:全部互比(all-vs-all)与新基因组 vs 参考集(query-vs-ref)
- **大数据集内存保护**:all-vs-all 基因组数超过阈值(默认 100)时自动切换逐轮 1-vs-all 遍历(每轮 1 个 query vs 全部),内存友好,断点续传按批次粒度
- 速度快(比对-free,MinHash 映射),成百上千基因组可跑

## 快速开始 | Quick Start

```bash
biopytools fastani -i genome_dir/ -o output_dir/
```

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗解释 |
|---|---|
| ANI | 两个基因组"拼写相似度平均分",满分 100。95-96% 以上通常同种 |
| 片段(fragment) | 工具把基因组切成 ~3000bp 小段,逐段找对应关系再平均 |
| aligned_fraction | 有多少比例的片段找到了对应:ANI 高但覆盖低 = 近缘但不完整 |
| NA | 两基因组 ANI 低于 ~80%,工具不报告——"远亲",需换蛋白层面(AAI)分析 |
| all-vs-all | 目录里所有基因组两两互比,出完整对称矩阵 |

## 输入 | Input

基因组 FASTA(`.fa/.fna/.fasta`,支持 `.gz`),三种给法:
1. 目录:工具自动收集所有 fasta 文件
2. 列表文件(`.txt`):每行一个基因组路径,`#` 开头为注释
3. 单个 fasta 文件

格式示例(列表文件):

```text
/path/R0590-6.fna
/path/R0601-2.fa.gz
# 注释行|comment
```

要求:同侧文件主名唯一(`x.fa` 与 `x.fna` 同放会报错);建议 N50 ≥ 10kb。

## 参数说明 | Parameters

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | — | Path | all-vs-all输入(目录/列表/单fasta;与-q/-r互斥)｜all-vs-all input (dir/list/single fasta) |
| `--query, -q` | — | Path | query侧输入(定向模式)｜Query side (directional) |
| `--reference, -r` | — | Path | reference侧输入(定向模式)｜Reference side (directional) |
| `--output-dir, -o` | `./fastani_output` | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Thread count |
| `--kmer, -k` | `16` | int | k-mer大小(<=16)｜K-mer size (<=16) |
| `--frag-len` | `3000` | int | 片段长度｜Fragment length |
| `--min-fraction` | `0.2` | float | 信任ANI的最小共享比例｜Min shared fraction to trust ANI |
| `--iterated/--no-iterated` | `True` |  | 大数据集自动切换逐轮1-vs-all(默认开)｜Auto switch to iterated 1-vs-all for large sets (default on) |
| `--iterated-threshold` | `100` | int | 触发遍历的基因组数阈值(默认100)｜Genome count threshold for iterated mode (default 100) |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  | all-vs-all输入(目录/列表文件/单fasta,与-q/-r互斥)｜all-vs-all input (dir/list/single fasta) |
| `-q, --query` | — |  | query侧输入(定向模式)｜Query side (directional mode) |
| `-r, --reference` | — |  | reference侧输入(定向模式)｜Reference side (directional mode) |
| `-o, --output-dir` | `./fastani_output` |  | 输出目录｜Output directory |
| `-t, --threads` | `12` | int | 线程数｜Thread count |
| `-k, --kmer` | `16` | int | k-mer大小(<=16)｜K-mer size (<=16) |
| `--frag-len` | `3000` | int | 片段长度｜Fragment length |
| `--min-fraction` | `0.2` | float | 信任ANI的最小共享比例｜Min shared fraction to trust ANI |
| `--iterated` | `True` | store_true |  |
| `--no-iterated` | — | store_false | 关闭大数据集自动遍历(强制all-vs-all)｜Disable auto iterated mode (force all-vs-all) |
| `--iterated-threshold` | `100` | int | 触发遍历的基因组数阈值(默认100)｜Genome count threshold for iterated mode (default 100) |
| `--log-level` | `INFO` | DEBUG/INFO/WARNING/ERROR | 日志级别｜Log level |

<!-- END PARAMS:auto -->

**通俗理解|In plain words:** 模式参数管"谁和谁比",数值参数是 fastANI 的内部精度旋钮——**一般都不用动**;`-t` 线程数按节点配额调大即可。

## 分析流程 | Pipeline

```
基因组收集 → 01_fastani(fastANI 运行+原始输出) → 02_results(矩阵+最近邻) → 99_logs
```

## 输出 | Output

```
output_dir/
├── 00_pipeline_info/software_versions.yml   # 软件版本与参数快照
├── 01_fastani/
│   ├── genome_list_ql.txt / genome_list_rl.txt  # 参与比对的基因组清单
│   ├── fastani.out        # 原始 5 列结果(query ref ANI 比上片段 总片段)
│   └── fastani.out.matrix     # phylip 下三角矩阵(官方原始格式)
├── 02_results/
│   ├── ani_matrix.tsv     # 行=genome 列=genome 的对称矩阵(<80% 记 NA)
│   └── nearest_genome.tsv # 每个基因组最像谁 + ANI + 覆盖度
└── 99_logs/               # 全量/输出/错误三份日志
```

## 结果解读 | Interpreting Results

- `ani_matrix.tsv`:行列同序;对角线 100;双向取平均(与官方 --matrix 一致)
  - 好判据:同种分离物之间 99%+;同属不同种 85-95%;NA = 低于 80% 的远缘
- `nearest_genome.tsv`:排序按相似度从高到低;`aligned_fraction` 高说明整体可比,低说明只有局部相似(如 shared plasmid)
- 全 NA 的行会有 WARNING 日志——检查是不是拿错了基因组(过于 divergent)

## 参数选择建议 | Parameter Guidance

- 一批分离物分型:`-i 目录` 默认参数即可
- 新基因组鉴定:`-q 新基因组 -r 参考集目录`
- 覆盖度低想看局部相似:不用调参,看 aligned_fraction 列即可

## 依赖 | Dependencies

- fastANI 1.34(conda env `pop`,`mamba install -n pop fastani`)
- biopytools(本模块)

## 常见问题 | FAQ

- **为什么 (A,B) 和 (B,A) 分数不完全一样?** fastANI 有方向性(query 切片段),矩阵已自动双向平均
- **为什么有 NA?** ANI < ~80% 工具不报告;请改用蛋白层面比较(AAI)或调整研究问题
- **多少基因组能跑?** all-vs-all 是平方级。**超过 100 个基因组(默认阈值)模块自动切换为逐轮 1-vs-all 遍历模式**——每轮只把 1 个基因组当 query 与全部比对,内存峰值从"全部基因组草图驻留"降到"单个基因组草图",884 个 ~1GB 基因组 800G 内存即可跑完;可用 `--iterated-threshold` 调整阈值,`--no-iterated` 强制 all-vs-all
- **重跑会重复计算吗?** 不会,`01_fastani/fastani.out` 已存在时自动跳过,只重算后处理
