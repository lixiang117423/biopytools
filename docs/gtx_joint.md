# GTX Joint Calling 命令生成器 | GTX Joint Calling Command Generator

一句话理解：**它自己不跑分析，而是根据你的一堆 GVCF 文件和参考基因组，自动生成一个「联合变异检测」的 shell 脚本**，你拿到脚本再去提交或运行，解决「多样本联合变异检测命令太长、要按染色体/区间拆分并行」的编排问题。

## 功能概述 | Overview { #overview }

- 输入一个 GVCF 目录(匹配 `*.g.vcf.gz`) + 参考基因组 + GTX 可执行文件
- 自动检查/重建参考基因组索引(samtools faidx)和每个 GVCF 的索引(tabix)
- 按样本数量自动选择三种模式生成脚本：单机模式(≤200 样本)、按染色体拆分(>200 样本)、按区间拆分(>200 样本且指定窗口)
- 只生成脚本不执行，输出可执行脚本 `run_gtx_joint.sh` + 运行日志

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools gtx-joint -r genome.fa -i ./gvcf -o ./output
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| GVCF(gVCF) | 含未变异位点的 VCF，每个位点都有记录，适合多样本联合 |
| joint calling | 把多个样本的变异放在一起重新判读，提高低频变异检出率 |
| faketime | 伪装系统时间的小工具，用于绕过 GTX 的 license 时间校验 |
| 染色体拆分 | 把「一次性处理整条染色体」的大命令切成「一条染色体一条命令」的小块，便于并行 |
| 区间拆分 | 在染色体内部再按固定长度切成更小的窗口(如每 10Mb 一块) |
| tabix 索引 | 给压缩的 VCF 建一个「快速定位」的辅助文件(.tbi) |

## 输入 | Input { #input }

### GVCF 目录

`-i` 指向一个目录，程序只认文件名匹配 `*.g.vcf.gz` 的文件。缺 `.tbi` 索引或索引过时会自动重建。

```text
gvcf/
├── sample1.g.vcf.gz
├── sample1.g.vcf.gz.tbi
├── sample2.g.vcf.gz
└── sample2.g.vcf.gz.tbi
```

### 参考基因组

`-r` 指定参考 FASTA；`.fai` 索引不存在时自动用 samtools faidx 生成。

### GTX 可执行文件

`-g` 指定 GTX 可执行文件(默认 `~/software/gtx/bin/gtx`)，必须是存在的可执行文件。

## 参数说明 | Parameters { #parameters }

### 必需参数 | Required

**通俗理解|In plain words:** `-r` 参考基因组、`-i` GVCF 目录、`-o` 输出目录、`-g` GTX 可执行文件，四个都要给。

### 运行参数 | Runtime

**通俗理解|In plain words:** `-t` 是生成的每条命令里的线程数(默认 12)；`-T` 临时目录(默认 `./tmp`)；`-s` 输出脚本文件名(默认 `run_gtx_joint.sh`)；`-f` 是 faketime 伪装的时间点(默认 `2020-10-20 00:00:00`)。**一般不用动。**

### 拆分参数 | Splitting

**通俗理解|In plain words:** `-p` 是染色体过滤正则(如 `^Chr[0-9]+$` 只保留 Chr 开头带数字的)，用来只对部分染色体生成命令；`-w` 是区间窗口大小(bp)，指定后大样本会按这个窗口切成小块。**样本不超过 200 时用不到；样本很多、单条染色体跑不完时才用 `-w` 切分。**

## 分析流程 | Pipeline { #pipeline }

**通俗理解|In plain words:** 先确认环境和索引都齐了，再扫描 GVCF、读染色体，最后按样本量选择一种模式把命令写进脚本。

```text
输入 GVCF 目录 + 参考基因组 + GTX
    |
    ▼
步骤1: 验证配置 + 参考索引(samtools faidx) + faketime 检测
    |
    ▼
步骤2: 扫描/验证 GVCF(自动重建 tabix 索引)
    |
    ▼
步骤3: 读取染色体列表(.fai) + 应用 -p 过滤
    |
    ▼
步骤4: 备份旧脚本
    |
    ▼
步骤5: 按样本量生成命令(单机 / 按染色体 / 按区间)
    |
    ▼
步骤6: 设置脚本可执行权限
```

## 输出 | Output { #output }

```text
output/
├── run_gtx_joint.sh            # 生成的命令脚本(可执行)
└── gtx_joint.log               # 运行日志
```

运行脚本后，`output/` 下会多出结果(取决于模式)：

```text
单机模式:  gtx_joint_raw.vcf.gz
按染色体:  Chr01.joint.vcf.gz  Chr02.joint.vcf.gz ...
按区间:    Chr01_1-10000000.joint.vcf.gz  Chr01_10000001-20000000.joint.vcf.gz ...
```

## 结果解读 | Interpreting Results { #interpreting }

### 1. run_gtx_joint.sh(生成的脚本)

**通俗理解|In plain words:** 这是核心产出，里面一行一条 `gtx joint` 命令。生成的命令会带上 faketime、`-r` 参考、`-o` 输出、`-t` 线程，以及所有 GVCF 的 `-v` 参数。

- 单机模式(≤200 样本)：脚本里只有 1 条命令，直接 `bash run_gtx_joint.sh`
- 批处理模式(>200 样本)：脚本里每条染色体一条命令，可 `bash run_gtx_joint.sh` 串行，或 `parallel -j 4 < run_gtx_joint.sh` 并行

### 2. 区间拆分结果的合并

**通俗理解|In plain words:** 按区间拆分时，一条染色体被切成多块、产出多个 VCF，用完要拼回去。

```bash
bcftools concat -o Chr01.merged.vcf.gz output/Chr01_*.joint.vcf.gz
```

## 参数选择建议 | Parameter Guidance { #guidance }

- **样本 ≤200**：不用管拆分参数，默认单机模式一条命令搞定
- **样本 >200、单机跑不动**：默认自动切按染色体拆分，用 `parallel` 并行
- **单条染色体还是太大**：加 `-w 10000000`(每 10Mb 一块)进一步切分，跑完记得 `bcftools concat` 合并
- **只想跑部分染色体**：用 `-p` 正则过滤

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--ref, -r` | 必填 |  | 参考基因组文件路径｜Reference genome file path |
| `--input, -i` | 必填 |  | GVCF文件目录｜GVCF files directory |
| `--output, -o` | 必填 | Path | 输出目录｜Output directory |
| `--threads, -t` | `12` | int | 线程数｜Number of threads |
| `--tmp-dir, -T` | `./tmp` | Path | 临时目录｜Temporary directory |
| `--script, -s` | `run_gtx_joint.sh` |  | 输出脚本文件名｜Output script filename |
| `--faketime, -f` | `2020-10-20 00:00:00` |  | faketime时间｜faketime time |
| `--pattern, -p` | — |  | 染色体过滤正则表达式｜Chromosome filter pattern (e.g., "^Chr[0-9]+$") |
| `--window, -w` | — | int | 区间大小(bp)｜Window size in bp (e.g., 10000000 for 10M) |
| `--gtx, -g` | `~/software/gtx/bin/gtx` |  | GTX可执行文件路径｜GTX executable path |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-g, --gtx` | `~/software/gtx/bin/gtx` |  | GTX可执行文件路径｜GTX executable path |
| `-r, --ref` | 必填 |  | 参考基因组文件路径｜Reference genome file path |
| `-i, --input` | 必填 |  | GVCF文件所在目录｜GVCF files directory |
| `-o, --output` | 必填 |  | 输出结果目录｜Output directory |
| `-t, --threads` | `88` | int | 线程数｜Number of threads |
| `-T, --tmp-dir` | `./tmp` |  | 临时目录｜Temporary directory |
| `-s, --script` | `run_gtx_joint.sh` |  | 输出脚本文件名｜Output script filename |
| `-f, --faketime` | `2020-10-20 00:00:00` |  | faketime时间｜faketime time |
| `-p, --pattern` | — |  | 染色体过滤正则｜Chromosome filter pattern |
| `-w, --window` | — | int | 区间大小(bp)｜Window size in bp (e.g., 10000000 for 10M) |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- GTX(生成的脚本会调用它，`-g` 指定)
- samtools(参考基因组索引构建，自动解析 align 域环境并经 conda run 调用，可用环境变量 SAMTOOLS_PATH 覆盖；域环境缺失时回退 PATH 直接调用)
- tabix(GVCF 索引重建，自动解析 align 域环境并经 conda run 调用，可用环境变量 TABIX_PATH 覆盖；域环境缺失时回退 PATH 直接调用)
- faketime(可选，检测到才使用；GTX 的 license 校验依赖它)

## 常见问题 | FAQ { #faq }

**Q1：为什么它不直接跑，只生成脚本？**
这是它的设计定位——把「生成可并行、可拆分的联合命令」和「真正跑」分开，方便你拿到脚本后在计算节点提交、用 parallel 并行，或人工检查命令无误后再执行。

**Q2：样本数超过 200 会怎样？**
自动从单机模式切换为按染色体拆分；若同时给了 `-w`，则按区间拆分。

**Q3：GVCF 文件没找到？**
程序只认文件名匹配 `*.g.vcf.gz` 的文件，其它命名不会被扫描到。缺 `.tbi` 索引或索引过时会自动用 tabix 重建。

**Q4：faketime 没装会报错吗？**
不会。faketime 可用则用，不可用会在日志里 WARNING 并不加 faketime 前缀；但实际跑 GTX 时若 license 校验失败，仍建议安装 libfaketime。

**Q5：参考基因组没有 .fai 索引？**
程序会自动用 samtools faidx 生成，无需手动建。
