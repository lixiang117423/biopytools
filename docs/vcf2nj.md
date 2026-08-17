# vcf2nj - VCF 邻接法建树 | vcf2nj Neighbor-Joining Tree from VCF

一句话理解：**从 VCF(SNP) 出发，先算两两样本的遗传距离矩阵(VCF2Dis)，再用「邻接法」(NJ) 建一棵进化树，可选外群重根。**

## 功能概述 | Overview

- VCF -> 距离矩阵(VCF2Dis) -> 邻接法(NJ)建树 -> 可选外群重根，三步一条龙
- 可 `--skip-vcf2dis` 跳过距离计算，直接用已有距离矩阵建树
- 输出 Newick 格式进化树 + 距离矩阵 + 总结报告
- 支持逗号分隔多个外群做重根化

## 快速开始 | Quick Start

```bash
biopytools vcf2nj -i wild.snp.vcf -o output_dir
```

最小输入：一个 VCF 文件(`-i`)和一个输出目录(`-o`)。结果文件(距离矩阵、树、日志)写在输出目录里。

## 零基础概念速览 | Concepts in plain words

| 术语<br>Term | 通俗理解<br>In plain words |
|------|------|
| VCF | 记录「每个样本在每个位点是什么基因型」的标准变异文件 |
| 距离矩阵 | 一张「两两样本有多不同」的对称表，数字越小亲缘越近 |
| 邻接法(NJ) | 一种「由距离反推树」的快速算法：谁俩离得近就先聚在一起 |
| Newick | 用括号串 `(A,(B,C));` 表示进化树的标准格式 |
| 外群(outgroup) | 明确「最远亲」的样本，用来给树定方向(重根) |

## 输入 | Input

### VCF 文件

标准 VCF 格式，需含 `#CHROM` 样本头行(前 9 列固定，之后每列一个样本)。示例：

```text
##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1  sample2  sample3
Chr1  100  .  A  G  .  PASS  .  GT  0/0  0/1  1/1
```

### 距离矩阵文件(可选)

用 `--skip-vcf2dis -d <矩阵文件>` 时，可跳过 VCF2Dis 直接给一个现成距离矩阵(首行为样本数，之后每行「样本名 + 距离值」)。

## 参数说明 | Parameters

### 输入 | Input

**通俗理解|In plain words:** `-i/--input` 是输入 VCF；`-d/--distance-matrix` 是「已有的距离矩阵」——只想跳过距离计算、复用现成矩阵时才用(需配合 `--skip-vcf2dis`)。

### 输出 | Output

**通俗理解|In plain words:** `-o/--output` 是输出目录(默认当前目录 `.`)。`-p/--prefix` 是输出文件名前缀，距离矩阵叫 `<前缀>.dis`、树叫 `<前缀>.nwk`。`-t/--tree-output` 可单独指定树的输出路径。

### 重根化 | Rerooting

**通俗理解|In plain words:** `--outgroup` 指定外群样本(逗号分隔多个)，程序会用 nw_reroot 把树「重新定根」到外群上，另存一个重根树。**不给外群就不重根。** `--nw-reroot-path` 是 nw_reroot 程序路径(默认已配好)，一般不用动。

### 工具与行为 | Tool paths and behavior

**通俗理解|In plain words:** `--vcf2dis-path` 是 VCF2Dis 程序路径(默认 `VCF2Dis`，需在 PATH 中)。`-w/--working-dir` 是执行外部命令的工作目录。`--skip-vcf2dis` 跳过距离计算步骤(此时必须给 `-d`)。

## 分析流程 | Pipeline

```text
VCF 文件
    │
    ▼
步骤1: VCF2Dis 计算两两遗传距离矩阵(<前缀>.dis)
    │
    ▼
步骤2: scikit-bio 用邻接法(NJ)建树(<前缀>.nwk)
    │
    ▼
步骤3: 验证结果 + 生成总结报告(<前缀>_summary.txt)
    │
    ▼
步骤4: (可选) 指定外群时用 nw_reroot 重根(<base>_<外群>_root.nwk)
```

## 输出 | Output

```text
output_dir/
├── <prefix>.dis                  # 两两遗传距离矩阵
├── <prefix>.nwk                  # NJ 进化树(Newick)
├── <prefix>.log                  # 运行日志
├── <prefix>_summary.txt          # 总结报告(写在当前工作目录)
└── <base>_<外群>_root.nwk        # 重根后的树(仅指定 --outgroup 时)
```

其中 `<prefix>` 是 `-p` 指定的前缀，`<base>` 是树文件 basename。

## 结果解读 | Interpreting Results

### 1. 进化树(<prefix>.nwk)

**通俗理解|In plain words:** 最终结果，Newick 格式。用 FigTree、iTOL 打开即见树形图。注意 NJ 树的分支长度反映「遗传距离」，**分叉越靠外的样本差异越大**。

### 2. 距离矩阵(<prefix>.dis)

**通俗理解|In plain words:** 两两样本的遗传距离表，对角线为 0。可以直接用来做聚类热图，或喂给其它建树软件。数值越大=两样本差异越大。

### 3. 重根树(<base>_<外群>_root.nwk)

**通俗理解|In plain words:** 只有指定 `--outgroup` 才会生成。它把树「重新定根」到外群上，让进化方向更明确(外群在最根部)。

## 参数选择建议 | Parameter Guidance

- **常规流程**：直接 `-i VCF -o 输出目录 -p 前缀` 即可
- **想复用已有距离矩阵**：`--skip-vcf2dis -d <矩阵文件>`
- **想给树定方向**：加 `--outgroup <外群样本名>`(多个用逗号分隔)
- **VCF2Dis 不在 PATH**：用 `--vcf2dis-path` 指定完整路径

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | — | Path | 输入VCF文件路径｜Input VCF file path |
| `--distance-matrix, -d` | — | Path | 已有距离矩阵文件路径｜Existing distance matrix file path |
| `--output, -o` | `.` | str | 输出目录｜Output directory |
| `--prefix, -p` | `vcf2nj` | str | 输出文件前缀｜Output file prefix |
| `--tree-output, -t` | — | Path | 系统发育树输出文件路径｜Phylogenetic tree output file path |
| `--outgroup` | — | str | 外群样本标签，多个用逗号分隔｜Outgroup sample labels, comma-separated |
| `--nw-reroot-path` | `~/miniforge3/envs/phylo/bin/nw_reroot` | str | nw_reroot程序路径｜nw_reroot program path |
| `--vcf2dis-path` | `VCF2Dis` | str | VCF2Dis程序路径｜VCF2Dis program path |
| `--working-dir, -w` | `.` | Path | 工作目录｜Working directory |
| `--skip-vcf2dis` | — |  | 跳过VCF2Dis步骤｜Skip VCF2Dis step |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | — |  | 输入VCF文件路径｜Input VCF file path |
| `-d, --distance-matrix` | — |  | 已有的距离矩阵文件路径（用于跳过VCF2Dis步骤）｜ Existing distance matrix file path (for skipping VCF2Dis step) |
| `-o, --output` | `.` |  | 输出目录｜Output directory |
| `-p, --prefix` | `phylo_analysis` |  | 输出文件前缀｜Output file prefix |
| `-t, --tree-output` | — |  | 系统发育树输出文件路径｜Phylogenetic tree output file path |
| `--outgroup` | — |  | 外群样本标签，多个用逗号分隔｜Outgroup sample labels, comma-separated |
| `--nw-reroot-path` | `~/miniforge3/envs/phylo/bin/nw_reroot` |  | nw_reroot程序路径｜nw_reroot program path |
| `--vcf2dis-path` | `VCF2Dis` |  | VCF2Dis程序路径｜VCF2Dis program path |
| `-w, --working-dir` | `.` |  | 工作目录｜Working directory |
| `--skip-vcf2dis` | — | store_true | 跳过VCF2Dis步骤，直接从距离矩阵构建树｜Skip VCF2Dis step, build tree directly from distance matrix |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

| 软件/库<br>Software/Library | 说明<br>Description |
|------|------|
| VCF2Dis | 计算遗传距离，默认程序名 `VCF2Dis`(需在 PATH)，可用 `--vcf2dis-path` 指定 |
| nw_reroot | 重根化工具(仅指定外群时需要)，默认 `~/miniforge3/envs/phylo/bin/nw_reroot` |
| numpy / pandas / scipy | Python 数值与表格处理库 |
| scikit-bio | 提供邻接法(NJ)建树(导入名 `skbio`) |

## 常见问题 | FAQ

### 1. 会断点续传吗？

**不会。** 每次运行都会重新执行 VCF2Dis 距离计算和 NJ 建树。中断后重跑会从头再来。

### 2. 报「缺少 Python 依赖 scikit-bio」？

NJ 建树依赖 scikit-bio(导入名 `skbio`)。请安装：`pip install scikit-bio` 或 `conda install -c conda-forge scikit-bio`。

### 3. 距离矩阵 / 树 / 日志各写在哪个目录？

距离矩阵(`<前缀>.dis`)、树(`<前缀>.nwk`)、日志(`<前缀>.log`)写在 `-o` 指定的输出目录；总结报告(`<前缀>_summary.txt`)写在**当前工作目录**。建议显式指定 `-o` 和 `-p` 避免混淆。

### 4. 没有指定前缀时输出文件叫什么？

注意：CLI 帮助里 `-p` 默认值显示为 `vcf2nj`，但经 `biopytools vcf2nj` 调用、不显式传 `-p` 时，实际生效的是模块默认前缀 `phylo_analysis`。**建议始终显式指定 `-p`** 以免输出文件名与预期不符。

### 5. VCF 必须只含 SNP 吗？

VCF2Dis 主要面向 SNP 数据。建议先用质控工具把 VCF 过滤干净(去低质量、去多等位)再传入，距离矩阵和树会更可靠。
