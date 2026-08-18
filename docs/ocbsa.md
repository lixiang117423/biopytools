# OcBSA 混池分离分析 | OcBSA BSA Analysis Suite

一句话理解：**一套 BSA（混池分离分析）工具箱，从 F1/F2 群体的 VCF 里定位控制性状的候选区域，并配套出图和候选区引物设计**，覆盖「定位 → 验证」的完整链路。

## 功能概述 | Overview

- 4 个子命令：`f1`（F1 群体 DHHP 分析）、`f2`（F2/RILs 群体 SNP-index/ED 分析）、`fig`（结果可视化）、`primer`（候选区域引物设计）
- `f1` 面向显性杂交 F1（outcross）群体，用 DHHP 算法（OcBSA，Molecular Plant 2024）
- `f2` 面向 F2/RILs 分离群体，支持 SNP-index 与 ED（欧氏距离）两种方法
- 亲本与混池各自独立的深度上下限过滤，按染色体多进程并行
- 滑窗平滑 + 自动绘图，免去手写 R 脚本

## 快速开始 | Quick Start

```bash
biopytools ocbsa f1 -i input.vcf -p1 1 -p2 2 -b1 3 -b2 4 -o ./output
```

`-p1/-p2` 是两个亲本在 VCF 样本区的列号，`-b1/-b2` 是两个混池的列号（均从 1 开始数）。

## 零基础概念速览 | Concepts in plain words

| 术语 | 通俗理解 |
|------|----------|
| BSA（混池分离分析） | 把性状相反的后代各自混成一个池测序，找两池差异最大的区域来定位基因 |
| 亲本/混池列号 | VCF 里每个样本占一列；`-p1` 等参数告诉程序「第几列是谁」 |
| F1（outcross） | 两个不同品种杂交的第一代，后代像「混合体」；DHHP 专为这种显性定位设计 |
| F2/RILs | 自交分离群体，个体趋近纯合；用 SNP-index 或 ED 找差异 |
| DHHP | 一种针对 F1 显性基因的统计量，相当于「两池支持/反对该位点的证据差」 |
| SNP-index / ED | 两种度量：SNP-index 是两池等位频率之差；ED 是两池频率向量的欧氏距离 |
| 滑窗平滑 | 把相邻位点的值取平均，抹平单点噪音看趋势 |

## 输入 | Input

### VCF 文件

标准 VCF（`.vcf` / `.vcf.gz`），FORMAT 列需含 AD（`GT:AD`），如 `0/1:12,30`。列号 `-p1/-p2/-b1/-b2` 从 1 开始数，指 VCF 样本区（`#CHROM` 行 FORMAT 之后）的第几个样本；程序会先打印检测到的样本名与列号映射，便于核对。

### 其他子命令的输入

- `fig`：`-i` 输入 `f1/f2` 生成的平滑结果文件（`*.smoothed`），`-o` 输出 `.png/.pdf`
- `primer`：`-g` 参考基因组 FASTA + `-i` OcValue 文件 + `--region chr,start,end`（逗号分隔）

## 分析流程 | Pipeline

```text
输入 VCF
    │
    ▼
读取 VCF(检测样本列号映射)
    │
    ▼
深度/基因型过滤(亲本、混池各用独立上下限)
    │
    ▼
计算指标(f1: DHHP / f2: SNP-index 或 ED)
    │
    ▼
滑窗平滑(按染色体多进程并行)
    │
    ▼
输出结果文件 + 自动绘图
```

## 输出 | Output

### f1 子命令（`ocbsa f1`）

```text
output/
├── ocbsa.OcValue          # 逐位点 DHHP 原始结果(核心)
├── ocbsa.smoothed         # 滑窗平滑后的 OcValue
├── ocbsa.summary.tsv      # 供绘图的汇总数据
├── ocbsa.OcValue.png      # 自动绘制的 OcValue 图
└── 99_logs/ocbsa_f1.log   # 运行日志
```

### f2 子命令（`ocbsa f2`）

```text
output/
├── f2bsa.snpindex 或 f2bsa.ED    # 逐位点结果(--method 决定)
├── f2bsa.smoothed                # 滑窗平滑结果
├── f2bsa.summary.tsv             # 供绘图的汇总数据
├── f2bsa.SNP-index.png 或 f2bsa.ED.png  # 自动绘制的图
└── 99_logs/ocbsa_f2.log
```

### primer 子命令（`ocbsa primer`）

```text
output/
├── primer_design_results.txt      # 引物设计结果(最终产物)
├── {chr}_{start}_{end}_indel.fasta  # 提取的 INDEL 侧翼序列
├── indel_genome.blast             # BLAST 比对结果
└── db/                            # BLAST 数据库
```

## 结果解读 | Interpreting Results

- **`*.OcValue / *.snpindex`**：逐位点的统计值，绝对值越大越像与性状连锁
- **`*.smoothed`**：滑窗平滑后的值，用于画图和找峰；图上明显隆起、超过阈值线的山头就是候选区域
- **`*.png 里的阈值线**：f1 的 OcValue 图在非指定区域模式下画红/绿/蓝三条虚线，分别对应 top95/99/99.9 分位数，点越靠近红色以上越显著
- **`primer_design_results.txt`**：每行一对引物（左右序列、Tm、GC、产物长度），用于候选区验证
- **好坏判据**：候选峰若由多个相邻位点共同支持（不是单点尖峰）、且落在目标性状已知区域，结论更可信

## 参数选择建议 | Parameter Guidance

**通俗理解|In plain words:** 亲本/混池深度阈值和窗口大小一般用默认即可；真正需要定的是「你的群体是 F1 还是 F2」以及四个列号。

- **`f1 还是 f2`**：显性杂交 F1（outcross）用 `f1`；自交分离群体（F2/RILs）用 `f2`
- **`f2 的 --method`**：默认 `snpindex` 最常用；`ED` 对某些数据更稳健，两者可都跑一遍互相印证
- **`-w/--window-size`（滑窗）**：默认 1Mb 一般不用动；标记稀疏可加大，标记极密可减小
- **深度阈值（--parent-min/max-dep、--pool-min/max-dep）**：默认亲本 10–100、混池 20–500，一般不用动；只在测序深度明显偏离时才调
- **`primer` 的引物参数**：默认（长度 18–24、Tm 50–65、产物 70–200）即可，一般不用动

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-vcf` | 必填 |  | VCF文件路径｜VCF file path |
| `-p1` | 必填 | int | 显性亲本列号｜Dominant parent column |
| `-p2` | 必填 | int | 隐性亲本列号｜Recessive parent column |
| `-b1` | 必填 | int | 显性表型混池列号｜Dominant pool column |
| `-b2` | 必填 | int | 隐性表型混池列号｜Recessive pool column |
| `-o, --output-dir` | `./output` |  | 输出目录｜Output directory |
| `-w, --window-size` | `1000000` | int | 滑窗大小｜Window size |
| `-p, --pvalue` | `99` | float | p值阈值｜P-value threshold |
| `--parent-min-dep` | `10` | int | 亲本最低覆盖度｜Parent min depth |
| `--parent-max-dep` | `100` | int | 亲本最高覆盖度｜Parent max depth |
| `--pool-min-dep` | `20` | int | 混池最低覆盖度｜Pool min depth |
| `--pool-max-dep` | `500` | int | 混池最高覆盖度｜Pool max depth |
| `--method` | `snpindex` | snpindex/ED | 分析方法｜Analysis method |
| `-i, --input-file` | 必填 |  | 输入结果文件｜Input result file |
| `-o, --output-file` | 必填 |  | 输出图片路径(.png/.pdf)｜Output figure path |
| `--plot-type` | `ocvalue` | ocvalue/snpindex/ed | 图表类型｜Plot type |
| `--position` | — |  | 特定区域(chr,start,end)｜Specific region |
| `--color` | `plasma_r` |  | 颜色方案｜Color scheme |
| `-g, --genome` | 必填 |  | 参考基因组路径｜Reference genome path |
| `-i, --ocvalue-file` | 必填 |  | OcValue文件路径｜OcValue file path |
| `--region` | 必填 |  | 目标区间(chr,start,end)｜Target region |
| `-n, --primer-num` | `10` | int | 引物对数量｜Number of primer pairs |
| `--flank-length` | `200` | int | INDEL侧翼长度｜INDEL flank length |
| `--primer-min-size` | `18` | int | 最短引物｜Min primer size |
| `--primer-opt-size` | `20` | int | 最适引物｜Opt primer size |
| `--primer-max-size` | `24` | int | 最长引物｜Max primer size |
| `--product-min` | `70` | int | 最短产物｜Min product size |
| `--product-max` | `200` | int | 最长产物｜Max product size |
| `--min-tm` | `50.0` | float | 最低Tm｜Min Tm |
| `--max-tm` | `65.0` | float | 最高Tm｜Max Tm |
| `--min-gc` | `35.0` | float | 最低GC｜Min GC% |
| `--max-gc` | `65.0` | float | 最高GC｜Max GC% |
| `--tm-diff` | `0.5` | float | Tm差异阈值｜Tm diff threshold |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-vcf` | 必填 |  | VCF文件路径｜VCF file path |
| `-p1, --p1` | 必填 | int | 显性亲本列号｜Dominant parent column |
| `-p2, --p2` | 必填 | int | 隐性亲本列号｜Recessive parent column |
| `-b1, --b1` | 必填 | int | 显性表型混池列号｜Dominant pool column |
| `-b2, --b2` | 必填 | int | 隐性表型混池列号｜Recessive pool column |
| `-o, --output-dir` | `./output` |  | 输出目录｜Output directory |
| `-w, --window-size` | `1000000` | int | 滑窗大小｜Window size |
| `-p, --pvalue` | `99` | float | p值阈值｜P-value threshold |
| `--parent-min-dep` | `10` | int | 亲本最低覆盖度｜Parent min depth |
| `--parent-max-dep` | `100` | int | 亲本最高覆盖度｜Parent max depth |
| `--pool-min-dep` | `20` | int | 混池最低覆盖度｜Pool min depth |
| `--pool-max-dep` | `500` | int | 混池最高覆盖度｜Pool max depth |
| `--method` | `snpindex` | snpindex/ED | 分析方法｜Analysis method |
| `-i, --input-file` | 必填 |  | 输入结果文件｜Input result file |
| `-o, --output-file` | 必填 |  | 输出图片路径｜Output figure path |
| `--plot-type` | `ocvalue` | ocvalue/snpindex/ed | 图表类型｜Plot type |
| `--position` | — |  | 特定区域(chr,start,end)｜Specific region |
| `--color` | `plasma_r` |  | 颜色方案｜Color scheme |
| `-g, --genome` | 必填 |  | 参考基因组路径｜Reference genome path |
| `-i, --ocvalue-file` | 必填 |  | OcValue文件路径｜OcValue file path |
| `--region` | 必填 |  | 目标区间(chr,start,end)｜Target region |
| `-n, --primer-num` | `10` | int | 引物对数量｜Number of primer pairs |
| `--flank-length` | `200` | int | INDEL侧翼长度｜INDEL flank length |
| `--primer-min-size` | `18` | int | 最短引物｜Min primer size |
| `--primer-opt-size` | `20` | int | 最适引物｜Opt primer size |
| `--primer-max-size` | `24` | int | 最长引物｜Max primer size |
| `--product-min` | `70` | int | 最短产物｜Min product size |
| `--product-max` | `200` | int | 最长产物｜Max product size |
| `--min-tm` | `50.0` | float | 最低Tm｜Min Tm |
| `--max-tm` | `65.0` | float | 最高Tm｜Max Tm |
| `--min-gc` | `35.0` | float | 最低GC｜Min GC% |
| `--max-gc` | `65.0` | float | 最高GC｜Max GC% |
| `--tm-diff` | `0.5` | float | Tm差异阈值｜Tm diff threshold |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 库：numpy、matplotlib（绘图）
- BLAST（`makeblastdb` / `blastn`，仅 `primer` 子命令需要；自动解析 annot 域环境并经 conda run 调用，可用环境变量 MAKEBLASTDB_PATH / BLASTN_PATH 覆盖；域环境缺失时回退 PATH 直接调用）
- primer3-py（仅 `primer` 子命令需要）

## 常见问题 | FAQ

**Q1：`-p1/-p2/-b1/-b2` 的列号怎么数？**
从 1 开始，指 VCF 样本区（`#CHROM` 行 FORMAT 字段之后）的第几个样本。程序启动会打印「检测到 N 个样本」和每个列号对应的样本名，照此核对即可；超出样本总数会直接报错。

**Q2：f1 和 f2 该用哪个？**
取决于群体类型：F1（两个品种杂交的后代，尚未自交）用 `f1`；F2/RILs（已自交分离、个体趋近纯合）用 `f2`。选错会得到无意义的结果。

**Q3：primer 子命令报「makeblastdb/blastn 找不到」？**
引物设计依赖 BLAST 去重（过滤非唯一序列）。请确认 `makeblastdb`、`blastn` 在 PATH 中，且已安装 primer3-py（`pip install primer3-py`）。

**Q4：为什么有的染色体没出现在结果里？**
f1 分析会跳过标记数 <1000 的染色体（数据量不足以可靠计算），日志会提示。
