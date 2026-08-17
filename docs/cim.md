# CIM - 复合区间作图分析 | Composite Interval Mapping (R/qtl)

一句话理解：**把「决定某个性状的基因在基因组哪个位置」这件事，用一套统计扫描自动找出来**。
输入一个群体的基因型数据(VCF)和性状记录(表型文件)，输出一张「嫌疑区域排行榜」(LOD 曲线 + 峰值表)。

## 功能概述 | Overview

- 支持三种图谱模式：`physical`（物理位置）、`estimate`（R/qtl `est.map()` 估算）、`mstmap`（MSTmap 构建连锁图谱，默认）
- **mstmap 模式下物理图谱与连锁图谱各跑一次 CIM**，互相补充验证
- 内置重组频率(RF)质控：杂合比例、局部平均 RF、短距离重组热点过滤（双亲群体大片交换假设）
- Pre-RF（基线）与 Post-RF（质控后）两轮分析，方便对比过滤效果
- 置换检验自动计算显著性阈值，LOD 扫描数据与 QTL 峰值表可直接供外部绘图
- 断点续传：已完成步骤自动跳过（换参数重跑需先删旧产物，见 FAQ）

## 快速开始 | Quick Start

```bash
biopytools cim -i input.vcf.gz -p phe.txt -o output_dir
```

最小输入：一个 VCF（.vcf/.vcf.gz）+ 一个两列表型文件（TSV，制表符分隔的 sample、value 两列）。

## 零基础概念速览 | Concepts in plain words

不熟悉生信术语的话，先花两分钟看这张表，后面的参数说明都会用到：

| 术语 | 通俗理解 |
|------|----------|
| 标记/SNP | 基因组上的一个「点位」，像高速公路上的里程桩，本工具扫描的就是这些桩 |
| 基因型 | 一个位点上两条同源染色体的构成；工具内部记为 A(父本纯合)/H(杂合)/B(母本纯合) |
| MAF | 「少数派基因型」在群体里占多少；太低=这个位点区分不了两个亲本，没信息量 |
| 缺失率 | 多少样本在这个位点「没测出来」；缺失太多=判断依据不足 |
| 杂合率 | 一个位点上有多少样本是「双亲各出一份」(H)；几乎全 H 多半是数据有问题 |
| 重组频率(RF) | 两个位点之间发生「交换」的比例；离得近的位点应当「一起走」，频繁分开就是可疑 |
| 连锁(LD) | 相邻位点往往「绑定遗传」，信息重复；降维=每个「团」留一个代表 |
| 物理位置(bp) | 位点在染色体上的实际坐标，像「门牌号」 |
| 遗传距离(cM) | 按「重组频率」折算的距离，像「开车时间」：物理距离一样，堵车路段(重组冷点)跑得慢 |
| LOD 得分 | 某个位置「存在 QTL」的证据强度打分，分越高越像真的 |
| QTL | 控制性状的基因组区域，本工具要找的东西 |
| 置换检验阈值 | 把表型标签随机打乱几千遍，看「纯靠运气」最高能拿多少分；超过运气上限才算真信号 |

## 输入 | Input

### VCF 文件

标准 VCF 格式（支持 `.vcf` / `.vcf.gz`），基因型编码 `0/0`、`0/1`、`1/1`，内部自动转为 ABH 编码（A=纯合父本，H=杂合，B=纯合母本）。

### 表型文件

TSV 格式，两列：`sample`（样本名，须与 VCF 一致）和性状值：

```text
sample    trait1
21-18     0
21-19     1
...
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 说明 |
|------|------|
| `-i, --input` | 输入 VCF 文件（.vcf / .vcf.gz） |
| `-p, --pheno` | 表型文件（TSV: sample, value） |
| `-o, --output` | 输出目录 |

### 群体类型与图谱模式 | Cross type & map mode

**通俗理解|In plain words:** 告诉程序两件事：你的群体是哪种杂交后代（这决定统计模型），以及「给位点排队」用什么尺子（实际染色体坐标，还是按重组率折算的遗传距离）。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-t, --type` | f2 | 群体类型：f2（F2 群体）/ bc（回交群体），按实验设计选 |
| `--map-mode` | mstmap | cM 来源：physical / estimate / mstmap |

- **physical**：直接用基因组物理位置（单位 Mb，pseudo-cM），不构建图谱，最快
- **estimate**：R/qtl `est.map()` 估算遗传距离
- **mstmap**：MSTmap 构建连锁图谱，physical 与 mstmap 两次 CIM 都跑

### 标记过滤 | Marker filtering

**通俗理解|In plain words:** 这两个参数决定「什么样的位点信不过、直接扔掉」。MAF 太低=这个位点在群体里几乎没有差异，区分不了它来自哪个亲本，留着只会添乱；缺失太多=很多样本在这个位点没数据，判断依据不足。**绝大多数项目用默认值即可，几乎不需要动。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--maf` | 0.05 | 最小等位基因频率(MAF)阈值，低于此值的位点删除 |
| `--missing` | 0.1 | 最大缺失率，缺失比例高于此值的位点删除 |

### 重组频率质控 | RF quality control

**通俗理解|In plain words:** 双亲的染色体是整段传给后代的，正常位点之间要么「一起走」（重组率低），要么远到随机（≈0.5）。如果一个位点和邻居的「分手率」异常高，大概率是它的基因型被测序/比对弄错了（典型原因：旁系同源——两条相似的染色体被当成一个位点）。这一步就是把这些「撒谎的位点」揪出来。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--max-het-rate` | 0.6 | 杂合基因型(H)最大比例，超过则删除该标记 |
| `--max-mean-rf` | 0.35 | K 近邻局部平均重组频率(RF)最大值，超过则删除 |
| `--rf-knn` | 10 | 局部 RF 窗口邻居数（物理位置前后各一半） |

### 短距离重组热点过滤 | Local recombination hotspot filter

**通俗理解|In plain words:** 两个位点物理上挨得极近（<1000bp）时，理论上几乎总是「同进同退」。如果它们频繁「分开」，现实中几乎不可能真发生，只能是数据错误。这一步在最近距离上专门抓这种错，比上一组更严格。**默认开启；正常数据被删掉的很少，若删得特别多，说明数据本身问题大。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--local-hotspot/--no-local-hotspot` | 开 | 总开关 |
| `--local-hotspot-dist` | 1000 | 物理距离阈值(bp)，距离更近的相邻对才检查 |
| `--local-hotspot-rf` | 0.20 | 相邻 RF 软阈值，达到即给两侧标记各记 1 分 |
| `--local-hotspot-hard-rf` | 0.30 | 硬阈值，单对 RF 达到即两侧标记都删 |
| `--local-hotspot-score` | 1 | 标记级删除线，软评分≥此值删除 |
| `--local-hotspot-relative` | 0.0 | 相对判据系数；>0 时阈值=max(绝对值, 系数×近距对 RF 中位数)，适合高重组背景染色体 |

诊断输出：`01_qc/local_hotspot_removed.tsv`（被删标记清单：marker_id/chr/pos/reason/score/max_adj_rf）。

### 物理距离降采样 | Physical thinning

**通俗理解|In plain words:** 几十万上百万个位点里，很多挨得极近、信息几乎一模一样（像同一段话被抄了很多遍）。只保留每个「密集区」的第一个代表，计算量骤降、结果几乎不变。位点特别多时才需要开。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--min-phys-gap` | 0 | 最小相邻物理距离(bp)，0=关闭；相邻距离小于此值的密集簇仅保留簇头 |

在热点过滤之后执行。66 万 SNP 配 `--min-phys-gap 500` 约剩 15–25 万标记，对 ~1cM 的 QTL 分辨率无损，LD/MSTmap/CIM 耗时与内存显著下降。

### LD 降维 | LD pruning

**通俗理解|In plain words:** 给长文章提取摘要——相邻位点往往「绑定遗传」，信息重复，每个「团」留一个代表就够。**默认关闭**，因为大家的输入数据通常已经做过这一步。

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--skip-ld` | true | **默认跳过 LD 降维**（输入数据通常已预降维） |
| `--no-skip-ld` | - | 启用 LD 降维 |
| `--ld-window` | 50 | LD 计算窗口（SNP 数） |
| `--ld-step` | 5 | LD 计算步长（SNP 数） |
| `--ld-r2` | 0.1 | LD r² 阈值 |

### CIM 扫描参数 | CIM scan

**通俗理解|In plain words:** 这是「扫描」本身的设置：window 是每个位置算分时参考多大邻域，step 是每隔多长距离算一次（越细曲线越平滑、越慢），协因子(n-marcovar)用来扣除其他区域的干扰让每个位置的信号更干净。**默认值经实践验证，一般不动。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--n-marcovar` | 10 | 协因子数量 |
| `--window` | 10.0 | 窗口大小（图谱单位：mstmap/estimate 为 cM，physical 为 Mb/pseudo-cM） |
| `--method` | hk | 扫描方法：hk / em / imp |
| `--step` | 1.0 | 伪标记步长（图谱单位） |

> 注意：物理模式坐标是 Mb（pseudo-cM），LOD 与阈值**只在同一图谱内比较**（见 FAQ）。

### 置换检验 | Permutation test

**通俗理解|In plain words:** 用来回答「这个峰是真的还是运气」——把表型标签随机打乱几千次重算，看纯运气最高能拿多少分，这就是阈值。**1000 次是常用精度；想快速预览可以先设 0 跳过（此时没有阈值，只有曲线形状可看）。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--n-perm` | 1000 | 置换检验次数（0=跳过，不生成显著性阈值） |

### MSTmap 参数（仅 mstmap 模式）| MSTmap options

**通俗理解|In plain words:** 用「重组频率」给位点重新排队、算出遗传距离(cM)，并把位点自动分成一组组「连锁群」。p 值决定分组松紧，**程序会自动反复调整到合适值，一般不用管。**

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--mstmap-pvalue` | 1e-6 | 聚类 p 值起始值，**自动调优**（目标核心 LG 数 ≤ 染色体数×3） |
| `--mstmap-distfun` | kosambi | 距离函数：kosambi / haldane |

### 环境参数 | Environment

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--r-env` | ~/miniforge3/envs/Rqtl | R conda 环境路径或名称 |

> 以下参数仅在直接调用 `python -m biopytools.cim` 时可用（click 包装器未暴露，用默认值）：`--mstmap-path`（默认 `~/miniforge3/envs/r/bin/mstmap`）、`--threads`（默认 1）。

## 分析流程 | Pipeline

**通俗理解|In plain words:** 先不过滤跑一遍(Pre-RF，当作「体检前」的基线)，再用重组质控把可疑位点删掉跑一遍(Post-RF，最终结果)，两轮对比着看过滤是否合理。

```text
输入 VCF + 表型文件
    │
    ▼
步骤1: 解析输入(VCF → ABH 编码矩阵)
    │
    ▼
步骤2: 标记过滤(MAF / 缺失率)
    │
    ▼
步骤2.5a: 短距离重组热点过滤(可选)
步骤2.5b: 物理距离降采样(可选)
    │
    ▼
步骤3: LD 降维(可选,默认跳过)
    │
    ▼
Pre-RF CIM 分析(基线,不做 RF 过滤)
    ├─ 构建 csvs 格式文件
    ├─ [mstmap] MSTmap 构建连锁图谱
    ├─ [mstmap] 物理位置 CIM + 连锁图谱 CIM
    ├─ [estimate] 估算遗传距离后 CIM
    └─ [physical] 直接用物理位置 CIM
    │
    ▼
步骤5: 重组频率(RF)质控
    ├─ H 比例过滤(杂合比例 > max_het_rate → 删)
    ├─ 平均 RF 过滤(局部平均 RF > max_mean_rf → 删)
    └─ 孤立重组检测(与相邻标记 RF 均 > 0.5 → 仅警告)
    │
    ▼
Post-RF CIM 分析(最终结果,同 Pre-RF 流程,用 RF 过滤后的标记)
    │
    ▼
输出结果
```

### 重组频率(RF)定义 | RF definition

双纯合口径（2026-08-15 修复）：只统计双纯合样本对，含 H 的样本对排除出分子和分母。旧「含 H 记 0.5」会把完全连锁标记对的 RF 基线抬到 ~0.25-0.3（F2 H 比例 50% 时），使所有阈值失效。

| 两个样本的基因型 | 判定 | 信号值 |
|---|---|---|
| A-A 或 B-B | 不重组 | 0 |
| A-B 或 B-A | 重组 | 1 |
| 含 H | 无确定信息 | 排除（分子分母均不计） |
| 含缺失 | 跳过 | - |

```text
RF(标记i, 标记j) = 双纯合样本对中的重组数 / 双纯合样本对数
局部平均RF(标记i) = 与物理位置前后各 rf_knn//2 个标记的RF之和 / 有效配对数
```

规则 2 用**局部**平均 RF 而非全染色体平均：全染色体平均在正确 RF 定义下无区分力（远距离对 RF≈0.5 占主导，正常标记的全局均值也 ≈0.45-0.5）。

## 输出 | Output

```text
output_dir/
├── 00_pipeline_info/
│   └── pipeline_params.txt              # 流程参数记录
├── 01_qc/
│   ├── marker_filter_stats.txt          # MAF/缺失率过滤统计
│   ├── hotspot_filter_stats.txt         # 短距离重组热点过滤统计
│   ├── phys_thin_stats.txt              # 物理降采样统计
│   ├── ld_prune_stats.txt               # LD 降维统计(启用时)
│   ├── rf_filter_stats.txt              # RF 质控过滤统计
│   ├── singleton_het_report.tsv         # 孤立重组标记 H 比例诊断
│   ├── local_hotspot_removed.tsv        # 热点过滤被删标记清单
│   ├── filtered_markers.vcf.gz          # 过滤/降维后保留的标记 VCF
│   └── rf_filtered_markers.vcf.gz       # RF 过滤后保留的标记 VCF
├── 02_cim/
│   ├── pre_rf/                          # Pre-RF 分析结果(基线)
│   │   ├── tidy_files/                  # csvs 格式输入文件
│   │   ├── physical/                    # 物理位置 CIM 结果(pdf/rds)
│   │   ├── mstmap/                      # 连锁图谱 CIM 结果(仅 mstmap 模式)
│   │   │   ├── linkage_map.csv          # 连锁图谱(标记→LG→cM)
│   │   │   ├── mstmap_map.csv / mstmap_gen.csv / mstmap_phe.csv
│   │   │   │                            #   (仅含每染色体最大LG,碎片LG已剔除;
│   │   │   │                            #    完整聚类结果在 mstmap_pvalue_*/)
│   │   │   └── marker_map_index.tsv     # 标记物理位置↔图谱位置索引
│   │   ├── plots/                       # 可视化数据(供外部绘图)
│   │   │   ├── cim_lod_data_physical.tsv    # LOD 扫描数据(物理,含threshold列)
│   │   │   ├── cim_lod_data_mstmap.tsv      # LOD 扫描数据(图谱cM,含threshold列)
│   │   │   │                                #   已标注 chr_bp/pos_bp 物理坐标列
│   │   │   ├── cim_peaks_physical.tsv       # QTL 峰值表(物理位置)
│   │   │   ├── cim_peaks_mstmap.tsv         # QTL 峰值表(连锁图谱)
│   │   │   ├── cim_peaks_mstmap_physical.tsv # MSTmap 峰值标注物理坐标
│   │   │   ├── cim_threshold_physical.txt   # 显著性阈值
│   │   │   └── cim_threshold_mstmap.txt     # 显著性阈值
│   │   └── cim_analysis.R              # 生成的 R 脚本
│   └── post_rf/                         # Post-RF 分析结果(最终,结构同 pre_rf)
└── 99_logs/
    └── cim.log                          # 运行日志
```

physical / estimate 模式下 plots/ 内为 `cim_lod_data.tsv`、`cim_peaks.tsv`（无 mstmap 子目录）。

## 结果解读 | Interpreting Results

### 1. LOD 扫描数据（`cim_lod_data_*.tsv`）

**通俗理解|In plain words:** 这是最核心的结果——全基因组每个位置的「嫌疑程度打分表」。把 `lod` 列画成曲线，超过 `threshold` 水平线的山头，就是要找的 QTL。

```text
chr    pos        lod        threshold
LG1    0.0        0.0787     16.2882
LG1    1.2        0.2134     16.2882
...
```

- `chr`：连锁群或染色体编号；`pos`：位置（physical 为 Mb，mstmap 为 cM）；`lod`：LOD 得分；`threshold`：置换检验阈值
- **LOD ≥ threshold 即显著 QTL**；绘制全基因组 LOD 曲线时把 threshold 画成水平线，越过的峰就是 QTL
- mstmap 模式的 `cim_lod_data_mstmap.tsv` 额外带 `chr_bp`/`pos_bp` 列：每个 LOD 网格点标注的物理染色体与位置（两侧标记线性插值），方便把图谱峰直接对应回基因组坐标
- **不同图谱的 LOD/阈值不可直接比较**：物理模式坐标是 Mb(pseudo-cM)，与 mstmap 的 cM 尺度不同

### 2. QTL 峰值表（`cim_peaks_*.tsv`）

**通俗理解|In plain words:** 上面曲线里「越过水平线的山头」，整理成了一张清单——每个候选 QTL 在哪、峰值多高。

LOD 曲线中超过阈值的峰的位置表。mstmap 模式下 `cim_peaks_mstmap_physical.tsv` 把每个峰标注回物理坐标（chr_bp/pos_bp），**优先用它定位候选区间**。

### 3. 连锁图谱（`linkage_map.csv`）

**通俗理解|In plain words:** 用「重组率」重新排的座位表——告诉你在遗传距离(cM)的世界里，每个位点排在第几组(LG)第几号。

```text
marker           chr    pos
Chr01_732093     LG1    0.0
Chr01_6422855    LG12   23.697
```

- MSTmap 按重组频率把标记聚类为连锁群(LG)并计算 cM
- **LG 数通常多于实际染色体数**：重组冷点（如着丝粒区）会把一条染色体拆成多个 LG，属正常现象
- **碎片 LG 不进 CIM**：CIM 只使用每染色体标记数最多的那个 LG（碎片 LG 标记稀疏，R/qtl hk 方法在其上会产生数值不稳定的虚假 LOD 峰），被剔除的碎片 LG 会以 WARNING 记录在日志里
- 结合 `marker_map_index.tsv` 可将 LG/cM 映射回物理染色体与 bp

### 4. 物理图谱 vs 连锁图谱 | Physical vs linkage map

**通俗理解|In plain words:** 物理图谱按「门牌号」排队，连锁图谱按「开车时间」排队。两把尺子各有优缺点，两边的嫌疑区域都对得上，结论才最可信。

| | 物理图谱 | 连锁图谱 |
|---|---|---|
| 位置单位 | Mb（pseudo-cM） | cM |
| 标记顺序 | 按物理位置 | 按遗传距离 |
| 是否断链 | 不会 | 可能（重组冷点） |
| CIM 优势 | 完整、直接对应基因组 | 遗传距离更准确、定位精度更高 |

**两个图谱在同一区域都检出 QTL 信号，结果更可信。**

> **注意|Warning: mstmap 块 LOD 系统性偏高**：mstmap 的遗传图谱 cM 会膨胀（MSTmap 对 F2 用 RIL2 近似，距离标尺放大），导致其 LOD 与阈值整体抬升，**与 physical 块的 LOD 量级不可直接比较**。正确用法：**以 physical 块为最终结果**（峰位置、显著性），mstmap 块只用于「峰物理位置互证」——两块的显著峰落在同一物理区间即互为印证，不必纠结 mstmap 块 LOD 更大或排位不同。

### 5. QC 统计文件 | QC stats

**通俗理解|In plain words:** 这些文件回答「数据质量如何、过滤删掉了多少位点」。删得少(<10%)说明数据干净；删得多(>20%)说明要么数据有问题，要么阈值太严。

- `marker_filter_stats.txt` / `rf_filter_stats.txt`：各过滤步骤删除/保留的标记数与比例。删除比例异常高（如 RF 过滤删掉 >20%）提示基因型质量或阈值设置问题
- `local_hotspot_removed.tsv`：被删标记及原因。**错误标记常成簇**（连续多对过硬阈值共享标记），实际删除数低于「对数×2」属正常
- `singleton_het_report.tsv`：孤立重组标记的 H 比例诊断，辅助判断杂合异常

## 参数选择建议 | Parameter Guidance

- **`--map-mode`**：有参考基因组先跑默认 `mstmap`（一次拿两个图谱互相验证）；仅要物理定位用 `physical` 最快；无参考/高杂合群体可用 `estimate`
- **`--n-perm`**：探索性快速扫描可先 `--n-perm 0` 跳过置换看曲线形状；正式结果用 1000（耗时与标记数/样本数成正比）
- **`--local-hotspot`**：BSA 双亲群体保持默认开；高质量、低错配数据可 `--no-local-hotspot` 关闭
- **`--local-hotspot-relative`**：高重组背景的染色体（绝对阈值误杀）可给 1.0~2.0，阈值变为 max(绝对值, 系数×近距对 RF 中位数)
- **`--min-phys-gap`**：标记 >30 万时建议 300~500，QTL 分辨率无损且大幅提速；稀疏数据(<10 万标记)保持 0
- **`--skip-ld`**：默认跳过——输入若未预降维且标记高度冗余（>50 万），用 `--no-skip-ld` 启用
- **`--mstmap-pvalue`**：无需手调，内置自动调优（迭代 p 值直到核心 LG 数 ≤ 染色体数×3）
- **`--threads`**：仅模块直调可用，R/MSTmap 步骤的并行度有限，默认 1 即可

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件｜Input VCF file (.vcf/.vcf.gz) |
| `-p, --pheno` | 必填 |  | 表型文件｜Phenotype file (TSV: sample, value) |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-t, --type` | `f2` | f2/bc | 群体类型｜Cross type (f2/bc) |
| `--map-mode` | `mstmap` | physical/estimate/mstmap | cM位置来源｜cM source (physical/estimate/mstmap) |
| `--maf` | `0.05` | float | 最小等位基因频率｜Minor allele frequency threshold |
| `--missing` | `0.1` | float | 最大缺失率｜Maximum missing rate |
| `--max-het-rate` | `0.6` | float | H基因型最大比例｜Max heterozygous genotype rate |
| `--max-mean-rf` | `0.35` | float | K近邻局部平均RF最大值｜Max local mean RF to k nearest markers |
| `--rf-knn` | `10` | int | 局部RF窗口邻居数｜Local RF window neighbor count |
| `--local-hotspot/--no-local-hotspot` | `True` |  | 短距离重组热点过滤｜Local recombination hotspot filter |
| `--local-hotspot-dist` | `1000` | int | 热点判定物理距离阈值bp｜Hotspot physical distance threshold (bp) |
| `--local-hotspot-rf` | `0.2` | float | 相邻RF软阈值｜Soft adjacent-RF threshold |
| `--local-hotspot-hard-rf` | `0.3` | float | 相邻RF硬阈值(两侧都删)｜Hard adjacent-RF threshold |
| `--local-hotspot-score` | `1` | int | 软评分删除线｜Soft-score deletion cutoff |
| `--local-hotspot-relative` | `0.0` | float | 相对判据系数(0=关)｜Relative criterion factor (0=off) |
| `--min-phys-gap` | `0` | int | 最小相邻物理距离bp(0=关)｜Min physical gap bp (0=off) |
| `--n-marcovar` | `10` | int | 协因子数量｜Number of marker covariates |
| `--window` | `10.0` | float | 窗口大小(cM)｜Window size in cM |
| `--method` | `hk` | hk/em/imp | 扫描方法｜Scanning method |
| `--step` | `1.0` | float | 伪标记步长(cM)｜Pseudomarker step in cM |
| `--n-perm` | `1000` | int | 置换检验次数(0=跳过)｜Permutation count (0=skip) |
| `--ld-window` | `50` | int | LD窗口(SNP数)｜LD window in SNP count |
| `--ld-step` | `5` | int | LD步长(SNP数)｜LD step in SNP count |
| `--ld-r2` | `0.1` | float | LD r2阈值｜LD r2 threshold |
| `--skip-ld` | `True` |  | 跳过LD降维｜Skip LD pruning (default: True) |
| `--no-skip-ld` | — |  | 启用LD降维｜Enable LD pruning |
| `--mstmap-pvalue` | `1e-06` | float | MSTmap聚类p值起始值(自动调优)｜MSTmap clustering p-value start, auto-tuned |
| `--mstmap-distfun` | `kosambi` | kosambi/haldane | MSTmap距离函数｜MSTmap distance function |
| `--r-env` | `~/miniforge3/envs/Rqtl` |  | R conda环境路径或名称｜R conda env path or name |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入VCF文件｜Input VCF file (supports .vcf and .vcf.gz) |
| `-p, --pheno` | 必填 |  | 表型文件｜Phenotype file (TSV with sample, value columns) |
| `-o, --output` | 必填 |  | 输出目录｜Output directory |
| `-t, --type` | `f2` | f2/bc | 群体类型｜Cross type (f2/bc, default: f2) |
| `--map-mode` | `mstmap` | physical/estimate/mstmap | cM位置来源｜cM position source (physical/estimate/mstmap, default: mstmap) |
| `--maf` | `0.05` | float | 最小等位基因频率阈值｜Minor allele frequency threshold (default: 0.05) |
| `--missing` | `0.1` | float | 最大缺失率｜Maximum missing rate (default: 0.1) |
| `--max-het-rate` | `0.6` | float | H基因型最大比例(规则1)｜Max heterozygous genotype rate (rule 1, default: 0.6) |
| `--max-mean-rf` | `0.35` | float | K近邻局部平均RF最大值(规则2)｜Max local mean RF to k nearest markers (rule 2, default: 0.35) |
| `--rf-knn` | `10` | int | 规则2局部RF窗口邻居数(前后各一半)｜Local RF window neighbor count (default: 10) |
| `--no-local-hotspot` | — | store_false | 关闭短距离重组热点过滤｜Disable local hotspot filter |
| `--local-hotspot-dist` | `1000` | int | 热点判定的物理距离阈值bp(默认1000)｜Hotspot physical distance threshold in bp (default: 1000) |
| `--local-hotspot-rf` | `0.2` | float | 相邻RF软阈值(默认0.20)｜Soft adjacent-RF threshold (default: 0.20) |
| `--local-hotspot-hard-rf` | `0.3` | float | 相邻RF硬阈值: 达到即两侧标记都删(默认0.30)｜Hard adjacent-RF threshold, delete both (default: 0.30) |
| `--local-hotspot-score` | `1` | int | 软评分删除线(默认1)｜Soft-score deletion cutoff (default: 1) |
| `--local-hotspot-relative` | `0.0` | float | 相对判据系数(0=只用绝对阈值)｜Relative criterion factor, 0=absolute only (default: 0) |
| `--min-phys-gap` | `0` | int | 最小相邻物理距离bp, 0=关闭(默认)｜Min physical gap in bp, 0=off (default: 0) |
| `--n-marcovar` | `10` | int | 协因子数量｜Number of marker covariates (default: 10) |
| `--window` | `10.0` | float | 窗口大小(cM)｜Window size in cM (default: 10) |
| `--method` | `hk` | hk/em/imp | 扫描方法｜Scanning method (default: hk) |
| `--step` | `1.0` | float | 伪标记步长(cM)｜Pseudomarker step in cM (default: 1) |
| `--n-perm` | `1000` | int | 置换检验次数(0=跳过)｜Number of permutations (0=skip, default: 1000) |
| `--ld-window` | `50` | int | LD计算窗口(SNP数)｜LD window in SNP count (default: 50) |
| `--ld-step` | `5` | int | LD计算步长(SNP数)｜LD step in SNP count (default: 5) |
| `--ld-r2` | `0.1` | float | LD r2阈值｜LD r2 threshold (default: 0.1) |
| `--skip-ld` | — | store_true | 跳过LD降维步骤｜Skip LD pruning |
| `--mstmap-pvalue` | `1e-06` | float | MSTmap聚类p值起始值(自动调优)｜MSTmap clustering p-value start value, auto-tuned (default: 1e-6) |
| `--mstmap-distfun` | `kosambi` | kosambi/haldane | MSTmap距离函数｜MSTmap distance function (default: kosambi) |
| `--mstmap-path` | `~/miniforge3/envs/r/bin/mstmap` |  | MSTmap二进制路径｜MSTmap binary path (default: ~/miniforge3/envs/r/bin/mstmap) |
| `--r-env` | `~/miniforge3/envs/Rqtl` |  | R conda环境路径或名称｜R conda env path or name (default: ~/miniforge3/envs/Rqtl) |
| `--threads` | `1` | int | 并行线程数｜Number of parallel threads (default: 1) |
| `--version` | — | version |  |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- Python 3
- R 环境（需安装 qtl 包；通过 `--r-env` 指定）
- MSTmap（仅 mstmap 模式；`--mstmap-path` 指定，默认 `~/miniforge3/envs/r/bin/mstmap`）
- bcftools（提取过滤后的 VCF）
- PLINK（LD 降维计算，仅启用 LD 时）

## 常见问题 | FAQ

**Q1：换参数重跑，为什么结果没变？**
断点续传按输出文件存在性判断。换过滤参数（如 `--maf`、`--min-phys-gap`）重跑旧输出目录前，先删除对应的旧产物（如 `01_qc/filtered_markers.vcf.gz`），否则会复用旧参数的结果。

**Q2：物理模式和 mstmap 模式的 LOD 能直接对比吗？**
不能。物理模式坐标是 Mb（pseudo-cM），window/step 与 mstmap 的 cM 尺度不同，LOD/阈值只在同一图谱内比较。

**Q3：为什么 LG 数比染色体数多？**
重组冷点（着丝粒等）使一条染色体被拆成多个 LG，正常现象。用 `marker_map_index.tsv` 把 LG 映射回物理染色体。

**Q4：RF 过滤删了很多标记正常吗？**
删除比例 <10% 通常正常；错误标记成簇删除时实际删除数低于「对数×2」。若 RF 过滤删 >20%，检查基因型错误率与 `--max-mean-rf`/`--local-hotspot` 阈值是否过严。

**Q5：`--threads`/`--mstmap-path` 在 `biopytools cim` 里没有？**
这两个参数只在直接调用 `python -m biopytools.cim` 时暴露，click 包装器用默认值（threads=1，mstmap 走 `--mstmap-path` 默认路径）。

**Q6：mstmap 的 LOD 文件里 chr_bp/pos_bp 列是什么？**
为 LOD 网格点标注的物理坐标（两侧标记线性插值），用于把图谱峰对应回基因组；单标记 LG 直接取该标记坐标。重复运行脚本不会重复追加这两列（幂等守卫）。

**Q7：mstmap 块的显著峰和 physical 块不在同一条染色体上？**
以 physical 块为准。mstmap 块 LOD 系统性偏高（见「结果解读 §4」），且碎片 LG（标记稀疏的小连锁群）不进 CIM；若两个块在**真 QTL 区间**（physical 块最高峰）附近都有信号，即互相印证。只有 mstmap 块某碎片区域出现孤立高 LOD 峰时，多为 hk 方法的数值伪峰（标记本身与表型无关联），可忽略。
