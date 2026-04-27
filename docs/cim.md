# CIM - 复合区间作图分析 (Composite Interval Mapping)

基于 R/qtl 的复合区间作图(CIM)分析工具，用于BSA群体的QTL定位。

## 用法

```bash
python -m biopytools.cim -i input.vcf.gz -p phe.txt -o output_dir
```

### 必需参数

| 参数 | 说明 |
|------|------|
| `-i, --input` | 输入VCF文件（支持 .vcf / .vcf.gz） |
| `-p, --pheno` | 表型文件（TSV格式，含 sample, value 两列） |
| `-o, --output` | 输出目录 |

### 群体类型

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `-t, --type` | f2 | 群体类型：f2（F2群体）或 bc（回交群体） |

### 遗传图谱模式

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--map-mode` | mstmap | cM位置来源，可选 physical / estimate / mstmap |

- **physical**：直接用基因组物理位置(bp)，不做图谱构建
- **estimate**：由R/qtl的 `est.map()` 估算遗传距离
- **mstmap**：由MSTmap软件构建连锁图谱，同时跑physical和mstmap两次CIM

### 标记过滤

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--maf` | 0.05 | 最小等位基因频率(MAF)阈值 |
| `--missing` | 0.1 | 标记最大缺失率 |

### 重组频率质控

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--max-het-rate` | 0.6 | 杂合基因型(H)最大比例，超过则删除该标记 |
| `--max-mean-rf` | 0.5 | 同染色体平均重组频率(RF)最大值，超过则删除该标记 |

### LD降维

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--ld-window` | 50 | LD计算窗口（SNP数） |
| `--ld-step` | 5 | LD计算步长（SNP数） |
| `--ld-r2` | 0.1 | LD r²阈值 |
| `--skip-ld` | false | 跳过LD降维 |

### CIM参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--n-marcovar` | 10 | 协因子数量 |
| `--window` | 10.0 | 窗口大小(cM) |
| `--method` | hk | 扫描方法：hk / em / imp |
| `--step` | 1.0 | 伪标记步长(cM) |

### 置换检验

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--n-perm` | 1000 | 置换检验次数（0=跳过） |

### MSTmap参数（仅mstmap模式）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--mstmap-pvalue` | 1e-6 | 聚类p值起始值（自动调优） |
| `--mstmap-distfun` | kosambi | 距离函数：kosambi / haldane |
| `--mstmap-path` | ~/miniforge3/envs/Rqtl/bin/mstmap | MSTmap二进制路径 |

### 环境参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--r-env` | Rqtl | R conda环境名 |
| `--threads` | 1 | 并行线程数 |

## 分析流程

整个流程分 Pre-RF 和 Post-RF 两轮，先跑一轮不做RF过滤的基线分析，再过滤后跑最终分析，方便对比。

```
输入VCF + 表型文件
    │
    ▼
步骤1: 解析输入
    VCF基因型 → ABH编码矩阵 (A=纯合父本, H=杂合, B=纯合母本)
    │
    ▼
步骤2: 标记过滤
    ├─ MAF过滤: 次等位基因频率太低 → 删除
    └─ 缺失率过滤: 缺失数据太多 → 删除
    │
    ▼
步骤3: LD降维（可选）
    滑动窗口内LD r²过高的标记 → 去冗余保留子集
    │
    ├──────────────────────────────┐
    ▼                              │
Pre-RF CIM分析（基线，不做RF过滤）  │
    ├─ 构建csvs格式文件             │
    ├─ [mstmap模式] MSTmap构建连锁图谱  │
    ├─ [mstmap模式] 用物理位置跑CIM    │
    ├─ [mstmap模式] 用连锁图谱跑CIM    │
    ├─ [estimate模式] 估算遗传距离后CIM │
    └─ [physical模式] 直接用物理位置CIM  │
                                   │
    ◄──────────────────────────────┘
    │
    ▼
步骤5: 重组频率(RF)质控
    ├─ H比例过滤: 杂合基因型比例 > max_het_rate → 删除
    ├─ 平均RF过滤: 同染色体平均RF > max_mean_rf → 删除
    └─ 孤立重组检测: 与相邻标记RF均 > 0.5 → 仅警告
    │
    ▼
Post-RF CIM分析（最终结果）
    （同Pre-RF流程，使用RF过滤后的标记）
    │
    ▼
输出结果
```

### 重组频率(RF)计算

对同一染色体内的标记两两计算重组频率：

| 两个样本的基因型 | 判定 | 信号值 |
|---|---|---|
| A-A 或 B-B | 不重组 | 0 |
| A-B 或 B-A | 重组 | 1 |
| 含 H | 不确定 | 0.5 |
| 含 NaN | 跳过 | - |

```
RF(标记i, 标记j) = 所有有效样本的信号之和 / 有效样本数
平均RF(标记i) = 与同染色体所有其他标记的RF之和 / (n-1)
```

### 物理图谱 vs 连锁图谱

| | 物理图谱 | 连锁图谱 |
|---|---|---|
| 位置单位 | bp（碱基对） | cM（厘摩） |
| 位置来源 | 参考基因组 | MSTmap根据重组频率计算 |
| 标记顺序 | 按物理位置 | 按遗传距离 |
| 是否断链 | 不会 | 可能（重组冷点导致） |
| CIM优势 | 完整性好，直接对应基因组位置 | 遗传距离更准确，QTL定位精度更高 |

mstmap模式下两种图谱各跑一次CIM，互相补充验证。两个图谱在同一区域都检测到QTL信号则结果更可信。

## 输入文件

### VCF文件

标准的VCF格式，基因型编码为0/0、0/1、1/1等，代码内部自动转为ABH编码。

### 表型文件

TSV格式，两列：

```
sample    trait1
21-18     0
21-19     1
...
```

## 输出目录结构

```
output_dir/
├── 00_pipeline_info/
│   └── pipeline_params.txt          # 流程参数记录
├── 01_qc/
│   ├── marker_filter_stats.txt      # MAF/缺失率过滤统计
│   ├── ld_prune_stats.txt           # LD降维统计
│   ├── rf_filter_stats.txt          # RF质控过滤统计
│   ├── singleton_het_report.tsv     # 孤立重组标记H比例诊断
│   ├── filtered_markers.vcf.gz      # LD降维后保留的标记VCF
│   └── rf_filtered_markers.vcf.gz   # RF过滤后保留的标记VCF
├── 02_cim/
│   ├── pre_rf/                      # Pre-RF分析结果（基线）
│   │   ├── tidy_files/              # csvs格式输入文件
│   │   ├── physical/                # 物理位置CIM结果
│   │   │   ├── cim_genome_plot.pdf
│   │   │   ├── cim_chr_*.pdf
│   │   │   ├── cim_results.rds
│   │   │   └── cim_perm_results.rds
│   │   ├── mstmap/                  # 连锁图谱CIM结果
│   │   │   ├── linkage_map.csv      # MSTmap连锁图谱（标记→LG→cM）
│   │   │   ├── mstmap_map.csv       # MSTmap标记位置
│   │   │   ├── mstmap_gen.csv       # MSTmap基因型矩阵
│   │   │   ├── mstmap_phe.csv       # MSTmap表型
│   │   │   ├── marker_map_index.tsv # 标记物理位置↔连锁图谱位置索引
│   │   │   ├── cim_genome_plot.pdf
│   │   │   ├── cim_chr_*.pdf
│   │   │   ├── cim_results.rds
│   │   │   └── cim_perm_results.rds
│   │   ├── plots/                   # 可视化数据（供外部绘图用）
│   │   │   ├── cim_lod_data_physical.tsv   # LOD扫描数据（含threshold列）
│   │   │   ├── cim_lod_data_mstmap.tsv     # LOD扫描数据（含threshold列）
│   │   │   ├── cim_peaks_physical.tsv      # QTL峰值表（物理位置）
│   │   │   ├── cim_peaks_mstmap.tsv        # QTL峰值表（连锁图谱）
│   │   │   ├── cim_peaks_mstmap_physical.tsv # MSTmap峰值标注物理坐标
│   │   │   ├── cim_threshold_physical.txt  # 显著性阈值
│   │   │   └── cim_threshold_mstmap.txt    # 显著性阈值
│   │   └── cim_analysis.R          # 生成的R脚本
│   └── post_rf/                     # Post-RF分析结果（最终，结构同pre_rf）
└── 99_logs/
    └── cim.log                     # 运行日志
```

## 关键输出文件说明

### LOD扫描数据 (cim_lod_data_*.tsv)

```
chr    pos        lod        threshold
LG1    0.0        0.0787     16.2882
LG1    1.2        0.2134     16.2882
...
```

- `chr`：连锁群或染色体编号
- `pos`：位置（物理模式为bp，MSTmap模式为cM）
- `lod`：LOD得分
- `threshold`：置换检验得到的显著性阈值（LOD超过此值为显著QTL）

### QTL峰值表 (cim_peaks_*.tsv)

LOD曲线中超过阈值的峰值位置。

### 连锁图谱 (linkage_map.csv)

```
marker           chr    pos
Chr01_732093     LG1    0.0
Chr01_6422855    LG12   23.697
```

MSTmap根据重组频率将标记聚类为连锁群(LG)，并计算遗传距离(cM)。连锁群数量通常大于实际染色体数，因为重组冷点（如着丝粒区域）会导致同一条染色体被拆分成多个LG。

## 依赖

- Python 3
- R环境（需安装 qtl 包）
- MSTmap（仅 mstmap 模式，默认路径 `~/miniforge3/envs/Rqtl/bin/mstmap`）
- bcftools（用于提取过滤后的VCF）
- PLINK（用于LD降维计算）
