# BSA分析套件 | OcBSA BSA Analysis Tool Suite

**基于DHHP/SNP-index/ED算法的F1/F2群体BSA分析, 含可视化和引物设计 | F1/F2 BSA analysis suite with visualization and primer design**

## 功能概述 | Overview

ocbsa 模块是完整的BSA(Bulked Segregant Analysis)混池分析套件, 提供四个子命令: f1 (F1群体DHHP算法分析)、f2 (F2/RILs群体SNP-index/ED分析)、fig (结果可视化绘图)、primer (候选区域引物设计)。DHHP算法适用于显性杂交F1群体, SNP-index和ED适用于F2或RILs分离群体。模块配套绘图和引物设计工具, 形成从定位到验证的完整流程。

## 快速开始 | Quick Start

```bash
# F1群体DHHP分析
biopytools ocbsa f1 -i input.vcf -p1 1 -p2 2 -b1 3 -b2 4 -o ./output

# F2群体SNP-index分析
biopytools ocbsa f2 -i input.vcf -p1 1 -p2 2 -b1 3 -b2 4 --method snpindex -o ./output

# 结果可视化
biopytools ocbsa fig -i result.txt -o output.png --plot-type ocvalue

# 候选区域引物设计
biopytools ocbsa primer -g genome.fa -i result.OcValue --region Chr01,100000,500000 -o ./output
```

## 子命令 | Subcommands

### f1 - F1群体DHHP分析 | F1 DHHP Analysis

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input-vcf` | 必需 | VCF文件路径 |
| `-p1`, `-p2` | 必需 | 显性/隐性亲本列号 |
| `-b1`, `-b2` | 必需 | 显性/隐性表型混池列号 |
| `-o, --output-dir` | `./output` | 输出目录 |
| `-w, --window-size` | `1000000` | 滑窗大小 |
| `-p, --pvalue` | `99` | p值阈值 |
| `--parent-min-dep/--max-dep` | `10/100` | 亲本最低/最高覆盖度 |
| `--pool-min-dep/--max-dep` | `20/500` | 混池最低/最高覆盖度 |

### f2 - F2/RILs群体分析 | F2/RILs Analysis

参数与 f1 基本一致, 额外支持:

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--method` | `snpindex` | 分析方法(snpindex/ED) |

### fig - 可视化 | Visualization

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input-file` | 必需 | 输入结果文件 |
| `-o, --output-file` | 必需 | 输出图片路径(.png/.pdf) |
| `--plot-type` | `ocvalue` | 图表类型(ocvalue/snpindex/ed) |
| `--position` | `None` | 特定区域(chr,start,end) |
| `--color` | `plasma_r` | 颜色方案 |

### primer - 引物设计 | Primer Design

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-g, --genome` | 必需 | 参考基因组路径 |
| `-i, --ocvalue-file` | 必需 | OcValue文件路径 |
| `--region` | 必需 | 目标区间(chr,start,end) |
| `-o, --output-dir` | `./output` | 输出目录 |
| `-n, --primer-num` | `10` | 引物对数量 |
| `--flank-length` | `200` | INDEL侧翼长度 |
| `--primer-min/opt/max-size` | `18/20/24` | 引物长度范围 |
| `--product-min/max` | `70/200` | 产物长度范围 |
| `--min-tm/--max-tm` | `50.0/65.0` | Tm范围 |
| `--min-gc/--max-gc` | `35.0/65.0` | GC含量范围 |

(运行 `biopytools ocbsa <subcommand> -h` 查看完整参数列表)

## 输出 | Output

- f1/f2: `{output}/result.OcValue` 或 `result.snpindex` 滑窗结果文件
- fig: PNG/PDF图片
- primer: `{output}/primers.txt` 引物列表 + `{output}/inDel.fa` INDEL序列

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

- **primer3**: 引物设计 (用于primer子命令)
- **Python库**: pandas, numpy, matplotlib

## 引用 | Citation

- Zou C., et al. OcBSA: an integrated bulked segregant analysis workflow for extreme phenotype. (基于DHHP算法)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
