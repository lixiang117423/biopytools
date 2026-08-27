# VCF 转 SplitsTree6 距离矩阵 | vcf2splitstree

把群体的变异数据(VCF)转成 SplitsTree6 能直接打开的距离矩阵文件——转换在超算上跑完,拿到的 CSV 用你自己电脑上的 SplitsTree6 图形界面打开,自动生成进化网络(NeighborNet),可以交互式调整展示。

## 功能概述 | Overview { #overview }

- **只做转换**:VCF 变异文件 → p-distance 距离矩阵 CSV,不依赖 SplitsTree6 本体(超算上不用装 GUI 软件)
- **numpy 向量化**:基因型展开为联合等位编码,矩阵乘法算距离,365 样本×几十万位点几分钟完成(旧实现约 1 小时)
- **缺失处理**:`./.` 缺失基因型在对应样本对比较中跳过,不中断
- **正确语义**:`0/1` 与 `1/0` 算位点差异(等价 GT 字符串比较),相位符 `|` 与 `/` 一视同仁
- **SplitsTree6 一键打开**:CSV 格式按官方 CSVReader 要求生成(首行含逗号),GUI 打开自动识别为距离矩阵并跑 NeighborNet(已实测 fit 99.9%)
- 断点续传:输出已存在则跳过重算

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools vcf2splitstree -i variants.vcf.gz -o vcf2s_out/
```

转换完成后,把 `vcf2s_out/variants.distances.csv` 拷回 Mac,用 SplitsTree6 的 File→Open 打开即可。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗解释<br>In plain words |
|------|------|
| p-distance | 两两样本间"不同变异位点比例",0=完全相同,1=全部不同;最朴素的遗传距离 |
| 距离矩阵 | 每对样本一个距离值的方阵,是 SplitsTree6/NeighborNet 的输入 |
| NeighborNet | SplitsTree6 招牌算法:把距离矩阵画成网状进化图,能展示重组/杂交信号 |
| 缺失基因型 | VCF 里 `./.` 表示该位点没测到,比较时跳过 |

## 输入 | Input { #input }

- **`.vcf` 或 `.vcf.gz`** 变异文件(双等位 SNP 效果最佳)
- 必需有 `#CHROM` 表头行,样本名列在表头第 10 列起
- 支持 GT 字段的 `/` 与 `|` 分隔;`./.`、`.|.`、`.` 视为缺失
- 非 VCF 文件会在参数校验时报错

## 参数说明 | Parameters { #parameters }

**通俗理解|In plain words:** 只有三个参数——`-i` 输入 VCF、`-o` 输出目录、`--log-level` 日志级别,没有别的可调项,转换逻辑是固定的。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | VCF 变异文件(.vcf/.vcf.gz)｜VCF variants file |
| `-o, --output-dir` | `./vcf2splitstree_output` |  | 输出目录(默认./vcf2splitstree_output)｜Output directory |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | VCF 变异文件(.vcf/.vcf.gz)｜VCF variants file |
| `-o, --output-dir` | `./vcf2splitstree_output` |  | 输出目录(默认./vcf2splitstree_output)｜Output directory |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

<!-- END PARAMS:auto -->

## 输出 | Output { #output }

```
output/
├── variants.distances.csv     # 距离矩阵 CSV(SplitsTree6 直接打开)
├── 00_pipeline_info/
│   └── software_versions.yml  # 输入/输出/样本数/耗时
└── 99_logs/
    └── vcf2splitstree.log     # 运行日志
```

**CSV 格式**:每行一个样本,`样本名,距离1,距离2,...,距离n`,无首行计数——这是 SplitsTree6 CSVReader 自动识别的格式。

## 结果解读 | Interpreting Results { #interpreting-results }

- **在 Mac 上打开**:SplitsTree6 → File → Open → 选 `variants.distances.csv` → 自动识别距离矩阵 → 自动跑 NeighborNet
- **fit 值**:网络对距离矩阵的拟合度(GUI 里 Splits Network 窗口可看),**>90 说明网络能很好表达数据**
- **矩阵质量检查**:对角线全 0;两样本距离接近 0 时警惕测序混样;CSV 行数 = 样本数
- **平行边(网状结构)**:NeighborNet 的网格状分支表示存在不一致的系统发育信号(重组/杂交/基因流),这是网络的特色不是错误

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

| 场景<br>Scenario | 建议<br>Recommendation |
|------|------|
| 常规 SNP 数据 | 默认即可,直接转换 |
| VCF 巨大(>50 万位点) | 建议先用筛选工具减位点再转(可显著加快) |
| 多等位位点多 | 转换只处理双等位,多等位自动跳过(不报错) |

**通俗理解|In plain words:** 不用调参,转完拿 CSV 去 Mac 上开。

## 依赖 | Dependencies { #dependencies }

| 依赖<br>Dependency | 说明<br>Note |
|------|------|
| Python 3.10 + numpy | 唯一依赖;转换在超算完成,不需要 SplitsTree6 本体 |
| SplitsTree6(可选) | 仅在你 Mac 上查看结果时需要(GUI 打开 CSV) |

## 常见问题 | FAQ { #faq }

- **SplitsTree6 打开 CSV 报"Unknown format"?** 确认文件是模块生成的(首行是 `样本名,数值...` 含逗号);手工改过的文件可能破坏格式
- **转换后样本对距离全是 0?** 检查 VCF 是否真的有多态位点(全部 `0/0` 会得到全 0 矩阵);或样本标签是否重复
- **为什么不在超算上直接出图?** SplitsTree6 是 GUI 软件(JavaFX),无头运行需要 Xvfb 虚拟屏且工作流格式复杂;转换后在你的 Mac 上打开最省事、交互最好
- **多等位位点会怎样?** 自动跳过,不影响其他位点计算
