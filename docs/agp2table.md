# AGP 转表格工具 | AGP to Table Converter

一句话理解：**把基因组组装里那张晦涩的 AGP 文件，翻译成一张人看得懂的表格**——每行说清「某个 scaffold 上哪一段、是真实序列还是缺口」，还能顺手统计每个 scaffold 的长度和缺口数量。

## 功能概述 | Overview { #overview }

- 把 AGP（9 列、Tab 分隔）转成带表头的易读表格，支持 `txt` / `tsv` / `csv` / `xlsx` 四种输出
- 自动识别组件类型：`W`（真实序列 contig）、`U` / `N`（缺口 gap）、`O`（其他）
- 可选添加统计信息：总长度、contig/gap 数量、gap 总长、前 20 长 scaffold 排行
- 默认按 scaffold 名排序后分组输出，可关掉恢复原始顺序；格式错误行自动跳过，不中断
- 纯 Python 实现，无外部生信软件依赖

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools agp2table -i assembly.agp -o assembly_table.txt
```

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗理解 |
|---|---|
| AGP 文件 | 组装结果的「施工清单」，逐行说明一条 scaffold 由哪些 contig 和缺口按什么顺序拼起来 |
| scaffold | 组装里较长的「大段序列」，由多个 contig 加缺口拼接而成，像一段拼好的墙 |
| contig | 没有缺口的「连续纯序列」，是拼墙时最小的完整砖块 |
| gap（缺口） | 拼接时不确定具体碱基、用一串 N 表示的空档，像墙砖之间的水泥缝 |
| 组件类型（W/U/N/O） | 每行说明这段是「砖」（W）还是「缝」（U/N），O 是其他特殊情况 |
| 方向（orientation） | 这段序列正着放（+）还是反过来放（-），像砖块的正反面 |

## 输入 | Input { #input }

标准 AGP 文件，Tab 分隔，每行至少 9 列（不足 9 列的行会被跳过并记警告），`#` 开头和空行会被忽略。

```text
group1	1	190127	1	W	ptg002519l	1	190127	+
group1	190128	190227	2	U	100	scaffold	yes	proximity_ligation
group1	190228	5516498	3	W	ptg000008l	1	5326271	-
```

9 列含义依次是：scaffold 名、起始位置、结束位置、部件编号、组件类型、组件 ID、组件起始、组件结束、方向。第 5 列的 `W` / `U` / `N` / `O` 决定这一行被当作 contig 还是 gap。

## 参数说明 | Parameters { #parameters }

### 输出格式 | Output format

**通俗理解|In plain words:** 决定结果文件长什么样。`txt` 和 `tsv` 其实都是 Tab 分隔、内容一模一样，选哪个都行；`csv` 用逗号分隔，方便 Excel 直接打开；`xlsx` 是真正的 Excel 工作簿，统计信息会放到独立的 `Statistics` 工作表里，最漂亮但需要装 `openpyxl`。**一般用默认 `txt` 即可；要交给不写代码的同事看就选 `xlsx`。**

相关参数：`-f, --format`（默认 `txt`）。

### 统计与表头 | Statistics and headers

**通俗理解|In plain words:** `--statistics` 打开后，输出文件末尾多出一段「总长度、contig/gap 数量、前 20 长 scaffold」的汇总，方便一眼看懂组装规模，默认关闭；`--no-headers` 关掉表头，适合下游脚本直接读纯数据；`--no-grouping` 关掉按 scaffold 分组，保留输入文件的原始行顺序（默认会按 scaffold 名排序后分组输出）。**三者都是可选开关，需要什么开什么，一般不用动。**

相关参数：`--statistics`、`--no-headers`、`--no-grouping`。

## 输出 | Output { #output }

本工具只输出一个表格文件（路径由 `-o` 指定），并在同一目录生成日志 `agp2table.log`。

```text
输出目录/
├── assembly_table.txt    # 转换后的表格（表头 + 逐行记录，可选统计段）
└── agp2table.log         # 运行日志
```

表格列（默认带表头）：

| 列名<br>Column | 含义 |
|---|---|
| Scaffold | scaffold / 染色体名 |
| Start / End | 该段在 scaffold 上的起止坐标 |
| Length | 该段长度（bp） |
| Gap_Count | 该 scaffold 内的缺口总数 |
| Part_Number | 该段在 scaffold 里的部件编号 |
| Type | 组件类型（Contig / Gap / Other） |
| Component_ID | contig 名称 |
| Component_Start / Component_End | 该段在 contig 上的起止坐标 |
| Orientation | 方向（+ / - / ?） |

## 结果解读 | Interpreting Results { #interpreting-results }

- 看 `Type` 列：`Contig` 行是真实序列，`Gap` 行是拼接缺口，两者交替出现是正常组装结构
- 看 `Gap_Count`：同一 scaffold 的所有行该值相同，等于该 scaffold 内的缺口数
- 开 `--statistics` 后，末尾「Gap 总长」占「总长度」的比例越低，说明拼接越完整
- 「前 20 个最长的 scaffold」排行里，第一名约等于组装最长的一条染色体 / scaffold

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

- 日常查看：默认 `-f txt`，不开统计
- 给同事看 / 归档：`-f xlsx --statistics`（统计进独立工作表）
- 喂给下游脚本：`--no-headers --no-grouping`，拿到干净有序的原始数据
- 输入行顺序很重要时：加 `--no-grouping`，保留原始顺序

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `--input, -i` | 必填 |  | AGP文件路径｜AGP file path |
| `--output, -o` | 必填 | Path | 输出表格文件路径｜Output table file path |
| `--format, -f` | `txt` | txt/tsv/csv/xlsx | 输出格式｜Output format |
| `--statistics` | — |  | 添加统计信息｜Add statistics |
| `--no-headers` | — |  | 不添加表头｜Do not add headers |
| `--no-grouping` | — |  | 不按scaffold分组｜Do not group by scaffold |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | AGP文件路径｜AGP file path |
| `-o, --output` | 必填 |  | 输出表格文件路径｜Output table file path |
| `-f, --format` | `txt` | txt/tsv/csv/xlsx | 输出格式｜Output format |
| `--statistics` | — | store_true | 添加统计信息｜Add statistics |
| `--no-headers` | — | store_true | 不添加表头｜Do not add headers |
| `--no-grouping` | — | store_true | 不按scaffold分组｜Do not group by scaffold |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies { #dependencies }

- Python 3
- `openpyxl`（仅输出 `xlsx` 时需要；缺了会报错并提示 `pip install openpyxl`）
- 无 conda 环境、无外部生信软件依赖

## 常见问题 | FAQ { #faq }

**Q1：输出 `xlsx` 报「openpyxl 未安装」？**
装一下即可：`pip install openpyxl`，或改用 `-f tsv` 输出文本表格。

**Q2：为什么有的行没出现在结果里？**
输入行少于 9 列会被跳过并在日志里记警告，`#` 开头和空行也会被忽略。用 `head` 检查输入文件的格式。

**Q3：`txt` 和 `tsv` 有什么区别？**
没有区别，两者都是 Tab 分隔、内容一致；仅 `csv`（逗号）和 `xlsx`（Excel）是不同格式。

**Q4：会断点续传吗？**
不会。这是一次性转换，重跑会直接覆盖输出文件，没有步骤跳过。
