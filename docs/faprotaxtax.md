# FAPROTAX功能注释 | FAPROTAX Functional Annotation

**用FAPROTAX对OTU/ASV表进行微生物群落功能(生态)注释 | Annotate OTU/ASV tables with FAPROTAX for microbial ecological functions**

## 功能概述 | Overview

faprotaxtax 模块封装了 [FAPROTAX](http://www.loucalab.com/archive/FAPROTAX/), 基于 OTU/ASV 丰度表(结合分类信息)推断群落生态功能(如碳/氮/硫循环、病原性等)。底层调用 FAPROTAX 的 `collapse_table.py` 将分类单元丰度聚合到功能组。

## 快速开始 | Quick Start

```bash
# 基础用法(BIOM 输入)
biopytools faprotaxtax -i otu_table.biom -o faprotaxtax_output

# 按元数据字段(如 taxonomy)折叠
biopytools faprotaxtax -i otu_table.tsv -o out/ --collapse-by-metadata taxonomy
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input-table` | 输入 OTU/ASV 表(BIOM 或 TSV) |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./faprotaxtax_output` | 输出目录 |
| `-g, --groups-file` | 内置 | FAPROTAX 功能组数据库文件 |
| `--collapse-table-path` | 内置 | collapse_table.py 脚本路径 |
| `--collapse-by-metadata` | — | 用于功能注释的 BIOM 元数据字段(如 taxonomy) |
| `--group-leftovers-as` | — | 未匹配功能组的记录归为此组名 |
| `-n, --normalize` | `none` | 标准化方式 |
| `--average` | `none` | 组内聚合方式 |
| `--output-format` | `auto` | 输出格式(auto/BIOM/classical) |
| `-t, --threads` | `1` | 线程数 |
| `-f, --force` | `False` | 覆盖已存在的输出 |

(运行 `biopytools faprotaxtax -h` 查看完整参数列表)

## 输出 | Output

- 功能组丰度表(collapse 后)及相关中间文件, 写入输出目录
- 运行日志

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-table` | 必填 |  | 输入OTU/ASV表（BIOM或TSV格式）｜Input OTU/ASV table (BIOM or TSV format) |
| `-o, --output-dir` | `./faprotaxtax_output` |  | 输出目录｜Output directory |
| `-g, --groups-file` | — |  | FAPROTAX功能组数据库文件路径｜Path to FAPROTAX groups database file |
| `--collapse-table-path` | — |  | collapse_table.py脚本路径｜Path to collapse_table.py script |
| `--collapse-by-metadata` | — |  | 用于功能注释的BIOM元数据字段名｜BIOM metadata field for functional annotation |
| `--group-leftovers-as` | — |  | 未匹配到功能组的记录归为此组名｜Group name for unmatched records |
| `-n, --normalize` | `none` | none/columns_before_collapsing/rows_before_collapsing/columns_after_collapsing/rows_after_collapsing/columns_before_collapsing_excluding_unassigned/rows_before_collapsing_excluding_unassigned | 标准化方式｜Normalization method |
| `--average` | `none` | none/across_records/across_group_members/across_used_group_members/maximum/minimum/minimum_across_records | 组内聚合方式｜Aggregation method within groups |
| `--row-names-are-in-column` | — |  | 包含行名的列名或索引｜Column name or index containing row names |
| `--output-format` | `auto` | auto/BIOM/classical | 输出格式｜Output format |
| `-t, --threads` | `1` | int | 线程数｜Number of threads |
| `-f, --force` | — |  | 覆盖已存在的输出文件｜Overwrite existing output files |
| `-v, --verbose` | — |  | 详细输出｜Verbose output |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input-table` | 必填 |  | 输入OTU/ASV表（BIOM或TSV格式）｜Input OTU/ASV table (BIOM or TSV format) |
| `-o, --output-dir` | `./faprotaxtax_output` |  | 输出目录｜Output directory (default: ./faprotaxtax_output) |
| `-g, --groups-file` | — |  | FAPROTAX功能组数据库文件路径｜Path to FAPROTAX functional groups database file |
| `--collapse-table-path` | — |  | collapse_table.py脚本路径｜Path to collapse_table.py script |
| `--collapse-by-metadata` | — |  | 用于功能注释的BIOM元数据字段名（如: taxonomy）｜BIOM metadata field for functional annotation (e.g., taxonomy) |
| `--group-leftovers-as` | — |  | 未匹配到功能组的记录归为此组名｜Group name for records not matching any functional group |
| `-n, --normalize` | `none` | none/columns_before_collapsing/rows_before_collapsing/columns_after_collapsing/rows_after_collapsing/columns_before_collapsing_excluding_unassigned/rows_before_collapsing_excluding_unassigned | 标准化方式｜Normalization method (default: none) |
| `--average` | `none` | none/across_records/across_group_members/across_used_group_members/maximum/minimum/minimum_across_records | 组内聚合方式｜Aggregation method within groups (default: none) |
| `--row-names-are-in-column` | — |  | 包含行名的列名或索引（仅经典TSV表格）｜Column name or index containing row names (for TSV tables only) |
| `--output-format` | `auto` | auto/BIOM/classical | 输出格式｜Output format (default: auto) |
| `-t, --threads` | `1` | int | 线程数｜Number of threads (default: 1) |
| `-f, --force` | — | store_true | 覆盖已存在的输出文件｜Overwrite existing output files |
| `-v, --verbose` | — | store_true | 详细输出｜Verbose output |

<!-- END PARAMS:auto -->

## 依赖 | Dependencies

- **FAPROTAX**: 功能注释数据库与 collapse_table.py 脚本
- **Python**: BIOM 格式支持(读取 BIOM 输入时)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
