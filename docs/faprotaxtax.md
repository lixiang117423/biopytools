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

## 依赖 | Dependencies

- **FAPROTAX**: 功能注释数据库与 collapse_table.py 脚本
- **Python**: BIOM 格式支持(读取 BIOM 输入时)

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
