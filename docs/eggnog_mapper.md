# eggNOG功能注释 | eggNOG Functional Annotation

**用eggnog-mapper对蛋白/CDS/基因组进行功能注释(GO/KEGG/COG/CAZy/Pfam) | Annotate proteins/CDS/genomes with eggnog-mapper (GO/KEGG/COG/CAZy/Pfam)**

## 功能概述 | Overview

eggnog_mapper 模块封装了 [eggnog-mapper](http://eggnog-mapper.embl.de/), 将输入序列映射到 eggNOG 直系同源组, 获得功能注释(GO、KEGG 通路、COG 功能类别、CAZy、Pfam 结构域等)。支持多种输入类型(蛋白/CDS/基因组/宏基因组)与多种搜索模式(mmseqs/diamond/hmmer), 并将原生产物重排版为整洁表格。

## 快速开始 | Quick Start

```bash
# 蛋白输入, 默认 mmseqs 模式
biopytools eggnog-mapper -i proteins.faa -o out/

# CDS 输入(翻译为蛋白再注释)
biopytools eggnog-mapper -i cds.fa -o out/ --itype CDS --translate

# 指定数据目录与 diamond 模式
biopytools eggnog-mapper -i proteins.faa -o out/ -m diamond --data-dir ~/database/eggnog
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入 FASTA(蛋白/CDS/基因组) |
| `-o, --output-dir` | 输出目录 |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--itype` | `proteins` | 输入类型(proteins/CDS/genome/metagenome) |
| `--translate` | `False` | CDS 翻译为蛋白(itype=CDS/genome/metagenome) |
| `-m, --mode` | `mmseqs` | 搜索模式(mmseqs/diamond/hmmer/no_search/cache) |
| `--cpu` | `12` | 线程数 |
| `--sensmode` | `sensitive` | 灵敏度 |
| `--seed-ortholog-evalue` | `0.001` | seed ortholog E 值 |
| `--data-dir` | `~/database/eggnog` | eggNOG 数据库目录 |
| `--prefix` | 输入文件名 | 输出前缀 |
| `--emapper-path` | 自动检测 | emapper.py 路径 |
| `--resume` | `False` | 续传 |
| `--override` | `False` | 覆盖已有输出 |
| `--no-format` | `False` | 跳过重排版, 只留原生产物 |

(运行 `biopytools eggnog-mapper -h` 查看完整参数列表)

## 输出 | Output

- eggnog-mapper 原始注释产物
- 重排版整洁注释表(默认; `--no-format` 时跳过)
- 运行日志

## 依赖 | Dependencies

- **eggnog-mapper (emapper.py)**: 功能注释主程序
- **eggNOG 数据库**: 直系同源组与功能注释数据

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
