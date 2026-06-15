# SNP Index 模块 | SNP Index Module

## 概述 | Overview

SNP Index模块是biopytools包的一个专门用于BSA（Bulked Segregant Analysis）分析的工具，可以从VCF文件计算SNP index和ΔSNP index，提供完整的统计分析、目标区域检测和可视化功能。

## 模块结构 | Module Structure

```
snp_index/
├── __init__.py          # 模块初始化和导出
├── config.py           # 配置管理
├── utils.py            # 工具函数
├── calculator.py       # SNP index计算核心
├── analyzer.py         # 结果分析
├── visualizer.py       # 数据可视化
├── main.py            # 主程序入口
└── README.md          # 本文件
```

## 核心类 | Core Classes

### SNPIndexConfig
配置管理类，用于管理所有分析参数。

### SNPIndexCalculator
SNP index计算器，负责从VCF文件计算SNP index和ΔSNP index。

### SNPIndexAnalyzer
结果分析器，提供统计分析和目标区域检测。

### SNPIndexVisualizer
可视化器，创建多种类型的分析图表。

### SNPIndexProcessor
主处理器，整合所有功能模块。

## 使用方法 | Usage

### 作为Python模块 | As Python Module

```python
from biopytools.snp_index import SNPIndexProcessor, SNPIndexConfig

# 创建配置
config = SNPIndexConfig(
    input_vcf="input.vcf.gz",
    output_dir="results",
    min_depth=10
)

# 运行分析
processor = SNPIndexProcessor(config)
processor.run_full_pipeline()
```

### 命令行使用 | Command Line

```bash
biopytools snp-index -i input.vcf.gz -o results/ -v
```

## 输出文件 | Output Files

- `{prefix}_results.tsv`: SNP index计算结果
- `{prefix}_extreme_sites.tsv`: 极端位点列表
- `{prefix}_potential_regions.tsv`: 潜在目标区域
- `{prefix}_comprehensive.png`: 综合分析图
- `{prefix}_manhattan.png`: 曼哈顿图

## 计算原理 | Calculation Principle

### SNP Index
```
SNP index = Alt_depth / (Ref_depth + Alt_depth)
```

### ΔSNP Index
```
ΔSNP index = SNP_index_pool1 - SNP_index_pool2
```

## 开发说明 | Development Notes

本模块遵循biopytools开发规范：

1. 统一的参数命名规范
2. 标准化的日志输出
3. 模块化的代码结构
4. 完整的错误处理
5. 详细的文档和注释

## 依赖项 | Dependencies

- numpy: 数值计算
- matplotlib: 数据可视化
- click: 命令行接口

## 作者 | Author

Xiang LI <xiang.li@your-domain.com>

## 版本 | Version

1.0.0 (2025-12-20)