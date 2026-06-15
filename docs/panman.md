# PanMAN 泛基因组分析模块

**专业的泛基因组突变注释网络工具 | Professional Pangenome Mutation-Annotated Network Tool**

## 功能概述 | Overview

PanMAN 泛基因组分析模块是基于 PanMAN (Pangenome Mutation-Annotated Network) 的强大工具，提供构建和分析泛基因组网络的完整功能。PanMAN 是一种新颖的泛基因组数据表示方法，通过突变注释树 (PanMATs) 组成的网络，实现极高的存储压缩率 (2.9-559倍) 和强大的代表性。

## 主要特性 | Key Features

### 核心功能
- **🧬 PanGraph 生成**: 从 FASTA 文件生成 PanGraph JSON 和 Newick 树
- **🔨 灵活的构建方式**: 支持 PanGraph/GFA/MSA 三种输入格式
- **📤 多格式数据提取**: 支持 16+ 种数据格式提取
- **🌳 高效的压缩存储**: 提供 2.9-559 倍的压缩率
- **🧬 完整的变异注释**: 支持核苷酸和结构变异

### 高级功能
- **🔍 子网络提取**: 基于节点列表提取 PanMAN 子网络
- **📝 节点注释**: 自定义注释 PanMAN 节点
- **🌲 重新扎根**: 基于参考序列重新扎根系统发育树
- **🌐 网络创建**: 合并多个 PanMAN 创建网络
- **📊 突变打印**: 打印详细的突变信息
- **🧬 范围查询**: 提取特定坐标范围的序列
- **🎯 多 ACR 算法**: 支持 Fitch 和 MPPA 祖先重建方法

### 技术特性
- **⚙️ 多后端支持**: 支持 Conda 和 Docker 两种运行后端
- **📊 详细日志记录**: 完整的处理过程追踪
- **🔧 易于集成**: 符合 biopytools 开发规范，易于使用和扩展

## 快速开始 | Quick Start

### 完整工作流程 | Complete Workflow

从原始 FASTA 序列到 PanMAN 分析的完整流程：

```bash
# 步骤 1: 安装依赖
# PanMAN (用于构建和分析)
conda create -n panman python=3.11 -y
conda activate panman
conda install panman -y

# PanGraph (用于生成 PanGraph JSON)
# 下载: https://github.com/neherlab/pangraph/releases
# 或设置环境变量指向已安装的 PanGraph
export PANGRAPH_PATH=/path/to/pangraph

# 步骤 2: 从 FASTA 生成 PanGraph JSON
biopytools panman generate-pangraph \
    -i input_sequences.fa \
    -o my_dataset \
    --output-dir ./pangraph_results \
    -t 8

# 步骤 3: 构建 PanMAN
biopytools panman build \
    -P pangraph_results/my_dataset.json \
    -N pangraph_results/my_dataset.nwk \
    -o my_panman

# 步骤 4: 提取数据进行分析
biopytools panman extract \
    -I build/my_panman.panman \
    --summary \
    --vcf \
    --reference "reference_sequence_name"
```

### 安装 PanMAN | Install PanMAN

```bash
# 创建 Conda 环境
conda create -n panman python=3.11 -y
conda activate panman

# 设置 channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# 安装 PanMAN
conda install panman -y
```

### 安装 PanGraph (可选) | Install PanGraph (Optional)

```bash
# 方法 1: 下载预编译版本 (推荐)
wget https://github.com/neherlab/pangraph/releases/latest/download/pangraph-linux-x86_64
chmod +x pangraph-linux-x86_64

# 方法 2: 设置环境变量
export PANGRAPH_PATH=/path/to/pangraph-linux-x86_64

# 方法 3: 或在命令中指定 --pangraph-path 参数
```

### 基本用法 | Basic Usage

```bash
# 模式 1: 生成 PanGraph (从 FASTA)
biopytools panman generate-pangraph -i input.fa -o output

# 模式 2: 构建 PanMAN (从 PanGraph/GFA/MSA)
biopytools panman build -P input.json -N tree.nwk -o my_panman

# 模式 3: 提取数据 (从 PanMAN)
biopytools panman extract -I data.panman --summary --vcf --reference="seq1"
```

## 参数说明 | Parameters

### PanGraph 生成模式参数 | PanGraph Generation Mode Parameters

| 参数 | 简写 | 类型 | 必需 | 说明 |
|------|------|------|------|------|
| `--fasta` | `-i` | PATH | **是** | 输入 FASTA 文件路径 |
| `--output-prefix` | `-o` | STR | 否 | 输出文件前缀 (默认: output) |
| `--output-dir` | | STR | 否 | 输出目录 (默认: ./panman_output) |
| `--threads` | `-t` | INT | 否 | 线程数 (默认: 8) |
| `--pangraph-path` | | PATH | 否 | PanGraph 可执行文件路径 |

**输出文件**:
- `{prefix}.json` - PanGraph JSON 文件
- `{prefix}.nwk` - Newick 系统发育树
- `{prefix}_pangraph.log` - PanGraph 运行日志

### 构建模式参数 | Build Mode Parameters

| 参数 | 简写 | 类型 | 必需 | 说明 |
|------|------|------|------|------|
| `--pangraph` | `-P` | PATH | 条件 | PanGraph JSON 文件 |
| `--gfa` | `-G` | PATH | 条件 | GFA 文件 |
| `--msa` | `-M` | PATH | 条件 | MSA 文件 (FASTA 格式) |
| `--newick` | `-N` | PATH | **是** | Newick 系统发育树文件 |
| `--output-prefix` | `-o` | STR | 否 | 输出文件前缀 (默认: output) |
| `--output-dir` | | STR | 否 | 输出目录 (默认: ./panman_output) |
| `--backend` | | STR | 否 | 后端选择 (conda/docker, 默认: conda) |
| `--threads` | `-t` | INT | 否 | 线程数 (默认: 8) |

**注意**: 构建模式需要提供 `-P`/`-G`/`-M` 中的至少一种输入格式

### 提取模式参数 | Extract Mode Parameters

#### 基础参数 | Basic Parameters

| 参数 | 简写 | 类型 | 必需 | 说明 |
|------|------|------|------|------|
| `--input-panman` | `-I` | PATH | **是** | PanMAN 文件路径 |
| `--output-prefix` | `-o` | STR | 否 | 输出文件前缀 (默认: output) |
| `--output-dir` | | STR | 否 | 输出目录 (默认: ./panman_output) |
| `--reference` | `-r` | STR | 条件 | 参考序列名称 (VCF/重新扎根需要) |
| `--backend` | | STR | 否 | 后端选择 (默认: conda) |
| `--threads` | `-t` | INT | 否 | 线程数 (默认: 8) |

#### 数据提取选项 | Data Extraction Options

| 参数 | 类型 | 说明 |
|------|------|------|
| `--summary` | FLAG | 提取摘要统计 |
| `--extract-fasta` | FLAG | 提取 FASTA 序列 |
| `--extract-msa` | FLAG | 提取 MSA 比对 |
| `--vcf` | FLAG | 提取 VCF 变异 (需要 --reference) |
| `--extract-gfa` | FLAG | 提取 GFA 格式 |
| `--extract-newick` | FLAG | 提取 Newick 树 |
| `--extended-newick` | FLAG | 提取扩展 Newick 格式 |
| `--maf` | FLAG | 提取 MAF 格式 |
| `--aa` | FLAG | 提取氨基酸翻译 |

#### 高级功能选项 | Advanced Feature Options

| 参数 | 类型 | 说明 |
|------|------|------|
| `--subnet` | FLAG | 提取子网络 (需要 --input-file) |
| `--annotate` | FLAG | 注释节点 (需要 --input-file) |
| `--reroot` | FLAG | 重新扎根树 (需要 --reference) |
| `--create-network` | FLAG | 创建网络 (需要 --input-file) |
| `--print-mutations` | FLAG | 打印突变信息 |

#### 高级参数 | Advanced Parameters

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--pangraph-path` | PATH | - | PanGraph 可执行文件路径 |
| `--tree-id` | STR | - | 树 ID (VCF 提取可选) |
| `--acr` | STR | fitch | ACR 方法 (fitch/mppa) |
| `--input-file` | PATH | - | 输入文件 (用于 subnet/annotate/create-network) |
| `--range-index` | STR | - | 范围查询 index 参数 |
| `--range-start` | INT | - | 范围查询起始坐标 |
| `--range-end` | INT | - | 范围查询结束坐标 |

## 使用示例 | Usage Examples

### 示例 1: 从 FASTA 生成 PanGraph | Example 1: Generate PanGraph from FASTA

```bash
# 基本用法
biopytools panman generate-pangraph \
    -i input_sequences.fa \
    -o my_dataset

# 指定输出目录和线程数
biopytools panman generate-pangraph \
    -i sars_20.fa \
    -o sars_20 \
    --output-dir ./pangraph_results \
    -t 16

# 使用自定义 PanGraph 路径
biopytools panman generate-pangraph \
    -i genomes.fa \
    -o output \
    --pangraph-path /share/org/YZWL/yzwl_lixg/software/pangraph/pangraph-x86_64-unknown-linux-gnu
```

### 示例 2: 从 PanGraph 构建 | Example 2: Build from PanGraph

```bash
biopytools panman build \
    -P pangraph.json \
    -N tree.nwk \
    -o sars_covid \
    --output-dir ./panman_results
```

### 示例 3: 从 GFA 构建 | Example 3: Build from GFA

```bash
biopytools panman build \
    -G assembly.gfa \
    -N phylogeny.nwk \
    -o bacterial_panman
```

### 示例 4: 从 MSA 构建 | Example 4: Build from MSA

```bash
biopytools panman build \
    -M alignment.msa.fa \
    -N tree.nwk \
    -o viral_panman
```

### 示例 5: 提取摘要和 FASTA | Example 5: Extract Summary and FASTA

```bash
biopytools panman extract \
    -I data.panman \
    --summary \
    --extract-fasta \
    -o extracted_data
```

### 示例 6: 提取 VCF (需要参考序列) | Example 6: Extract VCF (requires reference)

```bash
biopytools panman extract \
    -I data.panman \
    --vcf \
    --reference "Switzerland/SO-ETHZ-500145/2020|OU000199.2|2020-11-12" \
    -o variants
```

### 示例 7: 批量提取多种格式 | Example 7: Extract Multiple Formats

```bash
biopytools panman extract \
    -I data.panman \
    --summary \
    --extract-fasta \
    --extract-msa \
    --extract-newick \
    --vcf \
    --reference "reference_sequence" \
    -o full_analysis
```

### 示例 8: 高级功能 - 子网络提取 | Example 8: Advanced - Subnet Extraction

```bash
# 首先创建节点列表文件
cat > nodes.txt << EOF
seq1
seq5
seq10
EOF

# 提取子网络
biopytools panman extract \
    -I data.panman \
    --subnet \
    --input-file nodes.txt \
    -o subnet_data
```

### 示例 9: 高级功能 - 节点注释 | Example 9: Advanced - Node Annotation

```bash
# 创建注释文件 (TSV 格式)
cat > annotations.tsv << EOF
seq1\tAsia
seq2\tEurope
seq3\tNorth_America
EOF

# 注释节点
biopytools panman extract \
    -I data.panman \
    --annotate \
    --input-file annotations.tsv \
    -o annotated_data
```

### 示例 10: 高级功能 - 重新扎根树 | Example 10: Advanced - Reroot Tree

```bash
biopytools panman extract \
    -I data.panman \
    --reroot \
    --reference "new_reference_sequence" \
    --acr mppa \
    -o rerooted_data
```

### 示例 11: 高级功能 - 范围查询 | Example 11: Advanced - Range Query

```bash
biopytools panman extract \
    -I data.panman \
    --range-index no \
    --range-start 1000 \
    --range-end 2000 \
    --reference "reference_sequence" \
    -o range_data
```

### 示例 12: 使用 Docker 后端 | Example 12: Use Docker Backend

```bash
biopytools panman build \
    -P input.json \
    -N tree.nwk \
    --backend docker \
    -o output
```

### 示例 13: 完整工作流 | Example 13: Complete Workflow

```bash
# 步骤 1: 从 FASTA 生成 PanGraph
biopytools panman generate-pangraph \
    -i genomes.fa \
    -o my_analysis \
    --output-dir ./step1_pangraph

# 步骤 2: 构建 PanMAN
biopytools panman build \
    -P step1_pangraph/my_analysis.json \
    -N step1_pangraph/my_analysis.nwk \
    -o my_panman

# 步骤 3: 提取摘要统计
biopytools panman extract \
    -I build/my_panman.panman \
    --summary \
    -o summary

# 步骤 4: 提取 VCF 变异
biopytools panman extract \
    -I build/my_panman.panman \
    --vcf \
    --reference "reference_genome" \
    -o variants
```

## 输入文件格式 | Input File Formats

### PanGraph JSON 格式 | PanGraph JSON Format

PanGraph 是一种 JSON 格式的泛基因组图表示：

```json
{
  "sequences": {
    "seq1": {
      "id": "seq1",
      "name": "Sequence 1",
      "blocks": [
        {"blockId": "block1", "path": [["block1", 1]]}
      ]
    }
  },
  "blocks": {
    "block1": {
      "id": "block1",
      "length": 1000,
      "sequence": "ATCGATCG..."
    }
  }
}
```

### GFA 格式 | GFA Format

标准图格式组装 (Graphical Fragment Assembly) 格式：

```gfa
H\tVN:Z:1.0
S\tnode1\tATCGATCGATCG\tLN:i:12
S\tnode2\tGCTAGCTAGCTA\tLN:i:12
L\tnode1\t+\tnode2\t+\t12M
```

### MSA FASTA 格式 | MSA FASTA Format

多重序列比对 (Multiple Sequence Alignment) FASTA 格式：

```fasta
>seq1
ATCGATCG---ATCG
>seq2
ATCG---GATCGATCG
>seq3
ATCGATCGATCG---C
```

### Newick 树格式 | Newick Tree Format

系统发育树的标准文本表示：

```
((seq1:0.1,seq2:0.1):0.2,seq3:0.3);
```

## 输出结果 | Output Results

### 构建模式输出 | Build Mode Output

```
panman_output/
└── output.panman          # PanMAN 二进制文件
```

### 提取模式输出 | Extract Mode Output

根据选择的提取选项，输出不同的文件：

```
panman_output/
├── output_summary.txt          # 摘要统计
├── output.fa                   # FASTA 序列
├── output_msa.fa               # MSA 比对
├── output.vcf                  # VCF 变异
├── output.gfa                  # GFA 图
├── output.nwk                  # Newick 树
├── output_extended.nwk         # 扩展 Newick
├── output.maf                  # MAF 格式
├── output_aa.tsv               # 氨基酸翻译
└── panman_analysis.log         # 分析日志
```

## 输出格式说明 | Output Format Descriptions

### Summary (摘要统计)

包含 PanMAN 的统计信息：
- 节点数量
- 树的数量
- 突变统计
- 压缩率

### FASTA 序列

提取所有或指定序列的 FASTA 格式文件。

### MSA 比对

提取多重序列比对，包含插入缺失信息。

### VCF 变异

相对于指定参考序列的变异调用格式。

### GFA 图

将 PanMAN 转换为 GFA 图格式。

### Newick 树

提取系统发育树。

### MAF 格式

多基因组比对 (Multiple Whole Genome Alignment) 格式。

## 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **PanMAN** (v0.1.4 或更新)
  - Bioconda: `conda install panman`
  - 支持平台: linux-64, osx-64
- **Python** (3.7+)
  - `click` - 命令行界面
  - `pathlib` - 路径处理

### Conda 环境 | Conda Environment

```bash
# 默认环境名
panman_v.0.1.4

# 可在代码中修改
biopytools/panman/config.py 中的 conda_env 参数
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器 (推荐 4 核以上)
- **RAM**: 最少 4GB (大基因组推荐 16GB+)
- **存储**: 预留输入文件大小 2-3 倍的磁盘空间

## 注意事项 | Important Notes

1. **Conda 环境**: 使用前请确保已安装并配置好正确的 Conda 环境
2. **文件格式**: 确保输入文件符合标准格式要求
3. **参考序列**: VCF 提取必须指定参考序列名称
4. **内存使用**: 处理大型泛基因组时可能需要大量内存
5. **输出目录**: 程序会自动创建输出目录

## 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "Conda environment does not exist" 错误**

```bash
# 检查环境是否存在
conda env list

# 如果不存在，创建环境
conda create -n panman_v.0.1.4 python=3.11 -y
conda activate panman_v.0.1.4
conda install panman -y
```

**Q: "Command execution failed" 错误**

```bash
# 检查 PanMAN 是否正确安装
conda activate panman_v.0.1.4
panmanUtils --help

# 如果未安装
conda install panman -y
```

**Q: "Invalid input file" 错误**

```bash
# 验证输入文件格式
# 对于 PanGraph JSON:
python -m json.tool input.json

# 对于 FASTA:
head -n 1 input.fa  # 应该以 > 开头

# 对于 Newick:
tail -n 1 tree.nwk  # 应该以 ; 结尾
```

**Q: VCF 提取缺少参考序列**

```bash
# 必须指定 --reference 参数
biopytools panman extract \
    -I data.panman \
    --vcf \
    --reference "your_reference_sequence_name"
```

**Q: 内存不足错误**

```bash
# 减少线程数
biopytools panman build ... -t 4

# 或分批处理数据
```

## 高级用法 | Advanced Usage

### Python API 使用 | Python API Usage

```python
from biopytools.panman import PanMANBuildRunner, PanMANExtractRunner

# 构建 PanMAN
builder = PanMANBuildRunner(
    pangraph_file="input.json",
    newick_file="tree.nwk",
    output_prefix="my_panman",
    output_dir="./results"
)
success = builder.run_analysis()

# 提取数据
extractor = PanMANExtractRunner(
    panman_file="data.panman",
    output_prefix="output",
    extract_summary=True,
    extract_fasta=True,
    extract_vcf=True,
    reference="seq1"
)
success = extractor.run_analysis()
```

### 自定义配置 | Custom Configuration

```python
from biopytools.panman import PanMANConfig

config = PanMANConfig(
    mode="build",
    pangraph_file="input.json",
    newick_file="tree.nwk",
    conda_env="custom_panman_env",
    backend="conda",
    threads=16
)
config.validate()
```

## 性能优化建议 | Performance Optimization

1. **使用 SSD**: 将输入输出文件放在 SSD 上可显著提升速度
2. **调整线程数**: 根据 CPU 核心数调整 `--threads` 参数
3. **内存优化**: 处理大型数据集时增加可用内存
4. **批量提取**: 一次性提取多种格式比多次提取更高效

## 相关资源 | Related Resources

- [PanMAN 官方文档](https://turakhia.ucsd.edu/panman/)
- [PanMAN GitHub](https://github.com/TurakhiaLab/panman)
- [BioConda PanMAN Recipe](https://bioconda.github.io/recipes/panman/README.html)
- [Nature Genetics 论文](https://doi.org/10.1038/s41588-025-02478-7)

## 引用信息 | Citation

如果在学术研究中使用 PanMAN，请引用以下论文：

```
Walia S., Motwani H., Tseng YH., Smith K., Corbett-Detig R., Turakhia Y.,
Compressive pangenomics using mutation-annotated networks.
Nat Genet (2026). https://doi.org/10.1038/s41588-025-02478-7
```

## 许可证 | License

本项目采用 MIT 许可证 - 详见 [LICENSE](../LICENSE) 文件

**注意**: PanMAN 软件本身遵循其各自的许可证条款。

---

## 版本历史 | Version History

| 版本 | 日期 | 主要变更 |
|------|------|----------|
| 1.1.0 | 2026-01-15 | 新增 PanGraph 生成模式，支持从 FASTA 生成 PanGraph JSON；新增 7 项高级提取功能 (子网络/注释/重新扎根/创建网络/打印突变/ACR方法/范围查询)；修复 threads 参数传递；支持可配置 conda 路径 |
| 1.0.0 | 2026-01-15 | 初始版本，支持构建和提取功能 |

---

**开发团队 | Development Team**: BioPyTools Project
**最后更新 | Last Updated**: 2026-01-15
