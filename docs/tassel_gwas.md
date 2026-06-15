# 🌾 TASSEL GWAS 全基因组关联分析模块

**基于TASSEL的高效GWAS分析工具 | Efficient GWAS Analysis Tool Based on TASSEL**

## 📖 功能概述 | Overview

TASSEL GWAS模块是一个强大的全基因组关联分析工具，基于业界标准的TASSEL软件，支持GLM、MLM模型，提供自动表型识别、批量处理、群体结构校正等功能，适用于各种规模的GWAS研究项目。

## ✨ 主要特性 | Key Features

- **🎯 多模型支持**: 支持GLM、MLM和BOTH模型，满足不同分析需求
- **📋 自动表型处理**: 自动识别和批量处理表型文件中的所有表型
- **🧮 群体结构校正**: 内置PCA主成分分析和Kinship矩阵计算
- **⚡ 高性能并行**: 支持多表型并行处理，大幅提升分析效率
- **📊 质量控制**: 灵活的MAF、缺失率过滤参数
- **📈 标准输出**: 生成曼哈顿图输入文件和详细分析报告

## 🚀 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 基本GWAS分析 - 自动识别并处理所有表型
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results

# 指定模型和参数
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \
    --model MLM --memory 200g --maf 0.05 --miss 0.1

# 使用Q矩阵校正群体结构
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \
    --q-matrix Q.txt

# 并行处理多个表型
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \
    --parallel --workers 8
```

### 高级用法 | Advanced Usage

```bash
# 指定PCA主成分数量（默认5个）
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \
    --pca-components 10

# 运行GLM和MLM两种模型
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \
    --model BOTH --keep-temp

# 高性能分析配置
biopytools tassel-gwas -i large_dataset.vcf.gz -p traits.txt -o results \
    --model MLM --memory 400g --threads 32 --parallel --workers 16
```

## 📋 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --vcf` | VCF基因型文件路径 | `-i population.vcf.gz` |
| `-p, --pheno` | 表型文件路径 | `-p traits.txt` |
| `-o, --output` | 输出目录路径 | `-o gwas_results` |

### 分析参数 | Analysis Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--model` | `MLM` | 🧬 GWAS模型选择 (GLM/MLM/BOTH) |
| `--memory` | `100g` | 💾 Java最大内存设置 |
| `-t, --threads` | `4` | 🧵 并行线程数 |
| `--pca-components` | `5` | 🧮 PCA主成分数量（用作MLM协变量） |

### 质量控制参数 | Quality Control Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `--maf` | 最小等位基因频率过滤 | `--maf 0.05` |
| `--miss` | 最大缺失率过滤 | `--miss 0.1` |
| `--t, --threads` | 并行线程数 | `--threads 16` |

### 矩阵文件 | Matrix Files

| 参数 | 描述 | 示例 |
|------|------|------|
| `--q-matrix` | 群体结构Q矩阵文件 | `--q-matrix Q.txt` |
| `--kinship` | 亲缘关系K矩阵文件 | `--kinship K.txt` |

### 其他选项 | Other Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `--tassel-path` | `auto` | 🔧 TASSEL安装路径 |
| `--skip-sort` | `False` | ⏭️ 跳过VCF排序 |
| `--keep-temp` | `False` | 💾 保留临时文件 |
| `--parallel` | `False` | 🚀 并行处理多个表型 |
| `--workers` | `4` | 👥 并行工作线程数 |

## 📁 输入文件格式 | Input File Formats

### VCF文件要求 | VCF File Requirements

支持标准VCF格式文件（压缩或未压缩）：

```vcf
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3
1	14370	rs6054257	G	A	29	PASS	.	GT	0/0	1/0	1/1
1	17330	.	T	A	3	q10	.	GT	0/0	0/1	0/0
1	1110696	rs6040355	A	G,T	67	PASS	.	GT	1/2	2/1	2/2
```

### 表型文件格式 | Phenotype File Format

表型文件应为制表符分隔的文本文件：

```
<Trait>	Trait1	Trait2	Trait3
<Data>	<Data>	<Data>	<Data>
Sample1	12.5	0.85	1.23
Sample2	15.2	0.92	1.45
Sample3	11.8	0.78	1.12
```

**文件要求**:
- 第一列为样本ID
- 第一行为表头，第一列应为`<Trait>`
- 第二行为数据类型标识，第一列应为`<Data>`
- 支持任意数量的表型列
- 样本ID应与VCF文件中的样本顺序一致

## 💡 使用示例 | Usage Examples

### 示例1：基础单表型分析 | Example 1: Basic Single Trait Analysis

```bash
# 分析单一表型
biopytools tassel-gwas \
    -i genotype_data.vcf.gz \
    -p single_trait.txt \
    -o single_trait_gwas \
    --model MLM
```

### 示例2：多表型批量分析 | Example 2: Multiple Traits Batch Analysis

```bash
# 自动识别并处理所有表型
biopytools tassel-gwas \
    -i population_genotypes.vcf.gz \
    -p multiple_traits.txt \
    -o multi_trait_gwas \
    --model MLM \
    --maf 0.05 \
    --miss 0.1
```

### 示例3：高性能并行分析 | Example 3: High-Performance Parallel Analysis

```bash
# 使用并行处理加速多表型分析
biopytools tassel-gwas \
    -i large_cohort.vcf.gz \
    -p phenotypes.txt \
    -o parallel_gwas \
    --model MLM \
    --parallel \
    --workers 16 \
    --threads 32 \
    --memory 400g
```

### 示例4：群体结构校正分析 | Example 4: Population Structure Corrected Analysis

```bash
# 使用PCA协变量和Kinship矩阵
biopytools tassel-gwas \
    -i genotypes.vcf.gz \
    -p traits.txt \
    -o corrected_gwas \
    --model MLM \
    --pca-components 10 \
    --maf 0.02 \
    --miss 0.05
```

### 示例5：双模型对比分析 | Example 5: Dual Model Comparison Analysis

```bash
# 同时运行GLM和MLM模型
biopytools tassel-gwas \
    -i genotype_data.vcf.gz \
    -p phenotypes.txt \
    -o model_comparison \
    --model BOTH \
    --keep-temp \
    --maf 0.05
```

### 示例6：使用预计算矩阵 | Example 6: Using Pre-computed Matrices

```bash
# 使用外部计算的群体结构矩阵
biopytools tassel-gwas \
    -i genotypes.vcf.gz \
    -p traits.txt \
    -o external_matrix_gwas \
    --model MLM \
    --q-matrix population_structure.txt \
    --kinship kinship_matrix.txt
```

## 📁 输出文件说明 | Output Files Description

### 主要输出文件 | Main Output Files

每个表型都会在输出目录下创建单独的文件夹，包含：

```
trait_name/                           # 表型专用目录
├── trait_name_GWAS.mlm.manht_input   # 曼哈顿图输入文件
├── trait_name_GWAS.pipeline.log       # TASSEL运行日志
├── trait_name_GWAS.pheno.txt         # 提取的表型文件
├── trait_name.pheno.txt              # 原始表型数据
└── ...
```

### 批量处理报告 | Batch Processing Reports

```
output_directory/
├── batch_processing_report.txt       # 批量处理汇总报告
├── failed_traits.log                # 失败表型记录
├── gwas.log                         # 全局日志文件
├── temp_matrix_calculation.*        # 临时矩阵文件（如果启用）
└── [每个表型的结果目录]               # Individual trait results
```

### 曼哈顿图输入格式 | Manhattan Plot Input Format

生成的`.manht_input`文件格式如下：

```
Chromosome	Position	Marker	Allele1	Allele2	Freq1	Freq2	P-value
1	14370	rs6054257	G	A	0.85	0.15	2.1e-08
1	17330	.	T	A	0.92	0.08	1.5e-06
...
```

## 🔧 系统要求 | System Requirements

### 依赖软件 | Dependencies

- **TASSEL** (版本 5.2+)
  - 下载地址: https://www.maizegenetics.net/tassel
- **Java** (版本 8+)
  - 用于运行TASSEL
- **VCF2PCACluster** (版本 1.42+)
  - 用于PCA计算（可选，如使用外部Q矩阵则不需要）
- **Python** (版本 3.7+)
- **Python包**:
  - `pandas` - 数据处理
  - `numpy` - 数值计算

### 安装依赖软件 | Installing Dependencies

```bash
# 下载并安装TASSEL
wget https://www.maizegenetics.net/tassel/download/tassel-5-standalone.zip
unzip tassel-5-standalone.zip
export TASSEL_PATH=/path/to/tassel-5-standalone/run_pipeline.pl

# 安装VCF2PCACluster（如需要）
# 按照官方文档安装：http://faculty.washington.edu/browning/vcf2pca.html

# 安装Python包
pip install pandas numpy click
```

### 硬件建议 | Hardware Recommendations

- **CPU**: 多核处理器（推荐8核以上）
- **RAM**: 最少16GB（大数据集推荐64GB以上）
- **存储**: 至少预留数据集大小5倍的磁盘空间
- **Java内存**: 根据数据规模调整（通常50g-400g）

## ⚠️ 注意事项 | Important Notes

1. **数据质量**: 输入VCF和表型文件质量直接影响GWAS分析结果
2. **样本对应**: 确保VCF和表型文件中的样本顺序一致
3. **内存配置**: 根据数据规模适当调整Java内存设置
4. **模型选择**: MLM模型更适合存在群体结构的数据
5. **PCA数量**: PCA主成分数量应根据数据特性选择，通常3-10个
6. **并行处理**: 并行处理会占用更多系统资源，需权衡速度与资源使用

## 🐛 故障排除 | Troubleshooting

### 常见问题 | Common Issues

**Q: "TASSEL command not found" 错误**
```bash
# 检查TASSEL路径
which perl
ls /path/to/tassel/run_pipeline.pl

# 确保perl可执行
chmod +x /path/to/tassel/run_pipeline.pl
```

**Q: 内存不足错误**
```bash
# 减少Java内存分配
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \
    --memory 50g

# 或增加系统内存/使用更大内存的机器
```

**Q: 样本数量不匹配**
```bash
# 检查VCF和表型文件的样本数量
vcf_samples=$(bcftools query -l input.vcf.gz | wc -l)
pheno_samples=$(tail -n +3 traits.txt | wc -l)
echo "VCF samples: $vcf_samples, Phenotype samples: $pheno_samples"
```

**Q: 某些表型分析失败**
```bash
# 检查失败的表型日志
cat output_directory/failed_traits.log

# 查看具体错误信息
cat output_directory/trait_name/trait_name_GWAS.pipeline.log
```

**Q: PCA计算失败**
```bash
# 检查VCF2PCACluster安装
ls /path/to/VCF2PCACluster-1.42/bin/VCF2PCACluster

# 或使用外部PCA文件
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \
    --q-matrix external_pca.txt
```

**Q: 分析速度慢**
```bash
# 启用并行处理
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \
    --parallel --workers 8

# 或增加线程数
biopytools tassel-gwas -i input.vcf.gz -p traits.txt -o results \
    --threads 16
```

## 📊 结果解读指南 | Result Interpretation Guide

### 曼哈顿图解读 | Manhattan Plot Interpretation

- **显著阈值**: 通常使用P < 5×10⁻⁸作为全基因组显著性阈值
- **QQ图**: 检查群体膨胀，λ值接近1表示无显著群体结构
- **Peak区域**: 显著SNP聚集的区域可能是候选基因区域

### 模型比较 | Model Comparison

- **GLM**: 简单快速，但容易假阳性
- **MLM**: 控制群体结构，结果更可靠
- **BOTH**: 建议同时运行，比较结果的一致性

## 📚 相关资源 | Related Resources

- [TASSEL官方文档](https://www.maizegenetics.net/tassel/)
- [GWAS分析方法综述](https://www.nature.com/articles/nrg3769)
- [群体结构对GWAS的影响](https://www.nature.com/articles/ng0608-713)
- [PCA在GWAS中的应用](https://www.nature.com/articles/ng0211-548)

## 📄 许可证 | License

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 🔬 引用信息 | Citation

如果在学术研究中使用此工具，请引用相关方法学文献：

```
Bradbury PJ, Zhang Z, Kroon DE, Casstevens TM, Ramdoss Y, Buckler ES.
(2007) TASSEL: software for association mapping of complex traits in diverse samples.
Bioinformatics 23:2633-2635.

Zhang Z, Ersoz E, Lai CQ, et al.
(2010) Mixed linear model approach adapted for genome-wide association studies.
Nat Genet 42:355-360.
```