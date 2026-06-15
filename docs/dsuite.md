# Dsuite - D统计量基因渗入分析工具

## 功能简介

使用Dsuite计算D统计量（ABBA-BABA test），检测群体间的基因渗入事件。该工具基于四群体检验（D统计量）和f4-ratio分析，可用于识别物种间的基因流和杂交事件。

## 系统要求

- Python >= 3.8
- Dsuite (已安装)
- bcftools

## Dsuite安装

```bash
# 克隆Dsuite仓库
git clone https://github.com/millanek/Dsuite.git
cd Dsuite
make

# Dsuite可执行文件位于 Build/Dsuite
```

## 使用方法

### 查看帮助

```bash
biopytools dsuite -h
```

### 基本用法

```bash
# 运行Dsuite Dtrios分析
biopytools dsuite \
    -i variants.vcf.gz \
    -s sets.txt \
    -o output_dir
```

### 参数说明

| 短参数 | 长参数 | 类型 | 必需 | 默认值 | 说明 |
|--------|--------|------|------|--------|------|
| `-i` | `--vcf` | str | 是 | - | 输入VCF文件路径 |
| `-s` | `--sets` | str | 是 | - | SETS分组文件路径 |
| `-o` | `--output` | str | 是 | - | 输出目录 |
| `-p` | `--prefix` | str | 否 | dsuite | 输出文件前缀 |
| `--min-alleles` | - | int | 否 | 2 | 最小等位基因数 |
| `--max-alleles` | - | int | 否 | 2 | 最大等位基因数 |
| `--variant-type` | - | str | 否 | snps | 变异类型 (snps/indels/both/none) |
| `--dsuite-bin` | - | str | 否 | /share/.../Dsuite | Dsuite可执行文件路径 |
| `--bcftools` | - | str | 否 | bcftools | bcftools命令路径 |

### SETS文件格式

SETS.txt文件定义每个样本所属的群体，格式为两列：`样本ID` 和 `群体ID`

```
Ind1    Species1
Ind2    Species1
Ind3    Species2
Ind4    Species2
Ind5    Species3
Ind6    Outgroup
Ind7    Outgroup
```

**重要说明：**
- 至少需要一个样本标记为 `Outgroup`
- 如果要忽略某些样本，使用 `xxx` 作为群体ID

## 输出文件

工具会生成以下文件：

1. **dsuite_BBAA.txt** - D统计量和f4-ratio结果
2. **dsuite_Dmin.txt** - D最小值结果
3. **dsuite_tree.txt** - 按树结构排列的结果
4. **dsuite_analysis_YYYYMMDD_HHMMSS.log** - 运行日志

## 结果解读

### D统计量

- **D值范围**: -1 到 1
- **D ≈ 0**: 无基因渗入
- **D显著偏离0**: 存在基因渗入信号
- **Z分数**: 用于评估统计显著性

### f4-ratio

- 估算基因渗入的比例
- 范围: 0 到 1
- 值越大表示渗入比例越高

## 使用场景

### 1. 基因渗入检测

```bash
# 检测物种间的基因渗入
biopytools dsuite \
    -i filtered_variants.vcf.gz \
    -s population_sets.txt \
    -o introgression_analysis
```

### 2. 系统发育关系验证

```bash
# 使用D统计量验证系统发育树
biopytools dsuite \
    -i variants.vcf.gz \
    -s species_sets.txt \
    -o phylogeny_test \
    -p tree_test
```

### 3. 杂交事件检测

```bash
# 检测杂交事件
biopytools dsuite \
    -i hybrid_variants.vcf.gz \
    -s hybrid_sets.txt \
    -o hybrid_detection \
    --min-alleles 2 --max-alleles 2
```

## 工作流程示例

### 完整的基因渗入分析流程

```bash
# 1. VCF质量控制
# bcftools filter -e 'QUAL>30 && DP>10' input.vcf.gz -o filtered.vcf.gz

# 2. 运行Dsuite分析
biopytools dsuite \
    -i filtered.vcf.gz \
    -s sets.txt \
    -o dsuite_output

# 3. 查看结果
cat dsuite_output/dsuite_BBAA.txt | head -20

# 4. 进一步分析（使用Dsuite Dinvestigate）
# Dsuite Dinvestigate filtered.vcf.gz sets.txt test_trios.txt
```

## 分析步骤说明

### 步骤1: 环境检查
- 检查bcftools和Dsuite是否可用
- 验证输入文件存在
- 检查VCF文件格式

### 步骤2: VCF统计
- 统计样本数量
- 统计总变异数
- 统计染色体/scaffold数量

### 步骤3: 应用过滤条件
- 按等位基因数过滤
- 按变异类型过滤
- 统计过滤后的变异数

### 步骤4: 运行Dsuite Dtrios
- 使用bcftools过滤VCF
- 通过管道传递给Dsuite
- 计算所有可能的三元组合

### 步骤5: 结果验证
- 检查输出文件是否生成
- 统计分析的三元组数量
- 显示分析总结

## 注意事项

1. **VCF格式**: 输入必须是有效的VCF格式（支持gzip压缩）
2. **样本命名**: VCF中的样本名必须与SETS文件中的样本名一致
3. **Outgroup**: 必须指定至少一个外群样本
4. **双等位位点**: 默认只分析双等位SNP，这是D统计量的标准做法
5. **内存需求**: 对于大型数据集，可能需要较多内存

## 常见问题

### Q: 如何解读D值？
A: D值接近0表示无基因渗入，显著偏离0（通常|Z|>3）表示存在基因渗入。正负号表示渗入的方向。

### Q: f4-ratio和D值的区别？
A: D值检测是否存在渗入，f4-ratio估算渗入的比例。f4-ratio范围在0-1之间。

### Q: 可以分析全基因组数据吗？
A: 可以，但建议按染色体分割并行分析以提高效率。

### Q: 如何选择过滤参数？
A: 通常使用双等位SNP（默认设置）。如需分析indels，可修改--variant-type参数。

## 技术细节

### 实现原理

1. **VCF过滤**: 使用bcftools过滤VCF文件
2. **管道传输**: 通过Unix管道将过滤后的VCF传递给Dsuite
3. **Dtrios分析**: Dsuite计算所有可能的三元组合的D统计量
4. **Jackknife**: 使用Jackknife方法计算标准误和Z分数

### 性能优化

- 使用管道避免中间文件
- 流式处理减少内存占用
- 支持并行分析（使用DtriosParallel脚本）

## 参考文献

**Malinsky, M., Matschiner, M. and Svardal, H. (2021)**
Dsuite - fast D-statistics and related admixture evidence from VCF files.
*Molecular Ecology Resources* 21, 584–595.
doi: [https://doi.org/10.1111/1755-0998.13265](https://doi.org/10.1111/1755-0998.13265)

## 许可证

MIT License

## 作者信息

**李详 (Xiang Li)**
- Email: lixiang117423@gmail.com
- GitHub: [@lixiang117423](https://github.com/lixiang117423)

## 相关工具

- [VCF系统发育分析](./vcf_phylo.md) - VCF系统发育树构建
- [Admixture](./admixture.md) - 群体结构分析
- [TASSEL GWAS](./tassel_gwas.md) - 全基因组关联分析
