# NGenomeSyn Pipeline - 基因组共线性分析工具

## 功能简介

使用NGenomeSyn进行基因组共线性可视化分析的简化工具。从FASTA文件开始，自动完成基因组比对（Minimap2/MUMmer）和可视化（NGenomeSyn）。

## 与 genomesyn 命令的区别

- **ngenomesyn**（本工具）：简化版本，自动完成比对和可视化，不支持SyRI
- **genomesyn**：完整版本，支持SyRI结构变异分析、配置文件生成等高级功能

**建议**：
- 如果需要简单的快速可视化，使用 `biopytools ngenomesyn`
- 如果需要SyRI结构变异分析或配置文件，使用 `biopytools genomesyn`

## 系统要求

- Python >= 3.8
- NGenomeSyn (已安装)
- Minimap2 或 MUMmer (已安装)
- Biopython (pip install biopython)
- ImageMagick (可选，用于PNG输出)

## 使用方法

### 查看帮助

```bash
biopytools ngenomesyn -h
```

### 基本用法

```bash
# 使用Minimap2比对（默认）
biopytools ngenomesyn -s samples.txt -o output_dir

# 使用MUMmer比对
biopytools ngenomesyn -s samples.txt -o output_dir --aligner mummer
```

### 参数说明

| 短参数 | 长参数 | 类型 | 必需 | 默认值 | 说明 |
|--------|--------|------|------|--------|------|
| `-s` | `--sample-map` | str | 是* | - | 样本映射文件 |
| `-o` | `--output` | str | 是 | - | 输出目录 |
| `-a` | `--aligner` | str | 否 | minimap2 | 比对器类型 (minimap2/mummer) |
| `-t` | `--threads` | int | 否 | 16 | 线程数 |
| `--min-length` | - | int | 否 | 5000 | 最小比对长度 |
| `--minimap-preset` | - | str | 否 | asm5 | Minimap2预设模式 |
| `--mummer-match-type` | - | str | 否 | mumreference | MUMmer匹配类型 |
| `--chromosome` | - | str | 否 | - | 指定染色体 (如: "1,2,3") |
| `--output-formats` | - | list | 否 | svg png | 输出格式 |

* sample-map 和 config 参数必须提供一个

### sample_map 文件格式

```
genome_file.fa    GenomeName
reference.fa     Ref
query.fa         Query
```

格式为两列，用Tab分隔：
- 第一列：基因组FASTA文件路径
- 第二列：基因组名称

## 使用示例

### 示例1：基本用法

```bash
biopytools ngenomesyn -s samples.txt -o ./output
```

### 示例2：使用MUMmer比对器

```bash
biopytools ngenomesyn -s samples.txt -o ./output --aligner mummer
```

### 示例3：分析特定染色体

```bash
biopytools ngenomesyn -s samples.txt -o ./output --chromosome "1,2,3"
```

### 示例4：自定义参数

```bash
biopytools ngenomesyn -s samples.txt -o ./output \\
    --threads 32 --min-length 10000
```

### 示例5：仅生成SVG格式

```bash
biopytools ngenomesyn -s samples.txt -o ./output --output-formats svg
```

## 输出文件

工具会自动生成以下文件：

```
output_dir/
├── Ref.len                    # 基因组长度文件
├── Query.len                  # 基因组长度文件
├── Ref_vs_Query.paf          # Minimap2比对结果 (或.coords for MUMmer)
├── Ref_vs_Query.link         # LINK格式比对结果
├── ngenomesyn.conf           # NGenomeSyn配置文件
├── genome_synteny.svg        # SVG可视化图
├── genome_synteny.png        # PNG可视化图
├── ngenomesyn_pipeline.log   # 运行日志
└── alignment_commands.sh     # 执行的命令脚本
```

## 工作流程

1. **环境检查**：验证NGenomeSyn和比对器是否可用
2. **读取样本**：从sample_map文件读取基因组信息
3. **执行比对**：使用Minimap2或MUMmer进行基因组比对
4. **格式转换**：将比对结果转换为LINK格式
5. **生成配置**：自动生成NGenomeSyn配置文件
6. **可视化**：运行NGenomeSyn生成可视化图形

## 注意事项

1. **样本文件格式**：必须是两列（基因组文件、基因组名称），用Tab分隔
2. **基因组文件**：支持FASTA和GZIP压缩的FASTA格式
3. **至少两个基因组**：共线性分析至少需要两个基因组
4. **比对器选择**：
   - Minimap2：速度快，适合快速可视化
   - MUMmer：更精确，适合精细分析
5. **染色体编号**：染色体编号从1开始，连续编号

## 常见问题

### Q: 与 genomesyn 命令有什么区别？
A:
- `ngenomesyn`：简化版本，自动完成比对和可视化
- `genomesyn`：完整版本，支持SyRI、配置文件生成等高级功能

### Q: 如何选择比对器？
A:
- **Minimap2**：速度快，适合大多数情况
- **MUMmer**：更精确，适合需要高精度的分析

### Q: 染色体编号如何确定？
A: 染色体编号按照FASTA文件中序列的顺序，从1开始。

### Q: PNG转换失败怎么办？
A: 请确保安装了ImageMagick，或者使用 `--output-formats svg` 仅生成SVG格式。

## 许可证

MIT License

## 作者信息

**李详 (Xiang Li)**
- Email: lixiang117423@gmail.com
- GitHub: [@lixiang117423](https://github.com/lixiang117423)

## 相关工具

- [genomesyn](./ngenomesyn.md) - 完整的基因组共线性分析工具（含SyRI）
- [minimap2](https://github.com/lh3/minimap2) - 全基因组比对工具
- [MCScanX](https://github.com/wyp1125/MCScanX) - 共线性分析工具
