# HapHiC Pipeline模式使用指南

## 概述

本工具现在使用 **HapHiC Pipeline模式** 一步完成基因组scaffolding，相比之前的分步运行模式，具有以下优势：

- **更高效率**：减少文件I/O开销，避免多次进程启动
- **更好一致性**：使用HapHiC原生pipeline，确保参数传递的一致性
- **更简操作**：无需关心中间步骤，一键完成整个流程
- **更强鲁棒性**：原生pipeline有更好的错误处理和参数验证

## 工作流程

```
输入文件 (FASTA + FASTQ/BAM)
        ↓
    BWA比对 (如需要)
        ↓
  HapHiC Pipeline
    ├── cluster (聚类)
    ├── reassign (重新分配)
    ├── sort (排序和定向) - 默认使用FAST sorting，跳过ALLHiC优化
    └── build (构建scaffolds)
        ↓
   输出文件 (FASTA/AGP)
        ↓
  Juicebox文件生成 (可选)
        ↓
   最终结果 (.hic + .assembly)
```

## 基本用法

### 最简单的用法
```bash
python -m biopytools.haphic main genome.fa hic.bam 12
```

### 完整参数示例
```bash
python -m biopytools.haphic main \
    genome.fa \
    hic_R1.fastq.gz \
    12 \
    -o output_dir \
    --prefix sample_name \
    --threads 64 \
    --processes 32 \
    --correct-nrounds 2 \
    --RE GATC \
    --min-inflation 1.0 \
    --max-inflation 3.0 \
    --Nx 80 \
    --no-juicebox
```

## 主要参数说明

### 必需参数
- `asm_file`: 基因组组装文件 (FASTA格式)
- `hic_file`: Hi-C数据文件 (FASTQ或BAM格式)
- `nchrs`: 预期染色体数量

### 输出配置
- `-o/--output-dir`: 输出目录路径
- `--prefix`: 输出文件前缀

### 性能参数
- `--threads`: 读取BAM文件的线程数 (默认: 8)
- `--processes`: 并行进程数 (默认: 8)

### 组装校正
- `--correct-nrounds`: 组装校正轮数 (默认: 0, 禁用)

### 聚类参数
- `--min-inflation`: 最小膨胀参数 (默认: 1.0)
- `--max-inflation`: 最大膨胀参数 (默认: 3.0)
- `--inflation-step`: 膨胀参数步长 (默认: 0.2)
- `--Nx`: Nx参数 (默认: 80)

### 其他选项
- `--RE`: 限制性内切酶位点 (默认: GATC)
- `--no-juicebox`: 不生成Juicebox文件
- `--quick-view`: 快速查看模式
- `--verbose`: 详细输出模式

## 输出文件

### 主要输出
- `{prefix}.scaffolds.fa`: 最终scaffolds序列
- `{prefix}.scaffolds.agp`: scaffolds的AGP文件
- `{prefix}.scaffolds.raw.agp`: 原始AGP文件

### Juicebox输出 (如果启用)
- `{prefix}.hic`: Juicebox格式的Hi-C接触矩阵
- `scaffolds.assembly`: 用于Juicebox的assembly文件

### 日志文件
- `{prefix}_haphic.log`: 完整的运行日志

## 示例用法

### 1. 基本scaffolding
```bash
python -m biopytools.haphic main genome.fa hic.bam 24
```

### 2. 指定输出目录和前缀
```bash
python -m biopytools.haphic main genome.fa hic.bam 24 \
    -o /path/to/output \
    --prefix sample1
```

### 3. 使用FASTQ输入和高性能设置
```bash
python -m biopytools.haphic main genome.fa hic_R1.fastq.gz 24 \
    --threads 64 \
    --processes 32 \
    --correct-nrounds 2
```

### 4. Arima试剂盒 (多种酶切位点)
```bash
python -m biopytools.haphic main genome.fa hic.bam 24 \
    --RE "GATC,GANTC"
```

### 5. 快速模式 (跳过组装校正)
```bash
python -m biopytools.haphic main genome.fa hic.bam 24 \
    --correct-nrounds 0
```

## 性能建议

### 1. 线程设置
- `--threads`: 建议设置为CPU核心数的一半
- `--processes`: 建议设置为CPU核心数的1/4到1/8

### 2. 内存要求
- 小基因组 (<1 Gb): 最少16GB内存
- 大基因组 (>1 Gb): 建议64GB+内存

### 3. 存储空间
- 需要至少3倍于基因组大小的存储空间

### 4. ALLHiC优化
- **默认跳过ALLHiC优化**（使用`--skip_allhic`）
- FAST sorting算法通常足够好且更稳定
- 如需启用ALLHiC优化，使用`--allhic-optimization`参数
- ALLHiC优化对简单group（contig数<50）提升有限

## 故障排除

### 常见问题

1. **BWA比对失败**
   - 检查FASTQ文件格式是否正确
   - 确认BWA工具在PATH中
   - 检查磁盘空间是否充足

2. **HapHiC Pipeline失败**
   - 检查Hi-C数据质量
   - 调整`--nchrs`参数
   - 尝试不同的膨胀参数范围

3. **Juicebox文件生成失败**
   - 确认3D-DNA工具安装正确
   - 检查路径配置
   - 使用`--no-juicebox`跳过此步骤

### 日志分析
查看详细的日志文件来诊断问题：
```bash
tail -f {prefix}_haphic.log
```

## 与旧版本的兼容性

- 旧版本的分步运行参数仍然被接受，但会显示警告并使用pipeline模式
- 所有原有的输出文件格式保持不变
- 参数映射已自动处理

## 技术细节

### Pipeline模式 vs 分步模式

| 特性 | Pipeline模式 | 分步模式 (旧) |
|------|-------------|---------------|
| 执行方式 | 一步完成 | 分多步执行 |
| 效率 | 更高 | 较低 |
| 参数传递 | 内部优化 | 需手动指定 |
| 错误处理 | 统一管理 | 分步处理 |
| 断点续传 | 不支持 | 支持 |
| 调试友好性 | 一般 | 更好 |

### 参数映射

新版本会自动将旧参数映射到pipeline参数：

- `--min-links` → pipeline内部默认值
- `--min-link-density` → pipeline内部默认值
- `--min-RE-sites` → `--RE_site_cutoff`
- 其他参数保持一致的映射关系