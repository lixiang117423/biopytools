# DeepBSA批量分析模块

**专业的BSA分析批量处理工具 | Professional BSA Analysis Batch Processing Tool**

## 功能概述 | Overview

DeepBSA批量分析模块是对DeepBSA工具的封装和增强，提供了批量运行多种BSA分析方法、自动处理复杂文件名、智能结果合并等功能。该模块解决了DeepBSA原始工具在文件名识别、并行执行、结果整合等方面的局限性，大幅提升了BSA分析的效率和易用性。

## 主要特性 | Key Features

- **7种BSA方法批量运行**: DL、K、ED4、SNP、SmoothG、SmoothLOD、Ridit一键执行
- **智能文件名识别**: 支持复杂VCF文件名（如`variation.filtered.snp.biallelic.vcf`）
- **方法级并行执行**: 所有方法同时运行，大幅节省时间
- **独立工作目录**: 每个方法独立目录，避免文件冲突
- **自动结果合并**: CSV结果、PNG图片自动整合，标注来源方法
- **DL方法修复**: 自动创建Models符号链接，解决模型文件依赖问题
- **详细日志记录**: 每个方法独立日志，便于问题排查
- **完整统计报告**: 自动生成各方法QTL统计汇总

## 源码修改说明 | Source Code Modifications

为了解决DeepBSA原始工具的各种问题，我们对DeepBSA源码进行了以下修改：

### 1. Pandas兼容性修复

**文件位置**: `/share/org/YZWL/yzwl_lixg/software/DeepBSA/DeepBSA_linux_v1.4/bin/functions/vcf_handle.py`

**修改位置**: 第131行

**修改内容**:
```python
# 修改前（原代码）:
df.to_csv(save_path, header=False, index=False, sep=',', lineterminator='\n')

# 修改后（修复后）:
# 明确指定CSV格式（移除lineterminator以兼容旧版pandas）
# Specify CSV format (remove lineterminator for compatibility with older pandas)
df.to_csv(save_path, header=False, index=False, sep=',')
```

**修改原因**:
- `lineterminator`参数在较旧版本的pandas中不存在
- DeepBSA使用的pandas版本较旧，导致运行时报错`TypeError: to_csv() got an unexpected keyword argument 'lineterminator'`
- 移除此参数不影响CSV文件的正确性

**影响范围**: 所有使用`VCF2Excel`类的BSA方法

---

## 快速开始 | Quick Start

### 基本用法 | Basic Usage

```bash
# 运行所有7种BSA方法（并行）
biopytools deepbsa -i variant.vcf -o bsab_results

# 只运行指定的方法
biopytools deepbsa -i variant.vcf -m DL,K,ED4 -o bsab_results

# 串行运行（调试时使用）
biopytools deepbsa -i variant.vcf --no-parallel -o bsab_results
```

### 处理复杂文件名 | Handling Complex Filenames

```bash
# 自动处理复杂文件名（如 variation.filtered.snp.biallelic.vcf）
biopytools deepbsa -i variation.filtered.snp.biallelic.vcf -o bsab_results
# 工具会自动创建符号链接 deepbsa_input.vcf → 实际文件
```

## 参数说明 | Parameters

### 必需参数 | Required Parameters

| 参数 | 描述 | 示例 |
|------|------|------|
| `-i, --input-file` | 输入VCF文件路径 | `-i variant.vcf` |

### 可选参数 | Optional Parameters

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-m, --methods` | `全部` | 要运行的方法，逗号分隔（DL,K,ED4,SNP,SmoothG,SmoothLOD,Ridit） |
| `-o, --output-dir` | `deepbsa_results` | 输出目录 |
| `-p, --parallel` | `True` | 并行运行所有方法（默认） |
| `--no-parallel` | `False` | 串行运行所有方法 |
| `-n, --no-auto-clean` | `False` | 不自动清理VCF注释行（已弃用，DeepBSA自己处理） |
| `-k, --keep-clean` | `False` | 保留清理后的文件（已弃用） |
| `--deepbsa-path` | `~/software/DeepBSA/DeepBSA_linux_v1.4/bin/main.py` | DeepBSA主程序路径 |
| `--conda-env` | `/share/org/YZWL/yzwl_lixg/miniforge3/envs/DeepBSA` | Conda环境路径 |
| `-v, --verbose` | `False` | 详细输出 |

## 可用方法 | Available Methods

| 方法 | 全称 | 描述 | 是否需要模型 |
|------|------|------|-------------|
| **DL** | Deep Learning | 基于深度学习的BSA分析方法 | ✅ 是 |
| **K** | K-Method | 基于K值的BSA分析方法 | ❌ 否 |
| **ED4** | Euclidean Distance 4 | 基于欧氏距离的BSA分析方法 | ❌ 否 |
| **SNP** | SNP-Index | SNP指数法 | ❌ 否 |
| **SmoothG** | Smooth G-value | 平滑G值法 | ❌ 否 |
| **SmoothLOD** | Smooth LOD | 平滑LOD值法 | ❌ 否 |
| **Ridit** | Ridit Analysis | Ridit分析法 | ❌ 否 |

## 输出文件结构 | Output Structure

```
output_dir/
├── deepbsa.log                          # 主日志文件
├── merged_results/                      # 自动生成的合并结果 ⭐
│   ├── merged_results.csv               # 所有方法的QTL结果（含Method列）
│   ├── images/                          # 所有方法的图片（带方法名前缀）
│   │   ├── DL_0-DL-LOWESS-auto-0.0771.png
│   │   ├── ED4_0-ED4-LOWESS-auto-0.0267.png
│   │   ├── K_0-K-LOWESS-auto-0.1300.png
│   │   └── ...
│   └── summary_report.txt               # 汇总统计报告
├── DL/                                  # DL方法的独立目录
│   ├── Models → .../DeepBSA/bin/Models # 符号链接（自动创建）⭐
│   ├── Results/variant/
│   │   ├── 0-DL-LOWESS-auto-0.0771.csv
│   │   └── 0-DL-LOWESS-auto-0.0771.png
│   ├── Excel_Files/variant.csv
│   ├── NoPretreatment/
│   ├── Pretreated_Files/
│   └── DL.log                          # 方法独立日志
├── K/
├── ED4/
├── SNP/
├── SmoothG/
├── SmoothLOD/
└── Ridit/
```

⭐ = 本模块新增或自动处理的功能

## 合并结果说明 | Merged Results Details

### CSV结果文件 (merged_results.csv)

包含以下列：
- **Method**: 来源方法（新增）
- **QTL**: QTL编号
- **Chr**: 染色体
- **Left**: 左侧位置
- **Peak**: 峰值位置
- **Right**: 右侧位置
- **Value**: QTL值
- **Source_File**: 原始文件名（新增）

**注意**: 只保留Value != "-"的有效QTL结果

### PNG图片文件

所有PNG图片从各方法目录复制到`merged_results/images/`，文件名添加方法名前缀：
- 原始: `0-DL-LOWESS-auto-0.0771.png`
- 合并后: `DL_0-DL-LOWESS-auto-0.0771.png`

### 汇总报告 (summary_report.txt)

包含每个方法的统计信息：
- 总QTL数
- 有效QTL数（Value != "-"）
- 最高值
- 最低值

## 工作流程 | Workflow

```
1. 参数解析与验证
   ↓
2. 输入文件预处理
   - 检测文件名复杂度
   - 如需，创建符号链接（deepbsa_input.vcf → 实际文件）
   ↓
3. 为每个方法创建独立工作目录
   - 创建方法目录（DL/, K/, ED4/, ...）
   - 创建Models符号链接（DL方法需要）⭐
   ↓
4. 并行执行所有方法（默认）
   - 使用conda run隔离环境
   - 每个方法在独立目录运行
   - 禁用DeepBSA的pretreatment（--p 0）
   ↓
5. 等待所有方法完成
   - 监控进程退出状态
   - 记录成功/失败状态
   ↓
6. 清理临时文件
   - 删除临时符号链接（如deepbsa_input.vcf）
   ↓
7. 自动合并结果 ⭐
   - 合并CSV结果（添加Method列）
   - 复制PNG图片（添加方法名前缀）
   - 生成汇总报告
   ↓
8. 输出最终摘要
   - 显示各方法成功/失败状态
   - 显示输出目录结构
```

⭐ = 本模块新增或自动处理的功能

## 关键技术实现 | Key Technical Implementation

### 1. 智能文件名识别

**问题**: DeepBSA使用简单的`split('.')`解析文件名，无法处理`variation.filtered.snp.biallelic.vcf`

**解决方案**:
```python
# 使用pathlib.Path.suffix智能识别扩展名
file_ext = input_file.suffix.lower()  # .vcf

if file_ext in ['.vcf', '.vcf.gz', '.vcf.bz2']:
    if '.' in input_file.stem:  # 文件名包含多个点
        # 创建符号链接: deepbsa_input.vcf → 实际文件
        temp_link_name = self.config.output_path / f"deepbsa_input{file_ext}"
        temp_link_name.symlink_to(input_file.absolute())
```

### 2. 方法级并行执行

**特点**:
- 每个方法在独立工作目录运行
- 所有方法同时启动（并行模式）
- 独立日志文件，便于调试
- 避免文件冲突

**实现**:
```python
for method in methods_list:
    method_work_dir = self.config.output_path / method
    method_work_dir.mkdir(exist_ok=True)

    cmd = [
        "conda", "run", "-n", conda_env_name, "--no-capture-output",
        "python3", deepbsa_script,
        "--i", input_file,
        "--m", method,
        "--p", "0"  # 禁用pretreatment
    ]

    subprocess.Popen(cmd, cwd=method_work_dir, stdout=log_file)
```

### 3. DL方法模型文件修复

**问题**: DL方法需要预训练模型文件（.h5），但找不到`Models/row_finetune2pool.h5`

**解决方案**:
```python
# 为每个方法创建Models符号链接
models_src = self.config.deepbsa_script.parent / "Models"
models_link = method_work_dir / "Models"
if models_src.exists() and not models_link.exists():
    models_link.symlink_to(models_src)
```

### 4. 自动结果合并

**CSV合并**:
```python
for method in methods:
    csv_files = glob(f"{method}/Results/variant/*.csv")
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        df_filtered = df[df['Value'] != '-']  # 只保留有效结果
        df_filtered['Method'] = method         # 添加方法名列
        all_results.append(df_filtered)

merged_df = pd.concat(all_results)
merged_df.to_csv("merged_results/merged_results.csv")
```

**PNG复制**:
```python
for method in methods:
    png_files = glob(f"{method}/Results/variant/*.png")
    for png_file in png_files:
        new_name = f"{method}_{png_file.name}"  # 添加方法名前缀
        shutil.copy(png_file, f"merged_results/images/{new_name}")
```

## 使用示例 | Usage Examples

### 示例1: 标准BSA分析

```bash
# 运行所有7种方法，自动合并结果
biopytools deepbsa \
    -i /path/to/variant.vcf \
    -o /path/to/bsa_output

# 查看合并结果
cat /path/to/bsa_output/merged_results/merged_results.csv
cat /path/to/bsa_output/merged_results/summary_report.txt
ls /path/to/bsa_output/merged_results/images/
```

### 示例2: 只运行特定方法

```bash
# 只运行ED4和SNP-Index方法
biopytools deepbsa \
    -i variant.vcf \
    -m ED4,SNP \
    -o bsab_output
```

### 示例3: 调试模式（串行运行）

```bash
# 串行运行，便于查看每个方法的输出
biopytools deepbsa \
    -i variant.vcf \
    --no-parallel \
    -o bsab_output

# 查看每个方法的独立日志
tail -f bsab_output/DL/DL.log
tail -f bsab_output/ED4/ED4.log
```

### 示例4: 自定义DeepBSA路径

```bash
# 使用自定义的DeepBSA安装路径
biopytools deepbsa \
    -i variant.vcf \
    --deepbsa-path /custom/path/DeepBSA/bin/main.py \
    --conda-env /custom/path/conda/envs/deepbsa \
    -o bsab_output
```

## 故障排查 | Troubleshooting

### 1. DL方法失败

**症状**: `OSError: No file or directory found at .../Models/row_finetune2pool.h5`

**解决**:
- 本模块已自动修复，会创建Models符号链接
- 检查DeepBSA安装目录是否存在`Models/`文件夹

### 2. 复杂文件名无法识别

**症状**: `File type not recognized` 或 `ValueError: not enough values to unpack`

**解决**:
- 本模块已自动修复，会创建临时符号链接
- 确保输入文件扩展名为`.vcf`、`.vcf.gz`或`.vcf.bz2`

### 3. Pandas兼容性错误

**症状**: `TypeError: to_csv() got an unexpected keyword argument 'lineterminator'`

**解决**:
- 本模块已修改DeepBSA源码（`vcf_handle.py:131`）
- 移除了`lineterminator`参数以兼容旧版pandas

### 4. 某个方法失败

**症状**: 日志显示某个方法failed

**解决**:
```bash
# 查看该方法的独立日志
cat output_dir/METHOD_NAME/METHOD_NAME.log

# 常见原因：
# - 输入VCF文件格式不正确
# - VCF文件没有足够的SNP位点
# - Conda环境配置问题
```

### 5. 并行执行内存不足

**症状**: 系统内存耗尽

**解决**:
```bash
# 使用串行模式
biopytools deepbsa -i variant.vcf --no-parallel -o output
```

## 性能优化 | Performance Optimization

### 并行执行效率

- **7个方法并行**: 约7-8小时（基于小花糖芥数据集）
- **串行执行**: 约20-24小时
- **加速比**: 约3倍

### 磁盘空间需求

每个方法约产生20-25MB输出，7个方法总计约140-175MB

### 内存需求

- **并行模式**: 建议至少16GB内存
- **串行模式**: 建议8GB内存

## 与原始DeepBSA的对比 | Comparison with Original DeepBSA

| 特性 | 原始DeepBSA | biopytools deepbsa |
|------|-------------|-------------------|
| 文件名识别 | 仅支持简单文件名 | ✅ 支持复杂文件名 |
| 并行执行 | ❌ 需要手动启动 | ✅ 自动并行 |
| 工作目录 | 混在一起 | ✅ 独立目录 |
| DL方法 | ❌ 手动创建Models链接 | ✅ 自动创建 |
| 结果合并 | ❌ 手动整合 | ✅ 自动合并 |
| 日志管理 | 分散 | ✅ 统一管理 |
| Pandas兼容性 | ❌ 报错 | ✅ 已修复 |

## 参考文献 | References

1. **DeepBSA**: 原始DeepBSA工具文档
2. **BSA方法**: 参见各方法的原始文献

## 版本历史 | Version History

- **v1.0** (2026-03-30): 初始版本
  - 支持7种BSA方法批量运行
  - 智能文件名识别
  - 方法级并行执行
  - 自动结果合并
  - DL方法修复
  - Pandas兼容性修复

## 贡献者 | Contributors

- biopytools开发团队

## 许可证 | License

本模块遵循biopytools项目的许可证。DeepBSA原始工具请参考其各自的许可证。
