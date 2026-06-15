# HapHiC模块更新日志 | HapHiC Module Changelog

## v0.13.1 (2024-12-20) - Pipeline模式与断点续传

### 🚀 重大更新 | Major Updates

#### Pipeline模式集成 | Pipeline Mode Integration
- ✅ **原生Pipeline支持**: 采用HapHiC v1.0.7原生Pipeline模式
- ⚡ **一步执行**: 单命令完成所有scaffolding步骤
- 🛡️ **更好兼容性**: 完全兼容HapHiC原生参数和功能

#### 断点续传功能 | Resume Feature
- 🔄 **智能进度检测**: 自动检测各步骤完成状态
- 📍 **中断恢复**: 支持从任意中断点继续执行
- ⚡ **时间节省**: 避免重复计算，提高效率
- 🛠️ **强制重跑**: 提供`--force-rerun`选项完全重新执行

#### 自动可视化 | Automatic Visualization
- 📊 **默认生成**: Pipeline完成后自动生成Hi-C接触图
- 🎨 **多格式支持**: 输出PDF/PNG格式可视化文件
- 📁 **独立目录**: 存储在`05.plots/`目录
- ⚙️ **可配置**: 支持bin_size、min_len等参数调整

#### Juicebox集成优化 | Juicebox Integration
- 🥤 **自动生成**: 默认生成Juicebox兼容文件
- 📂 **标准格式**: 输出`.hic`和`.assembly`文件
- 🔧 **工具集成**: 自动调用matlock、agp2assembly、asm-visualizer
- 📁 **独立目录**: 存储在`06.juicebox/`目录

### 🔧 技术改进 | Technical Improvements

#### 文件路径修正 | File Path Corrections
- ✅ **AGP文件路径**: 修正为实际生成的位置`04.build/{prefix}.agp`
- 📁 **输出目录**: 更新文档反映真实的文件结构
- 🔍 **错误修复**: 修复Juicebox文件生成路径问题

#### 目录管理优化 | Directory Management
- 🗂️ **智能清理**: 只清理必要目录，保留已完成步骤
- 📁 **结构清晰**: 明确各步骤目录用途（01-06）
- 🔒 **数据安全**: 避免意外删除已完成的结果

#### 性能优化 | Performance Optimization
- ⚡ **减少I/O**: Pipeline模式减少文件读写开销
- 🔄 **并行处理**: 支持多进程并行执行
- 💾 **内存优化**: 智能内存管理，支持大基因组

### 📋 新增参数 | New Parameters

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--force-rerun` | flag | False | 强制重新运行所有步骤 |
| `--no-generate-plots` | flag | - | 禁用可视化图表生成 |
| `--no-generate-juicebox` | flag | - | 禁用Juicebox文件生成 |

### 🔄 行为变更 | Behavior Changes

#### 默认行为调整 | Default Behavior Adjustments
- 📊 **可视化**: 默认启用可视化生成（之前是可选）
- 🥤 **Juicebox**: 默认启用Juicebox文件生成
- 🔄 **断点续传**: 默认启用断点续传模式

#### 输出文件变更 | Output File Changes
```
新结构 | New Structure:
├── 04.build/{prefix}.fa          # 最终scaffolds
├── 04.build/{prefix}.agp         # AGP文件
├── 05.plots/*.pdf                # 可视化图表
├── 06.juicebox/{prefix}.hic      # Juicebox文件
└── 06.juicebox/{prefix}.assembly # Assembly文件
```

### 🐛 问题修复 | Bug Fixes
- ✅ 修复AGP文件路径错误
- ✅ 修复可视化生成失败问题
- ✅ 修复Juicebox文件生成路径
- ✅ 改进错误处理和日志记录

### 📚 文档更新 | Documentation Updates
- 📖 **全新文档**: 完全重写HapHiC模块文档
- 📝 **参数详解**: 详细说明所有参数选项
- 🎯 **使用示例**: 添加丰富的使用示例
- 🔧 **故障排除**: 常见问题和解决方案

### 🚨 重要提示 | Important Notes

#### 迁移指南 | Migration Guide
- ✅ **命令兼容**: 原有命令保持兼容
- 🔄 **行为变化**: 默认启用可视化和Juicebox生成
- 📁 **输出变化**: 文件存储位置有所调整

#### 建议配置 | Recommended Configuration
```bash
# 标准配置（推荐）
biopytools haphic -a assembly.fa -b hic.bam -c 24

# 高性能配置
biopytools haphic -a assembly.fa -b hic.bam -c 24 \
    --threads 32 --processes 16 --correct-nrounds 2

# 如需禁用新功能
biopytools haphic -a assembly.fa -b hic.bam -c 24 \
    --no-generate-plots --no-generate-juicebox
```

## 早期版本 | Earlier Versions

### v0.13.0
- 🧬 初始HapHiC模块集成
- ⚡ 支持基本的Hi-C scaffolding功能
- 📊 支持BWA比对和HapHiC过滤

---

**更新日期 | Update Date**: 2024-12-20
**维护者 | Maintainer**: BioPyTools Team