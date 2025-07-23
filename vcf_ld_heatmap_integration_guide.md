# VcfLdHeatmap模块集成指南

## 手动集成步骤

### 1. 更新主包 `biopytools/__init__.py`

添加以下内容:

```python
# 添加导入
try:
    from .vcf_ld_heatmap import VcfLdHeatmapProcessor, VcfLdHeatmapConfig
except ImportError:
    VcfLdHeatmapProcessor = None
    VcfLdHeatmapConfig = None

# 更新__all__列表
__all__ = [
    # ... 现有的
    "VcfLdHeatmapProcessor",
    "VcfLdHeatmapConfig",
]

# 更新list_tools函数
def list_tools():
    tools = []
    # ... 现有的工具
    if VcfLdHeatmapProcessor:
        tools.append(("vcf_ld_heatmap", "VcfLdHeatmap工具"))
    return tools
```

### 2. 更新CLI主入口 `biopytools/cli/main.py`

添加以下内容:

```python
# 导入命令
from .commands.vcf_ld_heatmap import vcf_ld_heatmap

# 添加命令到main group
main.add_command(vcf_ld_heatmap)
```

### 3. 更新 `pyproject.toml`

在 `[project.scripts]` 部分添加:

```toml
biopytools-vcf_ld_heatmap = "biopytools.vcf_ld_heatmap.main:main"
```

### 4. 验证集成

运行以下命令验证:

```bash
# 测试导入
python -c "from biopytools.vcf_ld_heatmap import VcfLdHeatmapProcessor; print('✅ 导入成功')"

# 测试CLI
python scripts/run_vcf_ld_heatmap.py --help

# 运行测试
pytest tests/test_vcf_ld_heatmap.py -v
```

## 文件位置总结

- 主模块: `biopytools/vcf_ld_heatmap/`
- CLI命令: `biopytools/cli/commands/vcf_ld_heatmap.py`
- 测试: `tests/test_vcf_ld_heatmap.py`
- 脚本: `scripts/run_vcf_ld_heatmap.py`
- 示例: `examples/vcf_ld_heatmap/`
