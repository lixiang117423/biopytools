# VcfFilter模块集成指南

## 手动集成步骤

### 1. 更新主包 `biopytools/__init__.py`

添加以下内容:

```python
# 添加导入
try:
    from .vcf_filter import VcfFilterProcessor, VcfFilterConfig
except ImportError:
    VcfFilterProcessor = None
    VcfFilterConfig = None

# 更新__all__列表
__all__ = [
    # ... 现有的
    "VcfFilterProcessor",
    "VcfFilterConfig",
]

# 更新list_tools函数
def list_tools():
    tools = []
    # ... 现有的工具
    if VcfFilterProcessor:
        tools.append(("vcf_filter", "VcfFilter工具"))
    return tools
```

### 2. 更新CLI主入口 `biopytools/cli/main.py`

添加以下内容:

```python
# 导入命令
from .commands.vcf_filter import vcf_filter

# 添加命令到main group
main.add_command(vcf_filter)
```

### 3. 更新 `pyproject.toml`

在 `[project.scripts]` 部分添加:

```toml
biopytools-vcf_filter = "biopytools.vcf_filter.main:main"
```

### 4. 验证集成

运行以下命令验证:

```bash
# 测试导入
python -c "from biopytools.vcf_filter import VcfFilterProcessor; print('✅ 导入成功')"

# 测试CLI
python scripts/run_vcf_filter.py --help

# 运行测试
pytest tests/test_vcf_filter.py -v
```

## 文件位置总结

- 主模块: `biopytools/vcf_filter/`
- CLI命令: `biopytools/cli/commands/vcf_filter.py`
- 测试: `tests/test_vcf_filter.py`
- 脚本: `scripts/run_vcf_filter.py`
- 示例: `examples/vcf_filter/`
