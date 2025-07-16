# Biohelpers模块集成指南

## 手动集成步骤

### 1. 更新主包 `biopytools/__init__.py`

添加以下内容:

```python
# 添加导入
try:
    from .biohelpers import BiohelpersProcessor, BiohelpersConfig
except ImportError:
    BiohelpersProcessor = None
    BiohelpersConfig = None

# 更新__all__列表
__all__ = [
    # ... 现有的
    "BiohelpersProcessor",
    "BiohelpersConfig",
]

# 更新list_tools函数
def list_tools():
    tools = []
    # ... 现有的工具
    if BiohelpersProcessor:
        tools.append(("biohelpers", "Biohelpers工具"))
    return tools
```

### 2. 更新CLI主入口 `biopytools/cli/main.py`

添加以下内容:

```python
# 导入命令
from .commands.biohelpers import biohelpers

# 添加命令到main group
main.add_command(biohelpers)
```

### 3. 更新 `pyproject.toml`

在 `[project.scripts]` 部分添加:

```toml
biopytools-biohelpers = "biopytools.biohelpers.main:main"
```

### 4. 验证集成

运行以下命令验证:

```bash
# 测试导入
python -c "from biopytools.biohelpers import BiohelpersProcessor; print('✅ 导入成功')"

# 测试CLI
python scripts/run_biohelpers.py --help

# 运行测试
pytest tests/test_biohelpers.py -v
```

## 文件位置总结

- 主模块: `biopytools/biohelpers/`
- CLI命令: `biopytools/cli/commands/biohelpers.py`
- 测试: `tests/test_biohelpers.py`
- 脚本: `scripts/run_biohelpers.py`
- 示例: `examples/biohelpers/`
