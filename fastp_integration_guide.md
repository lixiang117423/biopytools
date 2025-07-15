# Fastp模块集成指南

## 手动集成步骤

### 1. 更新主包 `biopytools/__init__.py`

添加以下内容:

```python
# 添加导入
try:
    from .fastp import FastpProcessor, FastpConfig
except ImportError:
    FastpProcessor = None
    FastpConfig = None

# 更新__all__列表
__all__ = [
    # ... 现有的
    "FastpProcessor",
    "FastpConfig",
]

# 更新list_tools函数
def list_tools():
    tools = []
    # ... 现有的工具
    if FastpProcessor:
        tools.append(("fastp", "Fastp工具"))
    return tools
```

### 2. 更新CLI主入口 `biopytools/cli/main.py`

添加以下内容:

```python
# 导入命令
from .commands.fastp import fastp

# 添加命令到main group
main.add_command(fastp)
```

### 3. 更新 `pyproject.toml`

在 `[project.scripts]` 部分添加:

```toml
biopytools-fastp = "biopytools.fastp.main:main"
```

### 4. 验证集成

运行以下命令验证:

```bash
# 测试导入
python -c "from biopytools.fastp import FastpProcessor; print('✅ 导入成功')"

# 测试CLI
python scripts/run_fastp.py --help

# 运行测试
pytest tests/test_fastp.py -v
```

## 文件位置总结

- 主模块: `biopytools/fastp/`
- CLI命令: `biopytools/cli/commands/fastp.py`
- 测试: `tests/test_fastp.py`
- 脚本: `scripts/run_fastp.py`
- 示例: `examples/fastp/`
