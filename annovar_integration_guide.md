# Annovar模块集成指南

## 手动集成步骤

### 1. 更新主包 `biopytools/__init__.py`

添加以下内容:

```python
# 添加导入
try:
    from .annovar import AnnovarProcessor, AnnovarConfig
except ImportError:
    AnnovarProcessor = None
    AnnovarConfig = None

# 更新__all__列表
__all__ = [
    # ... 现有的
    "AnnovarProcessor",
    "AnnovarConfig",
]

# 更新list_tools函数
def list_tools():
    tools = []
    # ... 现有的工具
    if AnnovarProcessor:
        tools.append(("annovar", "Annovar工具"))
    return tools
```

### 2. 更新CLI主入口 `biopytools/cli/main.py`

添加以下内容:

```python
# 导入命令
from .commands.annovar import annovar

# 添加命令到main group
main.add_command(annovar)
```

### 3. 更新 `pyproject.toml`

在 `[project.scripts]` 部分添加:

```toml
biopytools-annovar = "biopytools.annovar.main:main"
```

### 4. 验证集成

运行以下命令验证:

```bash
# 测试导入
python -c "from biopytools.annovar import AnnovarProcessor; print('✅ 导入成功')"

# 测试CLI
python scripts/run_annovar.py --help

# 运行测试
pytest tests/test_annovar.py -v
```

## 文件位置总结

- 主模块: `biopytools/annovar/`
- CLI命令: `biopytools/cli/commands/annovar.py`
- 测试: `tests/test_annovar.py`
- 脚本: `scripts/run_annovar.py`
- 示例: `examples/annovar/`
