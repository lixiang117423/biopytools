# MethylationPipeline模块集成指南

## 手动集成步骤

### 1. 更新主包 `biopytools/__init__.py`

添加以下内容:

```python
# 添加导入
try:
    from .methylation_pipeline import MethylationPipelineProcessor, MethylationPipelineConfig
except ImportError:
    MethylationPipelineProcessor = None
    MethylationPipelineConfig = None

# 更新__all__列表
__all__ = [
    # ... 现有的
    "MethylationPipelineProcessor",
    "MethylationPipelineConfig",
]

# 更新list_tools函数
def list_tools():
    tools = []
    # ... 现有的工具
    if MethylationPipelineProcessor:
        tools.append(("methylation_pipeline", "MethylationPipeline工具"))
    return tools
```

### 2. 更新CLI主入口 `biopytools/cli/main.py`

添加以下内容:

```python
# 导入命令
from .commands.methylation_pipeline import methylation_pipeline

# 添加命令到main group
main.add_command(methylation_pipeline)
```

### 3. 更新 `pyproject.toml`

在 `[project.scripts]` 部分添加:

```toml
biopytools-methylation_pipeline = "biopytools.methylation_pipeline.main:main"
```

### 4. 验证集成

运行以下命令验证:

```bash
# 测试导入
python -c "from biopytools.methylation_pipeline import MethylationPipelineProcessor; print('✅ 导入成功')"

# 测试CLI
python scripts/run_methylation_pipeline.py --help

# 运行测试
pytest tests/test_methylation_pipeline.py -v
```

## 文件位置总结

- 主模块: `biopytools/methylation_pipeline/`
- CLI命令: `biopytools/cli/commands/methylation_pipeline.py`
- 测试: `tests/test_methylation_pipeline.py`
- 脚本: `scripts/run_methylation_pipeline.py`
- 示例: `examples/methylation_pipeline/`
