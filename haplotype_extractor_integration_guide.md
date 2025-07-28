# HaplotypeExtractor模块集成指南

## 手动集成步骤

### 1. 更新主包 `biopytools/__init__.py`

添加以下内容:

```python
# 添加导入
try:
    from .haplotype_extractor import HaplotypeExtractorProcessor, HaplotypeExtractorConfig
except ImportError:
    HaplotypeExtractorProcessor = None
    HaplotypeExtractorConfig = None

# 更新__all__列表
__all__ = [
    # ... 现有的
    "HaplotypeExtractorProcessor",
    "HaplotypeExtractorConfig",
]

# 更新list_tools函数
def list_tools():
    tools = []
    # ... 现有的工具
    if HaplotypeExtractorProcessor:
        tools.append(("haplotype_extractor", "HaplotypeExtractor工具"))
    return tools
```

### 2. 更新CLI主入口 `biopytools/cli/main.py`

添加以下内容:

```python
# 导入命令
from .commands.haplotype_extractor import haplotype_extractor

# 添加命令到main group
main.add_command(haplotype_extractor)
```

### 3. 更新 `pyproject.toml`

在 `[project.scripts]` 部分添加:

```toml
biopytools-haplotype_extractor = "biopytools.haplotype_extractor.main:main"
```

### 4. 验证集成

运行以下命令验证:

```bash
# 测试导入
python -c "from biopytools.haplotype_extractor import HaplotypeExtractorProcessor; print('✅ 导入成功')"

# 测试CLI
python scripts/run_haplotype_extractor.py --help

# 运行测试
pytest tests/test_haplotype_extractor.py -v
```

## 文件位置总结

- 主模块: `biopytools/haplotype_extractor/`
- CLI命令: `biopytools/cli/commands/haplotype_extractor.py`
- 测试: `tests/test_haplotype_extractor.py`
- 脚本: `scripts/run_haplotype_extractor.py`
- 示例: `examples/haplotype_extractor/`
