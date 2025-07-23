# EnaDownloader模块集成指南

## 手动集成步骤

### 1. 更新主包 `biopytools/__init__.py`

添加以下内容:

```python
# 添加导入
try:
    from .ena_downloader import EnaDownloaderProcessor, EnaDownloaderConfig
except ImportError:
    EnaDownloaderProcessor = None
    EnaDownloaderConfig = None

# 更新__all__列表
__all__ = [
    # ... 现有的
    "EnaDownloaderProcessor",
    "EnaDownloaderConfig",
]

# 更新list_tools函数
def list_tools():
    tools = []
    # ... 现有的工具
    if EnaDownloaderProcessor:
        tools.append(("ena_downloader", "EnaDownloader工具"))
    return tools
```

### 2. 更新CLI主入口 `biopytools/cli/main.py`

添加以下内容:

```python
# 导入命令
from .commands.ena_downloader import ena_downloader

# 添加命令到main group
main.add_command(ena_downloader)
```

### 3. 更新 `pyproject.toml`

在 `[project.scripts]` 部分添加:

```toml
biopytools-ena_downloader = "biopytools.ena_downloader.main:main"
```

### 4. 验证集成

运行以下命令验证:

```bash
# 测试导入
python -c "from biopytools.ena_downloader import EnaDownloaderProcessor; print('✅ 导入成功')"

# 测试CLI
python scripts/run_ena_downloader.py --help

# 运行测试
pytest tests/test_ena_downloader.py -v
```

## 文件位置总结

- 主模块: `biopytools/ena_downloader/`
- CLI命令: `biopytools/cli/commands/ena_downloader.py`
- 测试: `tests/test_ena_downloader.py`
- 脚本: `scripts/run_ena_downloader.py`
- 示例: `examples/ena_downloader/`
