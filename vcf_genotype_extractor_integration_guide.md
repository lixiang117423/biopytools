# VcfGenotypeExtractor模块集成指南

## 手动集成步骤

### 1. 更新主包 `biopytools/__init__.py`

添加以下内容:

```python
# 添加导入
try:
    from .vcf_genotype_extractor import VcfGenotypeExtractorProcessor, VcfGenotypeExtractorConfig
except ImportError:
    VcfGenotypeExtractorProcessor = None
    VcfGenotypeExtractorConfig = None

# 更新__all__列表
__all__ = [
    # ... 现有的
    "VcfGenotypeExtractorProcessor",
    "VcfGenotypeExtractorConfig",
]

# 更新list_tools函数
def list_tools():
    tools = []
    # ... 现有的工具
    if VcfGenotypeExtractorProcessor:
        tools.append(("vcf_genotype_extractor", "VcfGenotypeExtractor工具"))
    return tools
```

### 2. 更新CLI主入口 `biopytools/cli/main.py`

添加以下内容:

```python
# 导入命令
from .commands.vcf_genotype_extractor import vcf_genotype_extractor

# 添加命令到main group
main.add_command(vcf_genotype_extractor)
```

### 3. 更新 `pyproject.toml`

在 `[project.scripts]` 部分添加:

```toml
biopytools-vcf_genotype_extractor = "biopytools.vcf_genotype_extractor.main:main"
```

### 4. 验证集成

运行以下命令验证:

```bash
# 测试导入
python -c "from biopytools.vcf_genotype_extractor import VcfGenotypeExtractorProcessor; print('✅ 导入成功')"

# 测试CLI
python scripts/run_vcf_genotype_extractor.py --help

# 运行测试
pytest tests/test_vcf_genotype_extractor.py -v
```

## 文件位置总结

- 主模块: `biopytools/vcf_genotype_extractor/`
- CLI命令: `biopytools/cli/commands/vcf_genotype_extractor.py`
- 测试: `tests/test_vcf_genotype_extractor.py`
- 脚本: `scripts/run_vcf_genotype_extractor.py`
- 示例: `examples/vcf_genotype_extractor/`
