# 路径管理详解|Path Management Details

> 📖 本文件为 CLAUDE.md 的按需参考文档,**改路径管理 / 迁移旧代码时读取**。核心规则(禁硬编码 / 优先级 / 必须展开~)见 CLAUDE.md §十一。
> On-demand reference for CLAUDE.md; read when touching path management or migrating legacy code.

---

## 1. `common/paths.py` 完整实现|Full Implementation

> 以下为 `common/paths.py` 的参考实现，所有模块直接从此导入|Reference implementation for all modules to import from:

```python
# common/paths.py
"""路径管理工具|Path management utilities"""

import os
import shutil
import logging
from typing import Optional

logger = logging.getLogger(__name__)


def expand_path(path: str) -> str:
    """
    展开路径中的~和环境变量|Expand ~ and environment variables in path

    Args:
        path: 原始路径（可包含~或$VAR）|Raw path (may contain ~ or $VAR)

    Returns:
        展开后的绝对路径|Expanded absolute path

    Examples:
        >>> expand_path("~/miniforge3/bin/tool")
        "/home/user/miniforge3/bin/tool"
        >>> expand_path("$SOFTWARE/bin/tool")
        "/opt/software/bin/tool"
    """
    return os.path.expandvars(os.path.expanduser(path))


def get_tool_path(
    tool_name: str,
    default_path: str,
    env_var: Optional[str] = None
) -> str:
    """
    按优先级获取工具路径|Get tool path by priority

    优先级（高→低）|Priority (high→low):
        1. 环境变量 (env_var)|Environment variable
        2. 用户配置文件 (~/.config/biopytools/config.yml)
        3. 代码默认值 (default_path)|Code default

    Args:
        tool_name: 工具名称（用于配置文件查找）|Tool name (for config file lookup)
        default_path: 默认路径（支持~展开）|Default path (supports ~ expansion)
        env_var: 环境变量名（可选）|Environment variable name (optional)

    Returns:
        展开后的工具路径|Expanded tool path
    """
    # 1. 环境变量优先|Env var takes priority
    if env_var and os.environ.get(env_var):
        path = expand_path(os.environ[env_var])
        logger.debug(f"使用环境变量路径|Using env var path [{env_var}]: {path}")
        return path

    # 2. 用户配置文件|User config file
    config_path = os.path.expanduser("~/.config/biopytools/config.yml")
    if os.path.exists(config_path):
        try:
            import yaml
            with open(config_path) as f:
                config = yaml.safe_load(f) or {}
            tool_path = config.get("tools", {}).get(tool_name)
            if tool_path:
                path = expand_path(tool_path)
                logger.debug(f"使用配置文件路径|Using config file path [{tool_name}]: {path}")
                return path
        except Exception as e:
            logger.warning(f"读取配置文件失败|Failed to read config file: {e}")

    # 3. 代码默认值|Code default
    path = expand_path(default_path)
    logger.debug(f"使用默认路径|Using default path [{tool_name}]: {path}")
    return path
```

### config.py 中使用|Use in config.py

```python
from dataclasses import dataclass, field
from common.paths import get_tool_path

@dataclass
class ModuleConfig:
    tool_path: str = field(
        default_factory=lambda: get_tool_path(
            'tool_name',
            '~/default/path/to/tool',
            'TOOL_PATH'  # 可选环境变量|Optional env var
        )
    )
```

### main.py 中使用|Use in main.py

```python
from common.paths import expand_path

# 展开用户提供的路径|Expand user-provided path
config.tool_path = expand_path(args.tool_path)
```

---

## 2. ⚠️ 必须展开~路径详解|Must Expand ~ Paths (Detail)

**问题|Problem:** Python不会自动展开字符串中的`~`符号|Python does NOT auto-expand `~` in strings:

```python
# ❌ 错误示例|Wrong Example
meme_path: str = "~/miniforge3/envs/meme_v.5.5.9/bin/meme"

# 在验证时会失败|Will fail in validation:
import shutil
shutil.which(meme_path)  # 返回None|Returns None ❌
os.path.exists(meme_path)  # 返回False|Returns False ❌
```

**解决方案|Solution:** 在`config.py`的`__post_init__`方法中必须展开所有包含`~`的路径:

```python
from dataclasses import dataclass
from ..common.paths import expand_path

@dataclass
class ModuleConfig:
    meme_path: str = '~/miniforge3/envs/meme_v.5.5.9/bin/meme'
    tool_path: str = '~/software/bin/tool'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # ⚠️ 关键：展开所有包含~的路径|CRITICAL: Expand all paths with ~
        self.meme_path = expand_path(self.meme_path)
        self.tool_path = expand_path(self.tool_path)

        # 现在可以正常检测工具|Now tool detection works:
        # shutil.which(self.meme_path)  # 返回完整路径|Returns full path ✓
```

**常见错误|Common Mistakes:**
1. ❌ 只在定义时使用`expand_path()`，不在`__post_init__`中展开
2. ❌ 忘记展开用户通过命令行参数传入的路径
3. ❌ 使用`os.path.expanduser()`而不是统一的`expand_path()`

---

## 3. 迁移现有代码|Migrating Existing Code

使用迁移脚本|Use migration script:

```bash
cd ~/software/biopytools
python scripts/migrate_paths.py --module module_name --dry-run  # 先试运行|Dry run first
python scripts/migrate_paths.py --module module_name --apply   # 确认后应用|Apply after review
```

验证迁移结果|Verify migration result:

```bash
grep -r "/share/org/YZWL/yzwl_lixg/" biopytools/module_name/
# 应该返回0个结果|Should return 0 results
```
