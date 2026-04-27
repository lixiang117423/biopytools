"""
版本信息模块|Version Information Module
"""

import sys
import os
import subprocess
import re
from datetime import datetime
from typing import Dict, Tuple, Any

# 从 pyproject.toml 动态读取版本号
def _get_version_from_pyproject() -> str:
    """从 pyproject.toml 读取版本号"""
    try:
        # 获取 pyproject.toml 路径（项目根目录）
        project_root = os.path.dirname(os.path.dirname(__file__))
        pyproject_path = os.path.join(project_root, 'pyproject.toml')

        if os.path.exists(pyproject_path):
            with open(pyproject_path, 'r', encoding='utf-8') as f:
                content = f.read()
                # 使用正则表达式提取 version = "x.x.x"
                match = re.search(r'^version\s*=\s*["\']([^"\']+)["\']', content, re.MULTILINE)
                if match:
                    return match.group(1)
    except Exception:
        pass

    # fallback 版本号
    return "0.0.0"

# 版本信息（从 pyproject.toml 读取）
_version_str = _get_version_from_pyproject()
__version__ = _version_str

# 解析版本号元组
def _parse_version(version_str: str) -> Tuple[int, ...]:
    """解析版本字符串为元组"""
    parts = version_str.split('.')
    return tuple(int(p) for p in parts[:3]) if len(parts) >= 3 else (0, 0, 0)

__version_info__ = _parse_version(__version__)

# 版本状态
VERSION_STATUS = "stable"  # alpha, beta, rc, stable

def get_git_revision() -> str:
    """
    获取Git提交信息
    Get Git commit information
    """
    try:
        # 获取最新commit的短hash
        commit = subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD'],
            stderr=subprocess.DEVNULL,
            universal_newlines=True,
            cwd=os.path.dirname(os.path.dirname(__file__))
        ).strip()
        
        # 检查是否有未提交的更改
        dirty = subprocess.call(
            ['git', 'diff-index', '--quiet', 'HEAD'],
            stderr=subprocess.DEVNULL,
            cwd=os.path.dirname(os.path.dirname(__file__))
        ) != 0
        
        return f"{commit}{'-dirty' if dirty else ''}"
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return "unknown"

def get_build_info() -> Dict[str, Any]:
    """
    获取构建信息
    Get build information
    """
    return {
        'version': __version__,
        'version_info': __version_info__,
        'status': VERSION_STATUS,
        'git_revision': get_git_revision(),
        'build_date': datetime.now().isoformat(),
        'python_version': f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        'python_implementation': sys.implementation.name,
        'platform': sys.platform
    }

def get_version_string() -> str:
    """
    获取完整版本字符串
    Get complete version string
    """
    git_rev = get_git_revision()
    if VERSION_STATUS != "stable":
        return f"{__version__}-{VERSION_STATUS}+{git_rev}"
    elif git_rev != "unknown":
        return f"{__version__}+{git_rev}"
    else:
        return __version__

def compare_versions(version1: str, version2: str) -> int:
    """
    比较两个版本号
    Compare two version numbers
    
    Returns:
        -1 if version1 < version2
         0 if version1 == version2
         1 if version1 > version2
    """
    def parse_version(version: str) -> Tuple[int, ...]:
        # 移除预发布标识符
        version = version.split('-')[0].split('+')[0]
        return tuple(map(int, version.split('.')))
    
    v1_parts = parse_version(version1)
    v2_parts = parse_version(version2)
    
    # 补齐版本号位数
    max_len = max(len(v1_parts), len(v2_parts))
    v1_parts = v1_parts + (0,) * (max_len - len(v1_parts))
    v2_parts = v2_parts + (0,) * (max_len - len(v2_parts))
    
    if v1_parts < v2_parts:
        return -1
    elif v1_parts > v2_parts:
        return 1
    else:
        return 0

def is_compatible_version(required_version: str, current_version: str = None) -> bool:
    """
    检查版本兼容性
    Check version compatibility
    
    Args:
        required_version: 要求的最低版本
        current_version: 当前版本（默认使用当前包版本）
    
    Returns:
        bool: 是否兼容
    """
    if current_version is None:
        current_version = __version__
    
    return compare_versions(current_version, required_version) >= 0

# 导出的公共接口
__all__ = [
    '__version__',
    '__version_info__',
    'VERSION_STATUS',
    'get_version_string',
    'get_build_info',
    'get_git_revision',
    'compare_versions',
    'is_compatible_version'
]