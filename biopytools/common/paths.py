"""路径管理工具|Path management utilities

提供统一的工具路径获取方式，支持环境变量、配置文件和默认路径
| Provide unified tool path resolution with env vars, config files, and defaults

优先级|Priority:
1. 环境变量|Environment variable
2. 用户配置文件|User config file (~/.config/biopytools/config.yml)
3. 默认路径（展开~）|Default path (expand ~)
"""

import os
import re
import logging
import yaml
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# 用户配置文件路径|User config file path
USER_CONFIG_FILE = Path("~/.config/biopytools/config.yml").expanduser()


def load_user_config() -> dict:
    """加载用户配置文件|Load user configuration file

    Returns:
        dict: 配置字典|Configuration dict, 空字典如果文件不存在|empty dict if file not found
    """
    if USER_CONFIG_FILE.exists():
        try:
            with open(USER_CONFIG_FILE, 'r') as f:
                return yaml.safe_load(f) or {}
        except Exception:
            # 配置文件解析失败，返回空字典|Config parse failed, return empty dict
            return {}
    return {}


def expand_path(path: str) -> str:
    """展开路径中的~和环境变量|Expand ~ and environment variables in path

    Args:
        path: 路径字符串|Path string, 支持~和$VAR|Supports ~ and $VAR

    Returns:
        str: 展开后的绝对路径|Expanded absolute path

    Examples:
        >>> expand_path("~/miniforge3/bin/fanc")
        '/home/user/miniforge3/bin/fanc'
        >>> expand_path("$SOFTWARE/fanc/bin/fanc")
        '/path/to/software/fanc/bin/fanc'
    """
    expanded = os.path.expandvars(path)
    expanded = os.path.expanduser(expanded)
    if '/' not in expanded and '.' not in expanded:
        return expanded
    return os.path.abspath(expanded)


def resolve_legacy_path(
    parent_dir,
    canonical_name: str,
    log=None,
    create: bool = False,
):
    """定位步骤目录/文件，优先下划线规范名，回退点号老名用于断点续传
    |Resolve a step dir/file: prefer the canonical underscore name, fall back
    to the legacy dot name (NN_xxx -> NN.xxx) so existing on-disk output from
    older runs is reused instead of being recomputed.

    Args:
        parent_dir: 父目录|parent directory
        canonical_name: 规范名(下划线)，如 "01_bam" / "01_run_gtx_jobs.sh"
                        |canonical underscore name
        log: 可选日志器，回退时发 WARNING|optional logger, warned on fallback
        create: 仅用于语义标记；本函数不创建目录(创建由调用方 makedirs 完成)
                |semantic hint only; this function never creates anything

    Returns:
        str: 实际应使用的路径。新名已存在→新名；新名不存在但老名存在→老名(回退)；
             都不存在→新名(供调用方创建)
             |Path to use: new if exists; else old if exists (fallback); else
             new for creation.

    Examples:
        >>> # 续跑老任务：磁盘上是 01.bam
        >>> resolve_legacy_path('/out', '01_bam')
        '/out/01.bam'   # 回退老名，复用已有结果|fallback, reuse old output
        >>> # 全新运行：两个都不存在
        >>> resolve_legacy_path('/out', '01_bam')
        '/out/01_bam'   # 返回新名供 makedirs|return new for creation
    """
    parent_dir = str(parent_dir)
    new_path = os.path.join(parent_dir, canonical_name)
    if os.path.exists(new_path):
        return new_path

    # 回退候选：把首个 NN_ 换成 NN. 得到老名(扩展名里的点不受影响)
    # |Fallback candidate: replace the leading NN_ with NN. to get legacy name
    old_name = re.sub(r'^(\d+)_', r'\1.', canonical_name)
    if old_name != canonical_name:
        old_path = os.path.join(parent_dir, old_name)
        if os.path.exists(old_path):
            if log is not None:
                log.warning(
                    f"检测到旧命名目录/文件，回退使用|Legacy name fallback in use: "
                    f"{old_path} (规范名|canonical: {canonical_name})"
                )
            return old_path

    return new_path


def resolve_legacy_path_chain(parent_dir, *canonical_names):
    """多层步骤目录的兼容定位|Resolve a chained multi-level step path.

    逐层应用 resolve_legacy_path，每层优先下划线规范名、回退点号老名，
    用于读取跨模块/多层嵌套的历史输出（如 hifi_hic 的 03_ngs_polish/04_reassembly/...）。
    |Applies resolve_legacy_path at each level so deeply nested legacy output
    (e.g. hifi_hic's 03_ngs_polish/04_reassembly/...) is reused on resume.

    Args:
        parent_dir: 顶层父目录|top-level parent directory
        *canonical_names: 逐层规范名(下划线)，可含非编号段(如样本名，原样拼接)
                          |canonical name per level (non-numbered segments pass through)

    Returns:
        str: 解析后的完整路径|resolved full path

    Examples:
        >>> resolve_legacy_path_chain('/out', 'sample1', '03_ngs_polish', '02_fasta')
        '/out/sample1/03.ngs_polish/02.fasta'  # 若老名存在则逐层回退|per-level fallback
    """
    path = str(parent_dir)
    for name in canonical_names:
        path = resolve_legacy_path(path, name)
    return path


def get_tool_path(
    tool_name: str,
    default_path: str,
    env_var: Optional[str] = None
) -> str:
    """获取工具路径|Get tool path

    优先级|Priority:
    1. 环境变量|Environment variable
    2. 用户配置文件|User config file (~/.config/biopytools/config.yml)
    3. 默认路径（展开~）|Default path (expand ~)

    Args:
        tool_name: 工具名称|Tool name, 用于在配置文件中查找|for lookup in config file
        default_path: 默认路径|Default path, 支持~和$VAR|Supports ~ and $VAR
        env_var: 环境变量名|Environment variable name, 如果指定优先使用|if specified, use first

    Returns:
        str: 展开后的绝对路径|Expanded absolute path

    Examples:
        >>> # 使用环境变量|Use environment variable
        >>> os.environ['FANC_PATH'] = '~/custom/fanc/bin/fanc'
        >>> get_tool_path('fanc', '~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc', 'FANC_PATH')
        '/home/user/custom/fanc/bin/fanc'

        >>> # 使用配置文件|Use config file
        >>> get_tool_path('fanc', '~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc')
        '/home/user/custom/fanc/bin/fanc'  # 如果配置文件中设置了|if set in config

        >>> # 使用默认值|Use default
        >>> get_tool_path('fanc', '~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc')
        '/home/user/miniforge3/envs/fanc_v.0.9.23b/bin/fanc'
    """
    # 1. 检查环境变量|Check environment variable
    if env_var:
        env_path = os.getenv(env_var)
        if env_path:
            return expand_path(env_path)

    # 2. 检查用户配置文件|Check user config file
    config = load_user_config()
    if 'tools' in config and tool_name in config['tools']:
        return expand_path(config['tools'][tool_name])

    # 3. 使用默认值|Use default
    return expand_path(default_path)


def get_software_path(
    software_name: str,
    default_path: str,
    env_var: Optional[str] = None
) -> str:
    """获取软件路径（别名函数）|Get software path (alias function)

    与get_tool_path相同，提供更语义化的名称
    | Same as get_tool_path, provides more semantic name

    Args:
        software_name: 软件名称|Software name
        default_path: 默认路径|Default path
        env_var: 环境变量名|Environment variable name

    Returns:
        str: 展开后的绝对路径|Expanded absolute path
    """
    return get_tool_path(software_name, default_path, env_var)


def validate_tool_path(path: str, tool_name: str = "") -> bool:
    """验证工具路径是否存在|Validate if tool path exists

    Args:
        path: 工具路径|Tool path
        tool_name: 工具名称|Tool name (用于错误消息|for error message)

    Returns:
        bool: 是否存在|Whether exists

    Raises:
        ValueError: 如果路径不存在且提供了工具名称|If path doesn't exist and tool_name provided
    """
    expanded_path = expand_path(path)
    if not Path(expanded_path).exists():
        if tool_name:
            raise ValueError(
                f"工具不存在|Tool not found: {tool_name}\n"
                f"路径|Path: {expanded_path}\n"
                f"请检查|Please check:\n"
                f"  1. 设置环境变量|Set environment variable\n"
                f"  2. 或配置~/.config/biopytools/config.yml|Or configure ~/.config/biopytools/config.yml\n"
                f"  3. 或安装工具到默认位置|Or install tool to default location"
            )
        return False
    return True


# 常用工具路径快捷方式|Common tool path shortcuts
# 这些函数提供常用工具的默认路径|These functions provide default paths for common tools

def get_bwa_path(custom_path: Optional[str] = None) -> str:
    """获取BWA路径|Get BWA path"""
    if custom_path:
        return expand_path(custom_path)
    return get_tool_path('bwa', '~/.local/bin/bwa', 'BWA_PATH')


def get_samtools_path(custom_path: Optional[str] = None) -> str:
    """获取samtools路径|Get samtools path"""
    if custom_path:
        return expand_path(custom_path)
    return get_tool_path('samtools', '~/.local/bin/samtools', 'SAMTOOLS_PATH')


def get_fanc_path(custom_path: Optional[str] = None) -> str:
    """获取FAN-C路径|Get FAN-C path"""
    if custom_path:
        return expand_path(custom_path)
    return get_tool_path('fanc', '~/miniforge3/envs/fanc_v.0.9.23b/bin/fanc', 'FANC_PATH')
