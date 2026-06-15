"""
EviAnn工具函数模块|EviAnn Utilities Module
"""

import logging
import sys
import os
from pathlib import Path
from typing import List


def expand_path(path: str) -> str:
    """
    展开路径中的~和环境变量|Expand ~ and environment variables in path

    Args:
        path: 原始路径（可包含~或$VAR）|Raw path (may contain ~ or $VAR)

    Returns:
        展开后的绝对路径|Expanded absolute path
    """
    return os.path.expandvars(os.path.expanduser(path))


class EviAnnLogger:
    """EviAnn日志管理器|EviAnn Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """
        设置日志|Setup logging

        Args:
            log_level: 日志级别|Log level
        """
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        # 配置handlers|Configure handlers
        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """
        获取日志器|Get logger

        Returns:
            logger对象|Logger object
        """
        return self.logger


def get_conda_env(command: str):
    """
    检测命令是否在conda环境中，返回环境名称|Detect if command is in conda environment

    Args:
        command: 命令名称或路径|Command name or path

    Returns:
        conda环境名称或None|Conda environment name or None
    """
    import shutil
    import re

    # 方法1: 从命令路径检测|Method 1: Detect from command path
    cmd_path = shutil.which(command)
    if cmd_path:
        # 检查路径中是否包含 envs|Check if path contains envs
        match = re.search(r'/envs/([^/]+)', cmd_path)
        if match:
            return match.group(1)

    return None


def build_conda_command(command: str, args: List[str]) -> List[str]:
    """
    构建conda run命令来运行conda环境中的软件|Build conda run command

    Args:
        command: 命令名称或完整路径|Command name or full path
        args: 命令参数列表|Command argument list

    Returns:
        完整命令列表|Full command list

    注意|Note:
        必须使用--no-capture-output避免conda缓冲输出导致内存问题
        Must use --no-capture-output to avoid conda buffering output causing memory issues
    """
    conda_env = get_conda_env(command)

    if conda_env:
        # 使用conda run调用|Use conda run
        full_cmd = ['conda', 'run', '-n', conda_env, '--no-capture-output', command] + args
    else:
        # 非conda环境，直接调用|Non-conda environment, call directly
        full_cmd = [command] + args

    return full_cmd
