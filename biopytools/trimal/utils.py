"""trimal 工具函数模块|trimal Utility Functions Module

trimal 是独立编译的 C++ 二进制,其动态库(libstdc++/libgcc_s)经 RPATH 解析,
无需 conda 环境激活即可直接调用绝对路径执行(实测 0.01s)。
|trimal is a standalone compiled C++ binary; its shared libs resolve via RPATH,
so it can be invoked directly by absolute path without `conda run` (verified 0.01s).
"""

import logging
import os
import subprocess
import sys
from typing import List, Tuple


class TrimalLogger:
    """trimal 日志管理器|trimal Logger Manager"""

    def __init__(self, log_file=None, verbose=False):
        self.log_file = log_file
        self.setup_logging(verbose)

    def setup_logging(self, verbose):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(log_format, datefmt=date_format)

        level = logging.DEBUG if verbose else logging.INFO

        logger = logging.getLogger("trimal")
        logger.setLevel(level)
        logger.handlers.clear()
        logger.propagate = False

        # stdout handler - INFO 级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING 及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件 handler - 所有级别|File handler - all levels
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        self.logger = logger

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def run_command(cmd: List[str], logger: logging.Logger, description: str = "") -> Tuple[str, str]:
    """
    执行命令(传入列表,shell=False)|Execute command (list, shell=False)

    trimal 为独立二进制,直接以绝对路径调用(不经 conda run 包装)。
    |trimal is a standalone binary, invoked directly by absolute path (no conda run wrapper).

    Args:
        cmd: 命令列表(首元素为 trimal 绝对路径)|Command list (first elem is trimal abs path)
        logger: 日志器|Logger
        description: 步骤描述|Step description
    Returns:
        (stdout, stderr)
    """
    if description:
        logger.info(f"执行|Executing: {description}")
    logger.info(f"命令|Command: {' '.join(cmd)}")

    result = subprocess.run(cmd, shell=False, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"命令执行失败|Command failed (exit={result.returncode}): {result.stderr.strip()}"
        )

    return result.stdout, result.stderr
