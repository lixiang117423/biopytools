"""NCBI datasets 下载工具模块|NCBI datasets download utility module

日志管理器(stdout/stderr/文件三分离) + datasets 工具自动安装 + 命令执行
|Logger manager (stdout/stderr/file separation) + datasets auto-install + command runner
"""

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

# datasets CLI 官方下载地址(Linux x86_64)|Official datasets CLI download URL (Linux x86_64)
DATASETS_DOWNLOAD_URL = (
    'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
)


class ModuleLogger:
    """NCBI datasets 日志管理器|NCBI datasets Logger Manager

    stdout(<=INFO) + stderr(>=WARNING) + 三个日志文件:
    {base}.ncbi_datasets.log(全量) .out.log(<=INFO) .err.log(>=WARNING)
    |stdout (<=INFO) + stderr (>=WARNING) + three log files
    """

    def __init__(self, log_file: str, log_level: str = 'INFO'):
        self.log_file = log_file
        self.out_log_file = log_file.replace('.log', '.out.log')
        self.err_log_file = log_file.replace('.log', '.err.log')
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        self.logger = self._setup_logging(log_level)

    def _setup_logging(self, log_level: str) -> logging.Logger:
        """设置日志|Setup logging (named logger, 不污染 root|no root pollution)"""
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        level = getattr(logging, log_level.upper(), logging.INFO)

        logger = logging.getLogger('biopytools.ncbi_datasets')
        logger.handlers.clear()
        logger.propagate = False
        logger.setLevel(logging.DEBUG)

        # stdout: <=INFO(受 --log-level 控制)|stdout: <=INFO (respects --log-level)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(level)
        stdout_handler.addFilter(lambda r: r.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr: >=WARNING|stderr: >=WARNING
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件 handler|File handlers
        full_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        full_handler.setLevel(logging.DEBUG)
        full_handler.setFormatter(formatter)
        logger.addHandler(full_handler)

        out_handler = logging.FileHandler(self.out_log_file, encoding='utf-8')
        out_handler.setLevel(level)
        out_handler.addFilter(lambda r: r.levelno <= logging.INFO)
        out_handler.setFormatter(formatter)
        logger.addHandler(out_handler)

        err_handler = logging.FileHandler(self.err_log_file, encoding='utf-8')
        err_handler.setLevel(logging.WARNING)
        err_handler.setFormatter(formatter)
        logger.addHandler(err_handler)

        return logger

    def get_logger(self) -> logging.Logger:
        """获取 logger 实例|Get logger instance"""
        return self.logger


def format_number(num: int) -> str:
    """大数字格式化(>1百万用 M 单位保留2位小数)|Format large numbers (M unit above 1M)

    Examples:
        >>> format_number(10000000)
        '10.00M'
        >>> format_number(1234)
        '1234'
    """
    if num >= 1_000_000:
        return f'{num / 1_000_000:.2f}M'
    return str(num)


def is_tool_ready(tool_path: str) -> bool:
    """检查工具是否已安装且可执行|Check tool is installed and executable"""
    path = Path(tool_path)
    return path.is_file() and os.access(path, os.X_OK)


def install_datasets_tool(tool_path: str, logger) -> Optional[str]:
    """下载安装 datasets 二进制|Download and install the datasets binary

    从 NCBI 官方地址 curl 到目标路径并 chmod +x(§13 无 conda 依赖,直接调用)
    |curl from the official NCBI URL to the target path, then chmod +x
    """
    dest = Path(tool_path)
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        logger.error(f"无法创建工具目录|Cannot create tool directory: {dest.parent} - {e}")
        return None

    cmd = ['curl', '-L', '-o', str(dest), DATASETS_DOWNLOAD_URL]
    logger.info("执行|Executing: 下载 datasets 工具|datasets tool download")
    logger.info(f"命令|Command: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        logger.error("未找到 curl 命令|curl command not found")
        return None
    except OSError as e:
        logger.error(f"curl 执行失败|curl execution failed: {e}")
        return None

    if result.returncode != 0:
        logger.error(f"datasets 下载失败|datasets download failed: "
                     f"{result.stderr.strip() or result.stdout.strip()}")
        return None
    try:
        os.chmod(dest, 0o755)
    except OSError as e:
        logger.error(f"chmod 失败|chmod failed: {e}")
        return None
    logger.info(f"datasets 已安装到|datasets installed to: {dest}")
    return str(dest)


def ensure_datasets_tool(tool_path: str, logger) -> Optional[str]:
    """确保 datasets 工具可用,缺失时自动下载安装|Ensure datasets is ready, auto-install if missing

    Returns:
        str: 可用的工具路径|usable tool path, None 表示失败|None on failure
    """
    if is_tool_ready(tool_path):
        logger.info(f"检测到 datasets 工具|datasets tool found: {tool_path}")
        return tool_path
    logger.warning(
        f"未找到 datasets 工具,尝试自动下载|datasets tool not found, "
        f"attempting auto-install: {tool_path}")
    installed = install_datasets_tool(tool_path, logger)
    if installed and is_tool_ready(installed):
        return installed
    return None


def run_command(cmd: List[str], logger, capture: bool = False):
    """执行外部命令并记录完整命令|Run external command with full command logging (§2.2.1)

    Args:
        cmd: 命令列表|command list
        logger: 日志器|logger
        capture: 捕获 stdout/stderr(供解析)|capture output for parsing

    Returns:
        subprocess.CompletedProcess 或 None(命令不存在/启动失败)|or None on failure
    """
    logger.info(f"命令|Command: {' '.join(cmd)}")
    try:
        return subprocess.run(cmd, capture_output=capture, text=True)
    except FileNotFoundError as e:
        logger.error(f"命令不存在|Command not found: {e}")
        return None
    except OSError as e:
        logger.error(f"命令执行失败|Command execution failed: {e}")
        return None
