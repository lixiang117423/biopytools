"""
CentIER工具函数模块|CentIER Utility Functions Module
"""

import logging
import sys
import os
from pathlib import Path
from typing import List
import subprocess


class CentIERLogger:
    """CentIER日志管理器|CentIER Logger Manager"""

    def __init__(self, log_file: Path, log_level: str = "INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.log_file = log_file
        self.log_level = log_level
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 删除已存在的日志文件|Remove existing log file
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, self.log_level.upper(), logging.INFO)

        # stdout handler - INFO级别|stdout handler - INFO level
        # → 超算系统捕获到 .out 文件|→ Captured by job scheduler to .out file
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_formatter = logging.Formatter(log_format, datefmt=date_format)
        stdout_handler.setFormatter(stdout_formatter)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        # → 超算系统捕获到 .err 文件|→ Captured by job scheduler to .err file
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_formatter = logging.Formatter(log_format, datefmt=date_format)
        stderr_handler.setFormatter(stderr_formatter)

        # 文件handler - 所有级别|File handler - all levels
        # → 本地完整日志|→ Local complete log
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(log_format, datefmt=date_format)
        file_handler.setFormatter(file_formatter)

        # 配置日志|Configure logging
        self.logger = logging.getLogger("CentIER")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False  # 避免重复|Avoid duplicates

        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(stderr_handler)
        self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def check_dependencies(config, logger) -> bool:
    """
    检查依赖工具|Check dependencies

    Args:
        config: CentIERConfig配置对象|CentIERConfig object
        logger: 日志器|Logger

    Returns:
        bool: 是否所有依赖都可用|Whether all dependencies are available
    """
    logger.info("检查依赖工具|Checking dependencies")

    all_ok = True

    # 检查centIER.py脚本|Check centIER.py script
    centier_script = config.get_centier_script_path()
    if not os.path.exists(centier_script):
        logger.error(f"未找到centIER.py脚本|centIER.py script not found: {centier_script}")
        all_ok = False
    else:
        logger.info(f"找到centIER.py|Found centIER.py: {centier_script}")

    # 检查bin目录下的工具|Check tools in bin directory
    bin_path = config.get_bin_path()

    required_tools = [
        ('hmmsearch', 'bin/hmmsearch'),
        ('ltr_finder', 'bin/ltr_finder/ltr_finder'),
        ('trf', 'bin/trf409.linux64'),
        ('REXdb.hmm', 'bin/REXdb.hmm'),
        ('Ty3_gypsy.hmm', 'bin/Ty3_gypsy.hmm')
    ]

    for tool_name, tool_rel_path in required_tools:
        tool_path = os.path.join(config.centier_path, tool_rel_path)
        if os.path.exists(tool_path):
            logger.info(f"找到|Found {tool_name}: {tool_path}")
        else:
            logger.warning(f"未找到|Warning: {tool_name} not found: {tool_path}")

    # 检查外部工具|Check external tools
    external_tools = {
        'gt': 'genometools (可选, 可跳过|optional, can skip)',
        'LTR_retriever': 'LTR_retriever (可选, 可跳过|optional, can skip)'
    }

    for tool, description in external_tools.items():
        tool_path = os.path.expanduser(f"~/miniforge3/envs/centier/bin/{tool}")
        if os.path.exists(tool_path):
            logger.info(f"找到|Found {tool}: {tool_path}")
        else:
            logger.warning(f"未找到|Warning: {tool} not found ({description})")

    if all_ok:
        logger.info("依赖检查完成|Dependency check completed")
    else:
        logger.error("依赖检查失败|Dependency check failed")

    return all_ok


def run_command(cmd: List[str], logger, check: bool = True) -> subprocess.CompletedProcess:
    """
    运行命令|Run command

    Args:
        cmd: 命令列表|Command list
        logger: 日志器|Logger
        check: 是否检查返回码|Whether to check return code

    Returns:
        subprocess.CompletedProcess: 命令执行结果|Command execution result
    """
    logger.debug(f"执行命令|Running command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=check
        )
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"命令执行失败|Command execution failed: {e}")
        logger.error(f"标准错误|Stderr: {e.stderr}")
        raise


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)
