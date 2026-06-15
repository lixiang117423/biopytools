"""
三代转录组比对工具函数模块|Long RNA-seq Alignment Utility Functions Module
"""

import logging
import sys
import subprocess
from pathlib import Path


class LongRNASeqLogger:
    """三代转录组比对日志管理器|Long RNA-seq Alignment Logger Manager"""

    def __init__(self, log_file: Path, log_name: str = "longrnaseq.log"):
        self.log_file = log_file
        self.log_name = log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def run_command(cmd: list, logger: logging.Logger = None,
                check: bool = True, capture_output: bool = True) -> subprocess.CompletedProcess:
    """
    运行命令|Run command

    Args:
        cmd: 命令列表|Command list
        logger: 日志器|Logger
        check: 是否检查返回码|Whether to check return code
        capture_output: 是否捕获输出|Whether to capture output

    Returns:
        subprocess.CompletedProcess: 命令执行结果|Command execution result
    """
    if logger:
        logger.debug(f"执行命令|Running command: {' '.join(cmd)}")

    return subprocess.run(
        cmd,
        capture_output=capture_output,
        text=True,
        check=check
    )


def check_tool(tool_path: str, tool_name: str, logger: logging.Logger = None) -> bool:
    """
    检查工具是否可用|Check if tool is available

    Args:
        tool_path: 工具路径|Tool path
        tool_name: 工具名称|Tool name
        logger: 日志器|Logger

    Returns:
        bool: 是否可用|Whether available
    """
    try:
        result = subprocess.run(
            [tool_path, "--version"],
            capture_output=True,
            text=True,
            check=False
        )
        version = result.stdout.split('\n')[0] if result.stdout else "unknown"
        if logger:
            logger.info(f"找到工具|Found tool: {tool_name} - {version}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        if logger:
            logger.error(f"未找到工具|Tool not found: {tool_name}")
            logger.error(f"请先安装{tool_name}|Please install {tool_name} first")
        return False


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)
