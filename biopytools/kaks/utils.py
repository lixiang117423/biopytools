"""
Ka/Ks Calculator工具函数模块|Ka/Ks Calculator Utility Functions Module
"""

import logging
import sys
import os
from pathlib import Path
from typing import List, Optional

# conda包装统一走公共层(§13): 同源conda绝对路径 + run -p <环境前缀>, 严禁裸调conda
from ..common.conda_runner import build_conda_command


class KaKsLogger:
    """Ka/Ks分析日志管理器|Ka/Ks Analysis Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "kaks_analysis.log", verbose: bool = False):
        """
        初始化日志器|Initialize logger

        Args:
            output_dir: 输出目录(日志自动放到99_logs/子目录)|Output directory (log auto-placed in 99_logs/)
            log_name: 日志文件名|Log file name
            verbose: 详细模式|Verbose mode
        """
        self.output_dir = output_dir
        self.log_dir = output_dir / "99_logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.log_file = self.log_dir / log_name
        self._verbose = verbose
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        logger = logging.getLogger("KaKsAnalyzer")
        logger.setLevel(logging.DEBUG if self._verbose else logging.INFO)
        logger.handlers.clear()
        logger.propagate = False

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        self.logger = logger

    def info(self, message: str):
        """信息日志|Info logging"""
        self.logger.info(message)

    def success(self, message: str):
        """成功日志|Success logging"""
        self.logger.info(message)

    def warning(self, message: str):
        """警告日志|Warning logging"""
        self.logger.warning(message)

    def error(self, message: str):
        """错误日志|Error logging"""
        self.logger.error(message)

    def debug(self, message: str):
        """调试日志|Debug logging"""
        self.logger.debug(message)

    def progress(self, message: str, current: int, total: int):
        """进度日志|Progress logging"""
        percentage = (current / total) * 100 if total > 0 else 0
        self.logger.info(f"{message} [{current}/{total}] ({percentage:.1f}%)")

    def separator(self, title: str = ""):
        """分隔符日志|Separator logging"""
        if title:
            self.logger.info(f"{title} " + "="*50)
        else:
            self.logger.info("="*60)

    def get_logger(self):
        """获取日志器对象|Get logger object"""
        return self.logger


def format_number(num: int) -> str:
    """格式化大数字|Format large number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return f"{num:,}"
