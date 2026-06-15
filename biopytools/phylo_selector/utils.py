"""
系统发育树样品选择工具函数模块|Phylogenetic Tree Sample Selector Utility Functions Module
"""

import logging
import sys
from pathlib import Path


class PhyloSelectorLogger:
    """系统发育树样品选择日志管理器|Phylogenetic Tree Sample Selector Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            # 确保日志目录存在|Ensure log directory exists
            log_path = Path(self.log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def format_number(num):
    """格式化数字|Format number

    大于1百万的数字使用M单位|Use M unit for numbers larger than 1 million

    Args:
        num: 数字|Number

    Returns:
        str: 格式化后的字符串|Formatted string
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)
