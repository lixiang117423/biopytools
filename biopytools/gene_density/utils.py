"""基因密度计算工具函数|Gene density utility functions

日志管理器与通用辅助函数|Logger manager and general helpers
"""

import logging
import sys


class GeneDensityLogger:
    """基因密度日志管理器|Gene density logger manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging (标准格式§2.1)"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

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
        """获取日志器|Get logger"""
        return self.logger


def format_number(num: int) -> str:
    """格式化大数字(>=1M用M单位)|Format large number (>=1M uses M unit)

    Args:
        num: 整数|Integer

    Returns:
        str: 格式化后的字符串|Formatted string (e.g. "10.00M")

    示例|Examples:
        >>> format_number(10_000_000)
        '10.00M'
        >>> format_number(999_999)
        '999999'
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    return str(num)
