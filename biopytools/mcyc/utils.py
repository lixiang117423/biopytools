"""
甲烷循环基因丰度分析工具函数模块|Methane Cycle Gene Abundance Analysis Utility Functions Module
"""

import logging
import sys
from pathlib import Path


class MCycLogger:
    """甲烷循环分析日志管理器|Methane Cycle Analysis Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志器|Initialize logger

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level
        """
        self.log_file = log_file
        self.log_level = log_level
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        handlers = [logging.StreamHandler(sys.stdout)]

        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file, encoding='utf-8'))

        logging.basicConfig(
            level=getattr(logging, self.log_level),
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger

    def info(self, msg):
        """记录信息|Log info"""
        self.logger.info(msg)

    def warning(self, msg):
        """记录警告|Log warning"""
        self.logger.warning(msg)

    def error(self, msg):
        """记录错误|Log error"""
        self.logger.error(msg)

    def debug(self, msg):
        """记录调试信息|Log debug"""
        self.logger.debug(msg)
