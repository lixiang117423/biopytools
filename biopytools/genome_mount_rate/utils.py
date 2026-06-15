"""
基因组挂载率统计工具函数模块|Genome Mount Rate Utility Functions Module
"""

import logging
import sys


class GenomeMountRateLogger:
    """基因组挂载率统计日志管理器|Genome Mount Rate Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path (可选|optional)
            log_level: 日志级别|Log level (默认|default: INFO)
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
            handlers.append(logging.FileHandler(self.log_file, encoding='utf-8'))

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
