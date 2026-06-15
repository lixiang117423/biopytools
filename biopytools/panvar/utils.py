"""泛基因组变异分析工具函数模块|Pan-genome Variant Analysis Utility Functions Module"""

import logging
import sys
from pathlib import Path


class PanvarLogger:
    """泛基因组变异分析日志管理器|Pan-genome Variant Analysis Logger Manager"""

    def __init__(self, log_file: str = None, log_level: str = "INFO"):
        """
        初始化日志|Initialize logger

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG, INFO, WARNING, ERROR)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level: str):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        # 创建logger|Create logger
        self.logger = logging.getLogger("Panvar")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers = []
        self.logger.propagate = False

        formatter = logging.Formatter(log_format, datefmt=date_format)

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler - 所有级别|File handler - all levels
        if self.log_file:
            Path(self.log_file).parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(self.log_file, mode='w', encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self) -> logging.Logger:
        """获取logger对象|Get logger object"""
        return self.logger
