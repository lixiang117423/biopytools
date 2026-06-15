"""
基因组Gap统计工具函数模块|Genome Gap Statistics Utility Functions Module
"""

import logging
import sys
from pathlib import Path


class GapStatLogger:
    """基因组Gap统计日志管理器|Genome Gap Statistics Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """初始化日志管理器|Initialize logger manager"""
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志系统（符合规范2.3）|Setup logging system (spec 2.3 compliant)"""
        # 删除旧日志|Delete old log
        if self.log_file and Path(self.log_file).exists():
            Path(self.log_file).unlink()

        # 标准日志格式|Standard log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = []

        # 文件handler - DEBUG级别|File handler - DEBUG level
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))
            handlers.append(file_handler)

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))
        handlers.append(stdout_handler)

        # stderr handler - WARNING及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(logging.Formatter(log_format, datefmt=date_format))
        handlers.append(stderr_handler)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=level,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger
