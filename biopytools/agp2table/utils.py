"""
AGP转表格工具函数模块|AGP to Table Utility Functions Module
"""

import logging
import sys
from pathlib import Path


class AGPLogger:
    """AGP转表格日志管理器|AGP to Table Logger Manager"""

    def __init__(self, log_file: Path):
        self.log_file = log_file
        self.setup_logging()

    def setup_logging(self):
        """设置日志系统（符合规范2.3）|Setup logging system (spec 2.3 compliant)"""
        # 删除旧日志|Delete old log
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler - DEBUG级别|File handler - DEBUG level
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler - INFO级别|stdout handler - INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # stderr handler - WARNING及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler, stderr_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger
