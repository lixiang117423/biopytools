"""
ps-gene-anno 工具函数|ps-gene-anno Utility Functions
日志管理器自写; 命令执行器/conda 检测复用 braker(已验证的通用实现, DRY)
|Logger is self-written; CommandRunner/conda detection reused from braker (DRY)
"""

import logging
import sys
from pathlib import Path

# 复用 braker 已验证的通用工具|Reuse braker's verified generic utilities (DRY)
from ..braker.utils import (
    CommandRunner, get_conda_env, format_number, check_step_completed
)

__all__ = [
    'PsGeneAnnoLogger', 'CommandRunner', 'get_conda_env',
    'format_number', 'check_step_completed',
]


class PsGeneAnnoLogger:
    """ps-gene-anno 日志管理器|Logger manager"""

    def __init__(self, log_file_path: str, log_level: str = "INFO"):
        self.log_file = Path(log_file_path)
        self.log_level = log_level
        self.setup_logging()

    def setup_logging(self):
        """设置日志(stdout INFO + file DEBUG)|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        level = getattr(logging, self.log_level.upper(), logging.INFO)
        formatter = logging.Formatter(log_format, datefmt=date_format)

        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(level)
        stdout_handler.setFormatter(formatter)

        logger = logging.getLogger('ps_gene_anno')
        logger.setLevel(logging.DEBUG)
        logger.handlers.clear()
        logger.propagate = False
        logger.addHandler(file_handler)
        logger.addHandler(stdout_handler)
        self.logger = logger

    def get_logger(self):
        return self.logger
