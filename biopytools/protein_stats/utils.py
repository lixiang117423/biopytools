"""
Protein Stats工具函数模块|Protein Stats Utility Functions Module
"""

import logging
import sys
from pathlib import Path


class ProteinStatsLogger:
    """Protein Stats日志管理器|Protein Stats Logger Manager"""

    def __init__(self, log_file: str = "protein_stats_analysis.log"):
        self.log_file = log_file
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger
