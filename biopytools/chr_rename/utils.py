"""
染色体重命名工具函数模块|Chromosome Rename Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path

from ..common.conda_runner import CommandRunner

from ..common.conda_runner import CommandRunner


class ChrRenameLogger:
    """染色体重命名日志管理器|Chromosome Rename Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "chr_rename.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()

        # 设置日志格式|Set log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # stdout handler|Stdout handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)

        # 配置日志|Configure logging
        logging.basicConfig(
            level=logging.DEBUG,
            handlers=[file_handler, stdout_handler]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


# 复用公共层 CommandRunner(自动记录完整命令到INFO, 支持conda环境包装)
# |Reuse common CommandRunner (auto full-command INFO logging, conda env wrapping)
