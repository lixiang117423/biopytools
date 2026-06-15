"""
VCF抽样工具模块|VCF Sampling Utilities Module
"""

import logging
import sys
from pathlib import Path


class VCFSamplerLogger:
    """VCF抽样日志管理器|VCF Sampler Logger Manager"""

    def __init__(self, log_file=None, log_level=logging.INFO):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径 (可选)|Log file path (optional)
            log_level: 日志级别|Log level
        """
        self.logger = logging.getLogger('vcf_sampler')
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()

        # 日志格式|Log format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # stdout handler - DEBUG到INFO级别|stdout handler - DEBUG to INFO level
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上级别|stderr handler - WARNING and above
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

        # 文件handler (如果指定)|File handler (if specified)
        if log_file:
            log_path = Path(log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(log_file, encoding='utf-8')
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_logger(self):
        """获取logger对象|Get logger object"""
        return self.logger
