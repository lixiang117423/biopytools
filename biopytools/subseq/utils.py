"""
序列子集提取工具函数模块|Sequence Subsequence Extraction Utility Functions Module
"""

import logging
import os
import sys
from pathlib import Path


class SubseqLogger:
    """序列子集提取日志管理器|Sequence Subsequence Extraction Logger Manager"""

    def __init__(self, output_dir, log_name: str = "subseq_extraction.log"):
        # 确保output_dir是Path对象|Ensure output_dir is Path object
        self.output_dir = Path(output_dir) if not isinstance(output_dir, Path) else output_dir
        self.log_file = self.output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 确保输出目录存在|Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

        if self.log_file.exists():
            self.log_file.unlink()

        # 日志格式|Log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        # 创建logger|Create logger
        self.logger = logging.getLogger(f"subseq_{id(self)}")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers.clear()
        self.logger.propagate = False  # 避免重复|Avoid duplicates

        formatter = logging.Formatter(log_format, datefmt=date_format)

        # 文件handler - 所有级别|File handler - all levels
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

        # stdout handler - INFO级别|stdout handler - INFO level
        # → 超算系统捕获到 .out 文件|→ Captured by job scheduler to .out file
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setLevel(logging.INFO)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

        # stderr handler - WARNING及以上|stderr handler - WARNING and above
        # → 超算系统捕获到 .err 文件|→ Captured by job scheduler to .err file
        stderr_handler = logging.StreamHandler(sys.stderr)
        stderr_handler.setLevel(logging.WARNING)
        stderr_handler.setFormatter(formatter)
        self.logger.addHandler(stderr_handler)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def validate_file_exists(file_path: str, description: str = "文件"):
    """验证文件是否存在|Validate file exists"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{description}不存在|{description} does not exist: {file_path}")
    return True


def validate_id_list(id_list_file: str, logger):
    """验证ID列表文件格式|Validate ID list file format"""
    try:
        with open(id_list_file, 'r') as f:
            id_list = [line.strip() for line in f if line.strip()]

        if not id_list:
            logger.warning("ID列表文件为空|ID list file is empty")
            return False

        logger.info(f"从文件中读取到 {len(id_list)} 个ID|Read {len(id_list)} IDs from file")
        return id_list

    except Exception as e:
        logger.error(f"读取ID列表文件失败|Failed to read ID list file: {e}")
        return False