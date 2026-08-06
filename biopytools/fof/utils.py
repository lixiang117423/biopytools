"""
FOF生成工具函数模块|FOF Generation Utility Functions Module
"""

import logging
import sys
from pathlib import Path
from typing import Optional


class FofLogger:
    """FOF生成日志管理器|FOF Generation Logger Manager"""

    def __init__(self, output_dir: Path, log_name: str = "fof.log"):
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


def extract_sample_name(filepath: Path) -> str:
    """
    从文件名提取样品名（去掉最后一层扩展名）|Extract sample name from filename (strip last extension)

    Args:
        filepath: 文件路径|File path

    Returns:
        样品名|Sample name

    Examples:
        >>> extract_sample_name(Path('/data/sample1_R1.fastq.gz'))
        'sample1_R1'
        >>> extract_sample_name(Path('/data/sample1.bam'))
        'sample1'
    """
    filename = filepath.name
    # 去掉最后一层扩展名|Strip the last extension
    # 例如: sample1_R1.fastq.gz -> sample1_R1
    dot_pos = filename.rfind('.')
    if dot_pos > 0:
        return filename[:dot_pos]
    return filename


def format_number(num: int) -> str:
    """格式化数字|Format number"""
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)
