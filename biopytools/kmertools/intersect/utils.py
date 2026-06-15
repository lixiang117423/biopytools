"""
工具函数模块|Utility Functions Module
"""

import logging
import sys


class KmerIntersectLogger:
    """Kmer交集日志管理器|Kmer Intersection Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level
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
            handlers.append(logging.FileHandler(self.log_file))

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


def reverse_complement(dna_sequence: str) -> str:
    """
    计算DNA序列的反向互补序列|Calculate reverse complement of DNA sequence

    Args:
        dna_sequence: DNA序列|DNA sequence

    Returns:
        str: 反向互补序列|Reverse complement sequence

    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'
        >>> reverse_complement("AANNN")
        'NNNTT'
    """
    # 互补碱基映射|Complement base mapping
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n'  # 模糊碱基|Ambiguous bases
    }

    # 反转序列并互补|Reverse sequence and complement
    try:
        reverse_comp = ''.join(complement.get(base, base) for base in reversed(dna_sequence))
        return reverse_comp
    except Exception as e:
        raise ValueError(f"计算反向互补序列失败|Failed to calculate reverse complement: {e}")


def format_number(num: int) -> str:
    """
    格式化数字|Format number

    Args:
        num: 数字|Number

    Returns:
        str: 格式化后的字符串|Formatted string

    Examples:
        >>> format_number(1000000)
        '1.00M'
        >>> format_number(1500)
        '1.50K'
        >>> format_number(999)
        '999'
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    else:
        return str(num)


def print_interval_progress(current: int, interval: int):
    """
    打印进度信息|Print progress information

    Args:
        current: 当前值|Current value
        interval: 间隔|Interval
    """
    print(f"\r已处理|Processed: {format_number(current)} 行|lines", end='', flush=True)
    if current % (interval * 10) == 0:
        print()  # 每10个间隔换行一次|New line every 10 intervals
