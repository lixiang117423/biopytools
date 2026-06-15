"""
基因组装配统计工具函数模块|Genome Assembly Statistics Utility Functions Module
"""

import logging
import sys
from pathlib import Path


class AssemblyStatsLogger:
    """基因组装配统计日志管理器|Genome Assembly Statistics Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            log_file: 日志文件路径|Log file path
            log_level: 日志级别|Log level (DEBUG/INFO/WARNING/ERROR)
        """
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 配置日志格式|Configure logging format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        # 设置日志级别|Set log level
        level = getattr(logging, log_level.upper(), logging.INFO)

        # 配置handlers|Configure handlers
        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        # 配置logging|Configure logging
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


def parse_fasta(file_path: str, min_length: int = 1):
    """
    解析FASTA文件|Parse FASTA file

    Args:
        file_path: 文件路径|File path
        min_length: 最小序列长度|Minimum sequence length

    Yields:
        tuple: (header, sequence)
    """
    current_header = None
    current_sequence = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                # 保存前一个序列|Save previous sequence
                if current_header is not None:
                    sequence = ''.join(current_sequence)
                    if len(sequence) >= min_length:
                        yield (current_header, sequence)

                # 开始新序列|Start new sequence
                current_header = line[1:].split()[0]  # 取第一个词作为header
                current_sequence = []
            else:
                current_sequence.append(line)

        # 保存最后一个序列|Save last sequence
        if current_header is not None:
            sequence = ''.join(current_sequence)
            if len(sequence) >= min_length:
                yield (current_header, sequence)


def parse_fastq(file_path: str, min_length: int = 1):
    """
    解析FASTQ文件|Parse FASTQ file

    Args:
        file_path: 文件路径|File path
        min_length: 最小序列长度|Minimum sequence length

    Yields:
        tuple: (header, sequence)
    """
    with open(file_path, 'r') as f:
        while True:
            # 读取四行|Read four lines
            header = f.readline().strip()
            if not header:
                break

            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()

            if not sequence or not plus or not quality:
                break

            if header.startswith('@'):
                header = header[1:].split()[0]

                if len(sequence) >= min_length:
                    yield (header, sequence)
