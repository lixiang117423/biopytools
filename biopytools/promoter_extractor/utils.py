"""
启动子提取器工具函数模块|Promoter Extractor Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional


class PromoterLogger:
    """启动子提取器日志管理器|Promoter Extractor Logger Manager"""

    def __init__(self, output_path: Path, log_name: str = "promoter_extractor.log",
                 log_level: str = "INFO", quiet: bool = False):
        """
        初始化日志管理器|Initialize logger manager

        Args:
            output_path: 输出目录|Output directory
            log_name: 日志文件名|Log file name
            log_level: 日志级别|Log level (DEBUG/INFO/WARNING/ERROR/CRITICAL)
            quiet: 静默模式|Quiet mode (only ERROR)
        """
        self.output_path = output_path
        self.log_file = output_path / log_name
        self.log_level = log_level
        self.quiet = quiet

        # 创建logger|Create logger
        self.logger = logging.getLogger("PromoterExtractor")
        self.logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))

        # 清除现有handlers|Clear existing handlers
        self.logger.handlers.clear()

        # 文件handler|File handler
        file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)

        # 控制台handler|Console handler
        if not quiet:
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(getattr(logging, log_level.upper(), logging.INFO))
            console_formatter = logging.Formatter('%(message)s')
            console_handler.setFormatter(console_formatter)
            self.logger.addHandler(console_handler)

    def get_logger(self):
        """获取logger对象|Get logger object"""
        return self.logger


class CommandRunner:
    """命令执行器|Command Runner"""

    def __init__(self, logger: logging.Logger, output_path: Path):
        """
        初始化命令执行器|Initialize command runner

        Args:
            logger: 日志对象|Logger object
            output_path: 输出目录|Output directory
        """
        self.logger = logger
        self.output_path = output_path

    def run_command(self, cmd: list, description: str = "") -> bool:
        """
        执行命令|Execute command

        Args:
            cmd: 命令列表|Command list
            description: 命令描述|Command description

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            self.logger.debug(f"执行命令|Executing: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode != 0:
                self.logger.error(f"{description}失败|{description} failed: {result.stderr}")
                return False

            if result.stdout:
                self.logger.debug(result.stdout)

            return True

        except Exception as e:
            self.logger.error(f"执行命令异常|Command execution error: {e}")
            return False


def parse_gene_list(gene_list_file: str) -> set:
    """
    解析基因列表文件|Parse gene list file

    Args:
        gene_list_file: 基因列表文件路径|Gene list file path
                       (每行一个基因ID|One gene ID per line)

    Returns:
        set: 基因ID集合|Gene ID set
    """
    genes = set()
    try:
        with open(gene_list_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    genes.add(line)
    except Exception as e:
        raise IOError(f"读取基因列表文件失败|Failed to read gene list file: {e}")

    return genes


def reverse_complement(sequence: str) -> str:
    """
    生成反向互补序列|Generate reverse complement sequence

    Args:
        sequence: DNA序列|DNA sequence

    Returns:
        str: 反向互补序列|Reverse complement sequence
    """
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n'
    }

    return ''.join(complement.get(base, base) for base in reversed(sequence))
