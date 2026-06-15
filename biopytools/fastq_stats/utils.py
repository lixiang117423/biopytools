"""
FASTQ文件统计工具模块|FASTQ File Statistics Utils Module
"""

import logging
import sys
import subprocess
import re
import glob
import os
from typing import Dict, List, Tuple, Optional


class FastqStatsLogger:
    """FASTQ文件统计日志管理器|FASTQ File Statistics Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        """初始化日志管理器|Initialize logger manager"""
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


class SeqkitChecker:
    """Seqkit工具检查器|Seqkit Tool Checker"""

    @staticmethod
    def check_seqkit() -> bool:
        """检查seqkit是否已安装|Check if seqkit is installed"""
        try:
            result = subprocess.run(
                ['seqkit', 'version'],
                capture_output=True,
                text=True
            )
            return result.returncode == 0
        except FileNotFoundError:
            return False


class FastqFileFinder:
    """FASTQ文件查找器|FASTQ File Finder"""

    # 支持的FASTQ扩展名|Supported FASTQ extensions
    FASTQ_EXTENSIONS = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')

    @staticmethod
    def pattern_to_regex(pattern: str) -> re.Pattern:
        """将用户提供的模式转换为正则表达式|Convert user pattern to regex"""
        escaped = re.escape(pattern)
        regex = escaped.replace(r'\*', '(.+)')
        return re.compile(regex)

    def find_paired_files(self, input_path: str, pattern: Optional[str] = None) -> Dict[str, Dict[str, str]]:
        """
        根据模式查找配对的FASTQ文件|Find paired FASTQ files by pattern

        Args:
            input_path: 输入文件或目录|Input file or directory
            pattern: 文件匹配模式|File matching pattern (e.g., "*_1.clean.fq.gz")

        Returns:
            字典格式: {sample_name: {'R1': r1_path, 'R2': r2_path}}
            Dict format: {sample_name: {'R1': r1_path, 'R2': r2_path}}
        """
        if not pattern:
            # 如果没有指定模式，查找所有fastq文件
            # If no pattern specified, find all fastq files
            return self._find_all_fastq_files(input_path)

        # 使用模式匹配|Use pattern matching
        return self._find_files_by_pattern(input_path, pattern)

    def _find_all_fastq_files(self, input_path: str) -> Dict[str, Dict[str, str]]:
        """查找所有FASTQ文件|Find all FASTQ files"""
        files = {}

        if os.path.isfile(input_path):
            # 单个文件|Single file
            basename = os.path.basename(input_path)
            files[basename] = {'R1': input_path}
        elif os.path.isdir(input_path):
            # 目录|Directory
            for f in os.listdir(input_path):
                if f.endswith(self.FASTQ_EXTENSIONS):
                    filepath = os.path.join(input_path, f)
                    files[f] = {'R1': filepath}

        return files

    def _find_files_by_pattern(self, input_path: str, pattern: str) -> Dict[str, Dict[str, str]]:
        """根据模式查找文件|Find files by pattern"""
        # 确定搜索路径|Determine search path
        if os.path.isfile(input_path):
            search_dir = os.path.dirname(input_path) or '.'
            search_pattern = os.path.basename(input_path)
        else:
            search_dir = input_path
            search_pattern = pattern

        # 查找匹配R1的文件|Find R1 matching files
        search_path = os.path.join(search_dir, search_pattern)
        r1_files = glob.glob(search_path)

        if not r1_files:
            return {}

        samples = {}
        regex = self.pattern_to_regex(pattern)

        for r1_file in r1_files:
            filename = os.path.basename(r1_file)
            match = regex.match(filename)

            if match:
                sample_name = match.group(1)

                # 尝试查找对应的R2文件|Try to find corresponding R2 file
                r2_pattern = self._generate_r2_pattern(pattern)
                r2_search = os.path.join(search_dir, r2_pattern.replace('*', sample_name))
                r2_files = glob.glob(r2_search)

                sample_data = {'R1': r1_file}
                # 只有当R2文件与R1文件不同时才添加R2
                # Only add R2 if it's different from R1
                if r2_files and r2_files[0] != r1_file:
                    sample_data['R2'] = r2_files[0]

                samples[sample_name] = sample_data

        return samples

    @staticmethod
    def _generate_r2_pattern(pattern: str) -> str:
        """生成R2文件模式|Generate R2 file pattern"""
        r2_pattern = pattern
        # 常见的R1到R2转换模式|Common R1 to R2 conversion patterns
        r2_pattern = r2_pattern.replace('_1.', '_2.')
        r2_pattern = r2_pattern.replace('_R1.', '_R2.')
        r2_pattern = r2_pattern.replace('.1.', '.2.')
        r2_pattern = r2_pattern.replace('.R1.', '.R2.')
        return r2_pattern


def format_number(num: int) -> str:
    """
    格式化数字|Format number

    大于1百万使用M单位显示|Display numbers > 1 million with M unit
    """
    if num >= 1_000_000:
        return f"{num / 1_000_000:.2f}M"
    elif num >= 1_000:
        return f"{num / 1_000:.2f}K"
    return str(num)


def format_size(size_bytes: int) -> str:
    """
    格式化文件大小|Format file size

    Args:
        size_bytes: 字节数|Size in bytes

    Returns:
        格式化的大小字符串|Formatted size string
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f}{unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f}PB"


def get_file_size(file_path: str) -> tuple:
    """
    获取文件大小（压缩和解压后）|Get file size (compressed and uncompressed)

    Args:
        file_path: 文件路径|File path

    Returns:
        (文件大小字节, 解压后大小字节)|(file_size_bytes, uncompressed_size_bytes)
        如果无法计算解压后大小，返回(file_size, None)|If cannot calculate uncompressed, return (file_size, None)
    """
    # 获取文件大小|Get file size
    file_size = os.path.getsize(file_path)

    # 如果不是gzip文件，返回相同值|If not gzip, return same value
    if not file_path.endswith('.gz'):
        return file_size, file_size

    # 对于gzip文件，尝试获取解压后大小|For gzip files, try to get uncompressed size
    try:
        import gzip
        with gzip.open(file_path, 'rb') as f:
            # 读取所有内容获取大小|Read all to get size
            uncompressed_size = len(f.read())
        return file_size, uncompressed_size
    except Exception:
        # 如果读取失败，只返回文件大小|If read fails, only return file size
        return file_size, None
