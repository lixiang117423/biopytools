"""
配置管理模块|Configuration Management Module
"""

import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .utils import open_input


@dataclass
class Kmer2VcfConfig:
    """Kmer转VCF配置类|Kmer to VCF Configuration Class"""

    # 输入文件|Input files
    input_matrix: str
    output_vcf: str

    # 模式选择|Mode selection
    fast_mode: bool = True  # 快速模式（默认）|Fast mode (default)
    standard_mode: bool = False  # 标准模式（3遍处理）|Standard mode (3-pass processing)

    # 快速模式参数|Fast mode parameters
    chr_length: int = 100000000  # 染色体长度|Chromosome length (default: 100M), 仅当chr_number未设置时有效|Only used when chr_number is not set
    chr_number: int = 0  # 染色体数量|Number of chromosomes (default: 0=auto), 如果设置则优先使用|If set, this takes priority over chr_length
    min_freq: int = 0  # 最小出现频次过滤（快速模式）|Minimum frequency filter (fast mode)
    kmer_length: int = 51  # kmer长度|Kmer length
    no_header: bool = False  # 输入文件无header行|Input file has no header line

    # 标准模式过滤参数|Standard mode filter parameters
    min_agg_count: int = 3  # 最小聚合频次|Minimum aggregated count (standard mode)

    # 多线程参数|Multithreading parameters
    threads: int = 12  # 线程数|Number of threads
    chunk_size: int = 100000  # 每块行数|Lines per chunk for streaming

    # 临时文件目录|Temporary file directory
    temp_dir: Optional[str] = None

    # 内部属性|Internal attributes
    samples: list = None
    num_samples: int = 0
    _temp_dir_obj = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建临时目录|Create temporary directory
        if self.temp_dir:
            self.temp_dir = os.path.normpath(os.path.abspath(self.temp_dir))
            os.makedirs(self.temp_dir, exist_ok=True)
        else:
            # 默认使用当前工作目录下的temp文件夹，避免/tmp空间不足
            # Use temp folder in current working directory by default to avoid /tmp space issues
            default_temp = os.path.join(os.getcwd(), 'temp')
            self.temp_dir = os.path.normpath(os.path.abspath(default_temp))
            os.makedirs(self.temp_dir, exist_ok=True)
            self._temp_dir_obj = None  # 不需要自动清理，保留临时文件供检查|Don't auto-clean, keep files for inspection

        # 标准化路径|Normalize paths
        self.input_matrix = os.path.normpath(os.path.abspath(self.input_matrix))
        self.output_vcf = os.path.normpath(os.path.abspath(self.output_vcf))

        # 创建输出目录|Create output directory
        output_dir = os.path.dirname(self.output_vcf)
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)

        # 解析样本名|Parse sample names
        self._parse_samples()

    def _parse_samples(self):
        """解析样本名|Parse sample names from first line"""
        if not os.path.exists(self.input_matrix):
            return

        with open_input(self.input_matrix) as f:
            first_line = f.readline().strip()

            # 尝试制表符分隔|Try tab delimiter
            parts = first_line.split('\t')

            # 如果制表符分隔结果少于2列，尝试空格分隔|If tab delimiter gives <2 columns, try space
            if len(parts) < 2:
                parts = first_line.split()

            # 检查第一列是否是KMER标识符|Check if first column is KMER identifier
            if len(parts) > 0 and parts[0] == 'KMER':
                # 第一列是kmer，后面都是样本|First column is kmer, rest are samples
                self.samples = parts[1:]
                self.num_samples = len(self.samples)
            elif len(parts) > 1:
                # 没有KMER列，所有列都是样本（旧格式兼容）|No KMER column, all columns are samples (legacy format compatibility)
                self.samples = parts
                self.num_samples = len(self.samples)
            else:
                self.samples = []
                self.num_samples = 0

    def get_output_txt_path(self):
        """获取TXT输出路径|Get TXT output path"""
        # 从VCF路径推导TXT路径|Derive TXT path from VCF path
        if self.output_vcf.endswith('.vcf.gz'):
            txt_path = self.output_vcf.replace('.vcf.gz', '.txt')
        elif self.output_vcf.endswith('.vcf'):
            txt_path = self.output_vcf.replace('.vcf', '.txt')
        else:
            txt_path = self.output_vcf + '.txt'
        return txt_path

    def get_log_path(self):
        """获取日志路径|Get log path"""
        # 日志放在输出文件同目录|Log file in same directory as output
        output_dir = os.path.dirname(self.output_vcf)
        if not output_dir:
            output_dir = '.'
        return os.path.join(output_dir, "kmer2vcf.log")

    def is_compressed_output(self):
        """判断是否输出压缩VCF|Check if output should be compressed"""
        return self.output_vcf.endswith('.gz')

    def get_temp_file_path(self, filename):
        """获取临时文件路径|Get temporary file path"""
        return os.path.join(self.temp_dir, filename)

    def cleanup(self):
        """清理临时文件|Clean up temporary files"""
        if self._temp_dir_obj:
            try:
                self._temp_dir_obj.cleanup()
            except Exception:
                pass

    def should_use_bucket_mode(self, threshold_lines: int = 100_000_000) -> bool:
        """
        判断是否应该使用分桶模式|Determine if bucket mode should be used

        对于超大数据集（>1亿行），使用分桶模式可以显著减少内存占用
        For ultra-large datasets (>100M rows), bucket mode significantly reduces memory usage

        Args:
            threshold_lines: 判断阈值行数|Threshold line count (default: 100M)

        Returns:
            bool: 是否使用分桶模式|Whether to use bucket mode
        """
        if not os.path.exists(self.input_matrix):
            return False

        # 获取文件大小|Get file size
        file_size = Path(self.input_matrix).stat().st_size
        size_gb = file_size / (1024**3)

        # 估算行数（假设每行约2KB）|Estimate line count (assume ~2KB per line)
        estimated_lines = file_size / 2048

        # 判断逻辑|Decision logic:
        # 1. 文件大小 > 10GB
        # 2. 或 估算行数 > threshold_lines
        use_bucket = size_gb > 10 or estimated_lines > threshold_lines

        return use_bucket

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.input_matrix):
            errors.append(f"输入矩阵文件不存在|Input matrix file does not exist: {self.input_matrix}")

        # 检查参数范围|Check parameter ranges
        if self.min_agg_count < 0:
            errors.append(f"最小聚合频次必须为非负整数|Minimum aggregated count must be non-negative: {self.min_agg_count}")

        if self.num_samples == 0:
            errors.append("未能解析样本名|Failed to parse sample names from input file")

        if errors:
            raise ValueError("\n".join(errors))

        return True
