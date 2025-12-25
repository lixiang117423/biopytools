"""
SNP Index工具函数模块 | SNP Index Utility Functions Module
"""

import gzip
import logging
import sys
import os
from typing import Tuple, List, Optional, TextIO


def setup_logger(name: str, log_file: Optional[str] = None,
                level: int = logging.INFO) -> logging.Logger:
    """
    设置标准化的日志系统 | Setup standardized logging system

    Args:
        name: logger名称 | Logger name
        log_file: 日志文件路径 | Log file path
        level: 日志级别 | Log level

    Returns:
        logging.Logger: 配置好的logger | Configured logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    # 日志格式 | Log format
    formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # stdout handler - INFO级别 | stdout handler - INFO level
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    stdout_handler.addFilter(lambda record: record.levelno <= logging.INFO)
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)

    # stderr handler - WARNING及以上级别 | stderr handler - WARNING and above
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(formatter)
    logger.addHandler(stderr_handler)

    # 文件handler | File handler
    if log_file:
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def open_file(file_path: str, mode: str = 'r') -> TextIO:
    """
    智能打开文件（支持gzip压缩）| Smart file opener (supports gzip compression)

    Args:
        file_path: 文件路径 | File path
        mode: 打开模式 | Open mode

    Returns:
        TextIO: 文件对象 | File object
    """
    if file_path.endswith('.gz'):
        return gzip.open(file_path, mode + 't')  # text mode
    else:
        return open(file_path, mode)


def parse_vcf_line(line: str) -> Tuple[str, int, str, str, List[str]]:
    """
    解析VCF文件的一行 | Parse a line from VCF file

    Args:
        line: VCF行 | VCF line

    Returns:
        tuple: (染色体, 位置, 参考碱基, 替代碱基, 样本信息) | (chrom, pos, ref, alt, samples)
    """
    parts = line.strip().split('\t')
    if len(parts) < 10:
        raise ValueError(f"Invalid VCF line format: {line[:50]}...")

    chrom = parts[0]
    pos = int(parts[1])
    ref = parts[3]
    alt = parts[4]
    samples = parts[9:]

    return chrom, pos, ref, alt, samples


def extract_ad_values(sample_info: str) -> Tuple[int, int]:
    """
    从样品信息中提取AD值（参考和替代等位基因深度）| Extract AD values from sample info

    Args:
        sample_info: 样品信息字段 | Sample info field

    Returns:
        tuple: (参考深度, 替代深度) | (reference depth, alternative depth)
    """
    # 解析GT:AD:DP:GQ:PL格式 | Parse GT:AD:DP:GQ:PL format
    fields = sample_info.split(':')
    if len(fields) < 2:
        return 0, 0

    try:
        ad_field = fields[1]
        # 处理多个ALT的情况，格式通常是ref_depth,alt1_depth,alt2_depth...
        # Handle multiple ALTs, format usually is ref_depth,alt1_depth,alt2_depth...
        ad_values = [int(x) for x in ad_field.split(',') if x != '.']
        ref_depth = ad_values[0] if len(ad_values) > 0 else 0
        alt_depth = sum(ad_values[1:]) if len(ad_values) > 1 else 0
        return ref_depth, alt_depth
    except (ValueError, IndexError):
        return 0, 0


def calculate_snp_index(ref_depth: int, alt_depth: int) -> float:
    """
    计算SNP index（替代等位基因频率）| Calculate SNP index (alternative allele frequency)

    Args:
        ref_depth: 参考等位基因深度 | Reference allele depth
        alt_depth: 替代等位基因深度 | Alternative allele depth

    Returns:
        float: SNP index值 | SNP index value
    """
    total_depth = ref_depth + alt_depth
    if total_depth == 0:
        return 0.0
    return alt_depth / total_depth


def parse_quality_filters(vcf_line: str) -> Tuple[Optional[int], Optional[int]]:
    """
    解析VCF行中的质量信息 | Parse quality information from VCF line

    Args:
        vcf_line: VCF行 | VCF line

    Returns:
        tuple: (质量值, mapping质量) | (quality, mapping quality)
    """
    parts = vcf_line.strip().split('\t')
    if len(parts) < 8:
        return None, None

    try:
        qual = float(parts[5]) if parts[5] != '.' else None
        if qual is not None:
            qual = int(qual)
    except (ValueError, IndexError):
        qual = None

    # 从INFO字段提取MQ | Extract MQ from INFO field
    info_field = parts[7]
    mq = None
    for item in info_field.split(';'):
        if item.startswith('MQ='):
            try:
                mq = int(item.split('=')[1])
            except (ValueError, IndexError):
                pass
            break

    return qual, mq


def get_sample_names_from_vcf(vcf_file: str) -> List[str]:
    """
    从VCF文件中获取样本名称 | Get sample names from VCF file

    Args:
        vcf_file: VCF文件路径 | VCF file path

    Returns:
        list: 样本名称列表 | Sample names list
    """
    with open_file(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                if len(parts) > 9:
                    return parts[9:]
                break
    return []


def validate_vcf_file(vcf_file: str, min_samples: int = 2) -> bool:
    """
    验证VCF文件格式和样本数 | Validate VCF file format and sample count

    Args:
        vcf_file: VCF文件路径 | VCF file path
        min_samples: 最少样本数 | Minimum sample count

    Returns:
        bool: 验证是否通过 | Whether validation passed
    """
    try:
        sample_names = get_sample_names_from_vcf(vcf_file)
        if len(sample_names) < min_samples:
            return False

        # 检查是否有数据行 | Check if there are data lines
        with open_file(vcf_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    return True
        return False
    except Exception:
        return False


def format_number(num: float, decimal_places: int = 4) -> str:
    """
    格式化数字显示 | Format number display

    Args:
        num: 数字 | Number
        decimal_places: 小数位数 | Decimal places

    Returns:
        str: 格式化字符串 | Formatted string
    """
    return f"{num:.{decimal_places}f}"


def format_large_number(num: int) -> str:
    """
    格式化大数字显示（添加千位分隔符）| Format large number display (add thousand separator)

    Args:
        num: 数字 | Number

    Returns:
        str: 格式化字符串 | Formatted string
    """
    return f"{num:,}"


def check_file_exists(file_path: str, description: str = "文件") -> None:
    """
    检查文件是否存在，不存在则抛出异常 | Check if file exists, raise exception if not

    Args:
        file_path: 文件路径 | File path
        description: 文件描述 | File description

    Raises:
        FileNotFoundError: 文件不存在 | File not found
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{description}不存在: {file_path}")


def ensure_dir_exists(dir_path: str) -> None:
    """
    确保目录存在，不存在则创建 | Ensure directory exists, create if not

    Args:
        dir_path: 目录路径 | Directory path
    """
    os.makedirs(dir_path, exist_ok=True)