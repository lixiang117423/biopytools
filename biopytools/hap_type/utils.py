"""
单倍型可视化工具函数模块|Haplotype Visualization Utility Functions Module
"""

import logging
import sys
from pathlib import Path
from typing import List, Dict, Any
from collections import Counter


class HapTypeLogger:
    """单倍型可视化日志管理器|Haplotype Visualization Logger Manager"""

    def __init__(self, log_file: Path, log_name: str = "hap_type.log"):
        self.log_file = log_file
        self.log_name = log_name
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


def parse_gff_attributes(attributes_str: str) -> Dict[str, str]:
    """解析GFF属性字符串|Parse GFF attribute string

    Args:
        attributes_str: GFF属性字段|GFF attribute field (e.g., "ID=gene001;Name=ABC")

    Returns:
        dict: 属性键值对|Attribute key-value pairs
    """
    attributes = {}
    for item in attributes_str.split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            attributes[key.strip()] = value.strip()
    return attributes


def format_genotype(gt: List[int]) -> str:
    """格式化基因型为字符串|Format genotype as string

    Args:
        gt: 基因型列表|Genotype list (e.g., [0, 1], [1, 1])

    Returns:
        str: 格式化的基因型|Formatted genotype (e.g., "0/1", "1/1")
    """
    if len(gt) == 2:
        return f"{gt[0]}/{gt[1]}"
    elif len(gt) == 1:
        return str(gt[0])
    else:
        return "./."


def count_haplotypes(haplotypes: List[List[int]]) -> Dict[tuple, int]:
    """统计单倍型频率|Count haplotype frequencies

    Args:
        haplotypes: 单倍型列表|Haplotype list (每个单倍型是SNP索引列表|Each haplotype is a list of SNP indices)

    Returns:
        dict: 单倍型频率字典|Haplotype frequency dictionary
    """
    # 将每个单倍型转换为元组作为字典的键|Convert each haplotype to tuple as dictionary key
    haplotype_tuples = [tuple(h) for h in haplotypes]
    return dict(Counter(haplotype_tuples))


def get_top_haplotypes(haplotype_counts: Dict[tuple, int], top: int) -> List[tuple]:
    """获取频率最高的top N单倍型|Get top N most frequent haplotypes

    Args:
        haplotype_counts: 单倍型频率字典|Haplotype frequency dictionary
        top: 返回前N个|Return top N

    Returns:
        list: 排序后的单倍型元组列表|Sorted list of haplotype tuples
    """
    # 按频率降序排序|Sort by frequency descending
    sorted_haplotypes = sorted(
        haplotype_counts.items(),
        key=lambda x: x[1],
        reverse=True
    )
    # 返回前top个单倍型的元组|Return top N haplotype tuples
    return [h[0] for h in sorted_haplotypes[:top]]


def calculate_maf(ref_count: int, alt_count: int) -> float:
    """计算最小等位基因频率|Calculate minor allele frequency

    Args:
        ref_count: 参考等位基因计数|Reference allele count
        alt_count: 变异等位基因计数|Alternative allele count

    Returns:
        float: MAF值|MAF value
    """
    total = ref_count + alt_count
    if total == 0:
        return 0.0
    return min(ref_count, alt_count) / total


def calculate_missing_rate(genotypes: List[List[int]]) -> float:
    """计算缺失率|Calculate missing rate

    Args:
        genotypes: 基因型列表|Genotype list

    Returns:
        float: 缺失率|Missing rate
    """
    if not genotypes:
        return 0.0

    missing_count = sum(1 for gt in genotypes if -1 in gt or gt.count(-1) > 0)
    return missing_count / len(genotypes)


def merge_intervals(intervals: List[tuple]) -> List[tuple]:
    """合并重叠的区间|Merge overlapping intervals

    Args:
        intervals: 区间列表|Interval list [(start1, end1), (start2, end2), ...]

    Returns:
        list: 合并后的区间列表|Merged interval list
    """
    if not intervals:
        return []

    # 按起始位置排序|Sort by start position
    sorted_intervals = sorted(intervals, key=lambda x: x[0])

    merged = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        last = merged[-1]
        # 如果当前区间与上一个区间重叠或相连|If current interval overlaps or connects with last
        if current[0] <= last[1]:
            # 合并区间|Merge intervals
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)

    return merged


def validate_interval(start: int, end: int, gene_start: int, gene_end: int) -> bool:
    """验证区间是否与基因区间重叠|Validate if interval overlaps with gene interval

    Args:
        start: 区间起始|Interval start
        end: 区间终止|Interval end
        gene_start: 基因起始|Gene start
        gene_end: 基因终止|Gene end

    Returns:
        bool: 是否重叠|Whether overlaps
    """
    return not (end < gene_start or start > gene_end)


def parse_interval(interval_str: str) -> tuple:
    """解析区间字符串|Parse interval string

    Args:
        interval_str: 区间字符串|Interval string (chr:start-end)

    Returns:
        tuple: (chrom, start, end)

    Raises:
        ValueError: 区间格式不正确时|When interval format is incorrect
    """
    import re

    interval_str = interval_str.strip()

    # 匹配 chr:start-end|Match chr:start-end
    match = re.match(r'^([^:]+):\s*(\d+)\s*[-]\s*(\d+)$', interval_str)
    if match:
        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        return chrom, start, end

    raise ValueError(
        f"区间格式不正确，应为 chr:start-end|"
        f"Incorrect interval format, should be chr:start-end: {interval_str}"
    )


def safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    """安全除法，避免除以零|Safe division to avoid division by zero

    Args:
        numerator: 分子|Numerator
        denominator: 分母|Denominator
        default: 除数为零时的默认值|Default value when denominator is zero

    Returns:
        float: 除法结果|Division result
    """
    if denominator == 0:
        return default
    return numerator / denominator


def format_number(num: float, precision: int = 2) -> str:
    """格式化数字|Format number

    Args:
        num: 数字|Number
        precision: 小数位数|Decimal places

    Returns:
        str: 格式化后的字符串|Formatted string
    """
    if num >= 1000000:
        return f"{num/1000000:.{precision}f}M"
    elif num >= 1000:
        return f"{num/1000:.{precision}f}K"
    else:
        return f"{num:.{precision}f}"
