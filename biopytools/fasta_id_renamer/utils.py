"""
FASTA序列ID重命名工具函数模块|FASTA ID Renamer Utility Functions Module
"""

import logging
import sys
from pathlib import Path

class RenamerLogger:
    """序列重命名日志管理器|Sequence Renamer Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
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
            handlers.append(logging.FileHandler(self.log_file, encoding='utf-8'))

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


def format_number_with_padding(num: int, width: int) -> str:
    """格式化数字，添加零填充|Format number with zero padding

    Args:
        num: 数字|Number
        width: 填充宽度|Padding width

    Returns:
        格式化后的字符串|Formatted string

    Examples:
        >>> format_number_with_padding(1, 2)
        '01'
        >>> format_number_with_padding(10, 2)
        '10'
        >>> format_number_with_padding(1, 3)
        '001'
    """
    return str(num).zfill(width)


def parse_sequence_header(header: str):
    """解析FASTA序列头信息|Parse FASTA sequence header

    Args:
        header: FASTA序列头（不包含>符号）|FASTA header (without > symbol)

    Returns:
        tuple: (序列ID, 描述信息)| (sequence ID, description)

    Examples:
        >>> parse_sequence_header("CM081539.1 Erysimum cheiranthoides chromosome 1")
        ('CM081539.1', 'Erysimum cheiranthoides chromosome 1')
        >>> parse_sequence_header("scaffold_1")
        ('scaffold_1', '')
    """
    parts = header.split(None, 1)
    seq_id = parts[0]
    description = parts[1] if len(parts) > 1 else ''
    return seq_id, description


def detect_sequence_type(header: str, description: str = "") -> str:
    """检测序列类型|Detect sequence type

    Args:
        header: 序列ID|Sequence ID
        description: 描述信息|Description

    Returns:
        str: 序列类型|Sequence type
            - 'chromosome': 染色体
            - 'scaffold': scaffold
            - 'contig': contig
            - 'plastid': 叶绿体
            - 'mitochondria': 线粒体
            - 'unknown': 未知类型
    """
    header_lower = header.lower()
    desc_lower = description.lower()

    # 优先检测描述中的明确类型关键词|Prioritize explicit type keywords in description
    # 检测染色体|Detect chromosome
    if 'chromosome' in desc_lower:
        return 'chromosome'

    # 检测叶绿体|Detect chloroplast (在描述中)
    if any(keyword in desc_lower
           for keyword in ['chloroplast', 'plastid', 'plastome']):
        return 'plastid'

    # 检测线粒体|Detect mitochondria (在描述中)
    if any(keyword in desc_lower
           for keyword in ['mitochondr', 'mitogenome', 'mitochondria']):
        return 'mitochondria'

    # 对于ID中包含关键词的情况，需要更精确的匹配|For keywords in ID, need more precise matching
    # 检测scaffold|Detect scaffold
    if any(keyword in header_lower or keyword in desc_lower
           for keyword in ['scaffold', 'unscaffold', 'scaf_', 'scaff', 'scf']):
        return 'scaffold'

    # 检测contig|Detect contig
    if any(keyword in header_lower or keyword in desc_lower
           for keyword in ['contig', 'ctg', 'ctg_', 'cont_']):
        return 'contig'

    # 检测scaffold|Detect scaffold
    if any(keyword in header_lower or keyword in desc_lower
           for keyword in ['scaffold', 'unscaffold', 'scaf_', 'scaff', 'scf']):
        return 'scaffold'

    # 检测contig|Detect contig
    if any(keyword in header_lower for keyword in ['contig', 'ctg', 'ctg_', 'cont_']):
        return 'contig'

    # 尝试从NCBI accession判断|Try to determine from NCBI accession
    # NCBI accession格式: XXXXXXXXXXXX_NNN (通常包含下划线的是未挂载序列)
    if '_' in header and any(c.isdigit() for c in header):
        # 可能是未挂载的scaffold或contig
        return 'scaffold'

    return 'unknown'


def extract_chromosome_number(description: str) -> int:
    """从描述信息中提取染色体编号|Extract chromosome number from description

    Args:
        description: 描述信息|Description

    Returns:
        int: 染色体编号|Chromosome number, 如果未找到返回0|returns 0 if not found

    Examples:
        >>> extract_chromosome_number("chromosome 1")
        1
        >>> extract_chromosome_number("chromosome X")
        0
    """
    import re

    # 查找 "chromosome 数字" 模式|Find "chromosome number" pattern
    match = re.search(r'chromosome\s+(\d+)', description, re.IGNORECASE)
    if match:
        return int(match.group(1))

    # 查找 "chr 数字" 模式|Find "chr number" pattern
    match = re.search(r'chr\s*\(?(\d+)\)?', description, re.IGNORECASE)
    if match:
        return int(match.group(1))

    return 0
