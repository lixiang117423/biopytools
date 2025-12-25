"""
GFF转换工具函数模块 | GFF Conversion Utility Functions Module
"""

import logging
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple

class GFFLogger:
    """GFF转换日志管理器 | GFF Conversion Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "gff_conversion.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger

class ChromosomeExtractor:
    """染色体编号提取器 | Chromosome Number Extractor"""
    
    @staticmethod
    def extract_chr_number(chr_name: str) -> str:
        """从染色体名称中提取编号并补0 | Extract number from chromosome name and pad with 0"""
        patterns = [
            r'Chr(\d+)',  # Chr1, Chr2, etc.
            r'chr(\d+)',  # chr1, chr2, etc.
            r'OV(\d+)_',  # OV12_RagTag
            r'LG(\d+)',   # LG1, LG2 (Linkage Group)
            r'scaffold_?(\d+)',  # scaffold1, scaffold_1
            r'(\d+)$',    # 以数字结尾
            r'[A-Za-z]+(\d+)',  # 字母+数字
        ]
        
        for pattern in patterns:
            match = re.search(pattern, chr_name)
            if match:
                chr_num = int(match.group(1))
                # 1-9补0，变成01-09
                return f"{chr_num:02d}"
        
        return "01"  # 默认值也补0

class AttributeParser:
    """GFF属性解析器 | GFF Attribute Parser"""
    
    @staticmethod
    def parse_attributes(attr_string: str) -> Dict[str, str]:
        """解析GFF属性字符串 | Parse GFF attribute string"""
        attributes = {}
        
        if not attr_string or attr_string == '.':
            return attributes
        
        # 分割属性
        pairs = attr_string.split(';')
        
        for pair in pairs:
            if '=' in pair:
                key, value = pair.split('=', 1)
                attributes[key.strip()] = value.strip()
        
        return attributes
    
    @staticmethod
    def format_attributes(attributes: Dict[str, str]) -> str:
        """格式化属性为GFF字符串 | Format attributes to GFF string"""
        if not attributes:
            return '.'
        
        pairs = []
        # 确保ID和Name在前面
        for key in ['ID', 'Name', 'Parent']:
            if key in attributes:
                pairs.append(f"{key}={attributes[key]}")
        
        # 添加其他属性
        for key, value in attributes.items():
            if key not in ['ID', 'Name', 'Parent']:
                pairs.append(f"{key}={value}")
        
        return ';'.join(pairs)

class GFFValidator:
    """GFF文件验证器 | GFF File Validator"""
    
    @staticmethod
    def is_valid_gff_line(line: str) -> bool:
        """验证GFF行格式 | Validate GFF line format"""
        if line.startswith('#') or not line.strip():
            return True
        
        fields = line.strip().split('\t')
        return len(fields) == 9
    
    @staticmethod
    def validate_coordinates(start: int, end: int) -> bool:
        """验证坐标有效性 | Validate coordinate validity"""
        return start > 0 and end >= start

def format_gene_id(species_name: str, species_prefix: str, chr_num: str, gene_num: int) -> str:
    """格式化基因ID | Format gene ID"""
    return f"{species_name}_{species_prefix}{chr_num}G{gene_num:06d}"

def get_file_stats(file_path: str) -> Dict[str, int]:
    """获取文件统计信息 | Get file statistics"""
    stats = {
        'total_lines': 0,
        'header_lines': 0,
        'feature_lines': 0,
        'genes': 0,
        'mrnas': 0,
        'exons': 0
    }
    
    with open(file_path, 'r') as f:
        for line in f:
            stats['total_lines'] += 1
            
            if line.startswith('#'):
                stats['header_lines'] += 1
            elif line.strip():
                stats['feature_lines'] += 1
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    feature_type = fields[2]
                    if feature_type == 'gene':
                        stats['genes'] += 1
                    elif feature_type in ['mRNA', 'transcript']:
                        stats['mrnas'] += 1
                    elif feature_type == 'exon':
                        stats['exons'] += 1
    
    return stats
