"""
FASTA ID分割工具函数模块 | FASTA ID Splitting Utility Functions Module
"""

import logging
import re
import sys
from pathlib import Path
from typing import List, Tuple, Optional

class SplitLogger:
    """分割分析日志管理器 | Splitting Analysis Logger Manager"""
    
    def __init__(self, output_file: str, log_name: str = "fasta_split.log"):
        self.output_dir = Path(output_file).parent
        self.log_file = self.output_dir / log_name
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

class FastaParser:
    """FASTA文件解析器 | FASTA File Parser"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def detect_delimiter(self, header_lines: List[str]) -> str:
        """自动检测分隔符类型 | Automatically detect delimiter type"""
        space_count = 0
        tab_count = 0
        
        for line in header_lines[:10]:  # 检查前10个序列名称行
            if ' ' in line:
                space_count += 1
            if '\t' in line:
                tab_count += 1
        
        if tab_count > space_count:
            self.logger.info("🔍 检测到制表符分隔 | Detected tab delimiter")
            return "tab"
        elif space_count > 0:
            self.logger.info("🔍 检测到空格分隔 | Detected space delimiter")
            return "space"
        else:
            self.logger.info("🔍 未检测到分隔符，使用空格作为默认 | No delimiter detected, using space as default")
            return "space"
    
    def split_header(self, header: str, delimiter: str, position: int) -> str:
        """分割序列名称行 | Split sequence header line"""
        # 移除开头的 > 符号 | Remove leading > symbol
        header = header.lstrip('>')
        
        # 根据分隔符类型进行分割 | Split based on delimiter type
        if delimiter == "auto":
            # 自动检测：先尝试制表符，再尝试空格 | Auto detect: try tab first, then space
            if '\t' in header:
                parts = header.split('\t')
            else:
                parts = re.split(r'\s+', header)
        elif delimiter == "tab":
            parts = header.split('\t')
        elif delimiter == "space":
            parts = re.split(r'\s+', header)
        elif delimiter == "both":
            # 同时按制表符和空格分割 | Split by both tab and space
            parts = re.split(r'[\s\t]+', header)
        else:
            # 使用用户指定的字符分隔符 | Use user-specified character delimiter
            # 转义特殊正则表达式字符 | Escape special regex characters
            escaped_delimiter = re.escape(delimiter)
            parts = re.split(escaped_delimiter, header)
        
        # 去除空字符串 | Remove empty strings
        parts = [part.strip() for part in parts if part.strip()]
        
        # 返回指定位置的元素 | Return element at specified position
        if position < len(parts):
            return parts[position]
        else:
            # 如果位置超出范围，返回第一个元素 | If position out of range, return first element
            return parts[0] if parts else header.strip()
    
    def get_fasta_stats(self, file_path: str) -> dict:
        """获取FASTA文件统计信息 | Get FASTA file statistics"""
        stats = {
            'total_sequences': 0,
            'header_lines': [],
            'file_size': 0
        }
        
        try:
            file_size = Path(file_path).stat().st_size
            stats['file_size'] = file_size
            
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        stats['total_sequences'] += 1
                        if len(stats['header_lines']) < 10:  # 保存前10个用于分析
                            stats['header_lines'].append(line)
            
        except Exception as e:
            self.logger.error(f"❌ 读取文件统计信息失败 | Failed to read file statistics: {e}")
        
        return stats
