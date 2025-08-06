"""
GFF3工具函数模块 | GFF3 Utility Functions Module
"""

import logging
import sys
from pathlib import Path
from typing import Dict, Any

class GFFLogger:
    """GFF3处理日志管理器 | GFF3 Processing Logger Manager"""
    
    def __init__(self, output_file: str, log_name: str = "gff_processing.log"):
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

class AttributeParser:
    """属性解析器 | Attribute Parser"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """
        解析GFF3第九列的属性字符串 | Parse GFF3 column 9 attributes string
        
        Args:
            attr_string: 属性字符串，如 "ID=gene1;Name=GeneA"
            
        Returns:
            Dict[str, str]: 属性字典，如 {'ID': 'gene1', 'Name': 'GeneA'}
        """
        attributes = {}
        
        if not attr_string or attr_string.strip() == '.':
            return attributes
        
        try:
            for part in attr_string.strip().split(';'):
                if '=' in part:
                    key, value = part.split('=', 1)
                    attributes[key.strip()] = value.strip()
        except Exception as e:
            self.logger.warning(f"解析属性字符串时出错 | Error parsing attributes: {attr_string}, {e}")
        
        return attributes

class GFFValidator:
    """GFF3文件验证器 | GFF3 File Validator"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def validate_gff_line(self, line: str) -> bool:
        """
        验证GFF3行格式 | Validate GFF3 line format
        
        Args:
            line: GFF3行
            
        Returns:
            bool: 是否有效
        """
        if line.startswith('#'):
            return False
        
        parts = line.strip().split('\t')
        if len(parts) != 9:
            return False
        
        # 检查必需字段 | Check required fields
        try:
            # 检查起始和结束位置是否为数字 | Check if start and end positions are numeric
            int(parts[3])  # start
            int(parts[4])  # end
            
            # 检查链方向 | Check strand
            if parts[6] not in ['+', '-', '.', '?']:
                return False
                
        except (ValueError, IndexError):
            return False
        
        return True
    
    def check_file_format(self, file_path: str) -> bool:
        """
        检查文件格式 | Check file format
        
        Args:
            file_path: 文件路径
            
        Returns:
            bool: 是否为有效的GFF3文件
        """
        try:
            with open(file_path, 'r') as f:
                # 检查前几行 | Check first few lines
                for i, line in enumerate(f):
                    if i > 100:  # 只检查前100行 | Only check first 100 lines
                        break
                    
                    if line.startswith('##gff-version'):
                        return True
                    
                    if not line.startswith('#') and line.strip():
                        # 检查数据行格式 | Check data line format
                        if self.validate_gff_line(line):
                            return True
                        else:
                            self.logger.warning(f"可能的格式错误 | Possible format error at line {i+1}: {line.strip()}")
                            return False
            
            return True
            
        except Exception as e:
            self.logger.error(f"文件格式检查失败 | File format check failed: {e}")
            return False
