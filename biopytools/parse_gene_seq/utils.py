"""
基因序列提取工具函数模块 🛠️ | Gene Sequence Extraction Utility Functions Module
"""

import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

class ExtractionLogger:
    """基因提取日志管理器 📝 | Gene Extraction Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "gene_extraction.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 📝 | Setup logging"""
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
        """获取日志器 📝 | Get logger"""
        return self.logger

class SequenceUtils:
    """序列处理工具类 🧬 | Sequence Processing Utilities"""
    
    @staticmethod
    def reverse_complement(seq: str) -> str:
        """计算反向互补序列 🔄 | Calculate reverse complement sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base.upper(), 'N') for base in reversed(seq))
    
    @staticmethod
    def extract_sequence_region(genome_seq: str, start: int, end: int, strand: str) -> str:
        """从基因组中提取序列片段 ✂️ | Extract sequence fragment from genome"""
        # GFF坐标是1-based，Python是0-based
        seq = genome_seq[start-1:end]
        
        # 如果是负链，需要反向互补
        if strand == '-':
            seq = SequenceUtils.reverse_complement(seq)
        
        return seq
    
    @staticmethod
    def format_fasta_sequence(seq: str, line_width: int = 60) -> str:
        """格式化FASTA序列 📄 | Format FASTA sequence"""
        lines = []
        for i in range(0, len(seq), line_width):
            lines.append(seq[i:i+line_width])
        return '\n'.join(lines)

class FileValidator:
    """文件验证器 🔍 | File Validator"""
    
    @staticmethod
    def validate_fasta_file(file_path: str) -> bool:
        """验证FASTA文件格式 🧬 | Validate FASTA file format"""
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                return first_line.startswith('>')
        except Exception:
            return False
    
    @staticmethod
    def validate_gff_file(file_path: str) -> bool:
        """验证GFF文件格式 📄 | Validate GFF file format"""
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    # 检查是否有9列
                    fields = line.split('\t')
                    return len(fields) == 9
            return True
        except Exception:
            return False
