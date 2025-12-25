"""
序列提取工具函数模块 | Sequence Extraction Utility Functions Module
"""

import logging
import sys
from pathlib import Path
from typing import Dict, List

class SequenceLogger:
    """序列提取日志管理器 | Sequence Extraction Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "sequence_extraction.log"):
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

class SequenceValidator:
    """序列验证器 | Sequence Validator"""
    
    def __init__(self, logger):
        self.logger = logger
        self.valid_bases = set('ATCGN-')
    
    def validate_sequence(self, sequence: str, sample_name: str = "") -> bool:
        """验证序列是否有效 | Validate if sequence is valid"""
        if not sequence:
            self.logger.warning(f"空序列 | Empty sequence for sample: {sample_name}")
            return False
        
        invalid_bases = set(sequence.upper()) - self.valid_bases
        if invalid_bases:
            self.logger.warning(f"序列包含无效碱基 | Sequence contains invalid bases: {invalid_bases} in sample: {sample_name}")
            return False
        
        return True
    
    def check_sequence_length(self, sequences: Dict[str, str], expected_length: int) -> bool:
        """检查序列长度一致性 | Check sequence length consistency"""
        inconsistent_samples = []
        
        for sample_name, sequence in sequences.items():
            if len(sequence) != expected_length:
                inconsistent_samples.append(f"{sample_name}: {len(sequence)}")
        
        if inconsistent_samples:
            self.logger.error(f"序列长度不一致 | Inconsistent sequence lengths (期望{expected_length}bp | expected {expected_length}bp):")
            for sample_info in inconsistent_samples:
                self.logger.error(f"  - {sample_info}")
            return False
        
        return True

class SequenceFormatter:
    """序列格式化器 | Sequence Formatter"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def format_sequence_name(self, sample_name: str, chrom: str, start: int, end: int) -> str:
        """格式化序列名称 | Format sequence name"""
        return f"{sample_name}_{chrom}:{start}-{end}"
    
    def wrap_sequence(self, sequence: str, width: int = 80) -> str:
        """序列换行 | Wrap sequence"""
        return '\n'.join(sequence[i:i+width] for i in range(0, len(sequence), width))
    
    def calculate_sequence_stats(self, sequences: Dict[str, str]) -> Dict[str, Dict[str, int]]:
        """计算序列统计信息 | Calculate sequence statistics"""
        stats = {}
        
        for sample_name, sequence in sequences.items():
            base_counts = {
                'A': sequence.count('A'),
                'T': sequence.count('T'),
                'C': sequence.count('C'),
                'G': sequence.count('G'),
                'N': sequence.count('N'),
                '-': sequence.count('-'),
                'total': len(sequence)
            }
            
            # 计算GC含量 | Calculate GC content
            gc_count = base_counts['G'] + base_counts['C']
            total_valid = base_counts['total'] - base_counts['N'] - base_counts['-']
            gc_percent = (gc_count / total_valid * 100) if total_valid > 0 else 0
            
            base_counts['GC_percent'] = round(gc_percent, 2)
            stats[sample_name] = base_counts
        
        return stats

def check_dependencies():
    """检查依赖 | Check dependencies"""
    try:
        import pysam
    except ImportError:
        print("错误 | Error: 需要安装 pysam 库 | pysam library required")
        print("请运行 | Please run: pip install pysam")
        sys.exit(1)
