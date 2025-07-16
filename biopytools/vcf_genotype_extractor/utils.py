"""
VCF基因型提取工具函数模块 | VCF Genotype Extraction Utility Functions Module
"""

import gzip
import logging
import sys
from pathlib import Path
from typing import Optional

class VCFLogger:
    """VCF提取日志管理器 | VCF Extraction Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "vcf_extraction.log"):
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

class FileUtils:
    """文件工具类 | File Utility Class"""
    
    @staticmethod
    def is_gzipped(file_path: str) -> bool:
        """检查文件是否为gzip压缩 | Check if file is gzip compressed"""
        try:
            with gzip.open(file_path, 'rt') as f:
                f.readline()
            return True
        except:
            return False
    
    @staticmethod
    def open_file(file_path: str, mode: str = 'r'):
        """智能打开文件（自动检测是否压缩） | Smart file opening (auto-detect compression)"""
        if FileUtils.is_gzipped(file_path):
            return gzip.open(file_path, mode + 't', encoding='utf-8')
        else:
            return open(file_path, mode, encoding='utf-8')

class GenotypeUtils:
    """基因型工具类 | Genotype Utility Class"""
    
    @staticmethod
    def parse_genotype(gt_field: str) -> Optional[str]:
        """解析基因型字段 | Parse genotype field"""
        if gt_field in ['.', './.' , '.|.']:
            return None
        
        # 移除相位信息，只保留基因型 | Remove phasing info, keep only genotype
        gt = gt_field.split(':')[0] if ':' in gt_field else gt_field
        
        return gt
    
    @staticmethod
    def is_biallelic(ref: str, alt: str) -> bool:
        """检查是否为双等位变异 | Check if variant is biallelic"""
        # 如果ALT字段不包含逗号，则为双等位 | If ALT field doesn't contain comma, it's biallelic
        return ',' not in alt
    
    @staticmethod
    def calculate_genotype_stats(genotypes: list[Optional[str]]) -> tuple:
        """计算基因型统计信息 | Calculate genotype statistics"""
        valid_gts = [gt for gt in genotypes if gt is not None]
        
        if not valid_gts:
            return 0.0, 0.0
        
        homozygous_count = 0
        heterozygous_count = 0
        
        for gt in valid_gts:
            # 处理不同的分隔符 | Handle different separators
            if '/' in gt:
                alleles = gt.split('/')
            elif '|' in gt:
                alleles = gt.split('|')
            else:
                continue
            
            if len(alleles) == 2:
                if alleles[0] == alleles[1]:
                    homozygous_count += 1
                else:
                    heterozygous_count += 1
        
        total = len(valid_gts)
        homozygous_ratio = homozygous_count / total if total > 0 else 0.0
        heterozygous_ratio = heterozygous_count / total if total > 0 else 0.0
        
        return homozygous_ratio, heterozygous_ratio

def check_dependencies():
    """检查依赖 | Check dependencies"""
    optional_deps = {}
    
    try:
        import cyvcf2
        optional_deps['cyvcf2'] = True
    except ImportError:
        optional_deps['cyvcf2'] = False
    
    try:
        import pandas
        optional_deps['pandas'] = True
    except ImportError:
        optional_deps['pandas'] = False
    
    return optional_deps
