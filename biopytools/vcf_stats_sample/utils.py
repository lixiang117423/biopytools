"""
VCF基因型统计工具函数模块 | VCF Genotype Statistics Utility Functions Module
"""

import logging
import sys
from pathlib import Path
import gzip

class VCFStatsLogger:
    """VCF基因型统计日志管理器 | VCF Genotype Statistics Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "vcf_stats_analysis.log"):
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

class GenotypeParser:
    """基因型解析器 | Genotype Parser
    
    支持的基因型格式示例 | Supported genotype format examples:
    - 未定相 | Unphased: 0/0, 0/1, 1/1, 0/2, 1/2, ./.
    - 已定相 | Phased: 0|0, 0|1, 1|1, 0|2, 1|2, .|.
    - 带FORMAT信息 | With FORMAT info: 0/1:20:10,10:30:99, 0|1:25:12,13:35:99
    
    所有格式都会被标准化处理为统一的统计类别。
    All formats will be normalized to unified statistical categories.
    """
    
    @staticmethod
    def parse_genotype(gt_str: str) -> str:
        """解析基因型字符串 | Parse genotype string
        处理定相(0|1)和未定相(0/1)基因型 | Handle both phased (0|1) and unphased (0/1) genotypes
        """
        if not gt_str or gt_str in ['.', './.', '.']:
            return 'missing'
        
        # 统一处理定相(|)和未定相(/)基因型，移除其他FORMAT字段信息
        # Handle both phased (|) and unphased (/) genotypes, remove other FORMAT field info
        gt_str = gt_str.replace('|', '/').split(':')[0]
        
        # 标准化基因型 | Normalize genotype
        if '/' in gt_str:
            alleles = gt_str.split('/')
            if len(alleles) == 2:
                try:
                    a1, a2 = int(alleles[0]), int(alleles[1])
                    if a1 == a2:
                        if a1 == 0:
                            return 'hom_ref'  # 0/0
                        else:
                            return 'hom_alt'  # 1/1, 2/2, etc.
                    else:
                        return 'het'  # 0/1, 0/2, 1/2, etc.
                except ValueError:
                    return 'missing'
        
        return 'missing'
    
    @staticmethod
    def get_genotype_category(gt_str: str) -> str:
        """获取基因型类别的详细描述 | Get detailed genotype category description
        处理定相(0|1)和未定相(0/1)基因型 | Handle both phased (0|1) and unphased (0/1) genotypes
        """
        if not gt_str or gt_str in ['.', './.', '.']:
            return './.'
        
        # 统一处理定相(|)和未定相(/)基因型，移除其他FORMAT字段信息
        # Handle both phased (|) and unphased (/) genotypes, remove other FORMAT field info
        gt_str = gt_str.replace('|', '/').split(':')[0]
        
        if '/' in gt_str:
            alleles = gt_str.split('/')
            if len(alleles) == 2:
                try:
                    a1, a2 = int(alleles[0]), int(alleles[1])
                    # 标准化表示 (较小的等位基因在前)
                    return f"{min(a1, a2)}/{max(a1, a2)}"
                except ValueError:
                    return './.'
        
        return './.'

class VCFReader:
    """VCF文件读取器 | VCF File Reader"""
    
    def __init__(self, vcf_file: str, logger):
        self.vcf_file = vcf_file
        self.logger = logger
        self.is_gzipped = vcf_file.endswith('.gz')
    
    def open_file(self):
        """打开VCF文件 | Open VCF file"""
        if self.is_gzipped:
            return gzip.open(self.vcf_file, 'rt', encoding='utf-8')
        else:
            return open(self.vcf_file, 'r', encoding='utf-8')
    
    def parse_header(self):
        """解析VCF头部信息 | Parse VCF header"""
        sample_names = []
        
        with self.open_file() as f:
            for line in f:
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    # 解析样本名称 | Parse sample names
                    fields = line.strip().split('\t')
                    if len(fields) > 9:
                        sample_names = fields[9:]  # 从第10列开始是样本名
                    break
        
        self.logger.info(f"检测到 {len(sample_names)} 个样本 | Detected {len(sample_names)} samples")
        if sample_names:
            self.logger.info(f"样本名称 | Sample names: {', '.join(sample_names[:5])}{'...' if len(sample_names) > 5 else ''}")
        
        return sample_names
