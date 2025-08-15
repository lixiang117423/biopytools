"""
🧬 基因型数据处理模块 | Genotype Data Processing Module
"""

from typing import List, Dict, Any, Optional
from collections import defaultdict
from .vcf_parser import VCFRecord
from .utils import GenotypeUtils

class GenotypeProcessor:
    """🧬 基因型处理器 | Genotype Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.results = defaultdict(list)  # 🗂️ 按染色体分组的结果 | Results grouped by chromosome
    
    def process_record(self, record: VCFRecord) -> Optional[Dict[str, Any]]:
        """🔄 处理单个VCF记录 | Process single VCF record"""
        # 🔍 检查是否只要双等位位点 | Check if only biallelic sites are needed
        if self.config.biallelic_only and not GenotypeUtils.is_biallelic(record.ref, record.alt):
            return None
        
        # 🏗️ 构建输出行 | Build output row
        row = {
            'CHROM': record.chrom,
            'POS': record.pos,
            'ID': record.id,
            'REF': record.ref,
            'ALT': record.alt,
            'QUAL': record.qual
        }
        
        # ➕ 添加基因型信息 | Add genotype information
        genotype_values = []
        for sample in record.genotypes:
            gt = record.genotypes[sample]
            row[sample] = gt if gt is not None else './.'
            genotype_values.append(gt)
        
        # 📊 计算纯合/杂合比例 | Calculate homozygous/heterozygous ratios
        homo_ratio, hetero_ratio = GenotypeUtils.calculate_genotype_stats(genotype_values)
        row['Homozygous_Ratio'] = f"{homo_ratio:.4f}"
        row['Heterozygous_Ratio'] = f"{hetero_ratio:.4f}"
        
        return row
    
    def add_record(self, row: Dict[str, Any]):
        """➕ 添加处理后的记录 | Add processed record"""
        chrom = row['CHROM']
        self.results[chrom].append(row)
    
    def get_results_by_chromosome(self) -> Dict[str, List[Dict[str, Any]]]:
        """🧬 按染色体获取结果 | Get results by chromosome"""
        return dict(self.results)
    
    def get_all_results(self) -> List[Dict[str, Any]]:
        """📋 获取所有结果（合并所有染色体） | Get all results (merged from all chromosomes)"""
        all_results = []
        for chrom in sorted(self.results.keys()):
            all_results.extend(self.results[chrom])
        return all_results
    
    def get_summary_stats(self) -> Dict[str, Any]:
        """📊 获取汇总统计 | Get summary statistics"""
        total_variants = sum(len(records) for records in self.results.values())
        chromosomes = list(self.results.keys())
        
        return {
            'total_variants': total_variants,
            'chromosomes': sorted(chromosomes),
            'chromosome_counts': {chrom: len(records) for chrom, records in self.results.items()}
        }