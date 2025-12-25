"""
📊 VCF基因型统计数据处理模块 | VCF Genotype Statistics Data Processing Module
"""

import os
from collections import defaultdict, Counter
from .utils import VCFReader, GenotypeParser

class VCFProcessor:
    """🔄 VCF处理器 | VCF Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.vcf_reader = VCFReader(config.vcf_file, logger)
        self.genotype_parser = GenotypeParser()
        
        # 📋 统计数据结构 | Statistics data structures
        self.sample_stats = defaultdict(lambda: defaultdict(int))
        self.detailed_stats = defaultdict(lambda: defaultdict(int))
        self.total_snps = 0
        self.processed_snps = 0
    
    def process_vcf(self):
        """🚀 处理VCF文件 | Process VCF file"""
        self.logger.info("🚀 开始处理VCF文件 | Starting VCF file processing")
        
        # 📋 首先解析头部获取样本名称 | First parse header to get sample names
        sample_names = self.vcf_reader.parse_header()
        if not sample_names:
            self.logger.error("❌ 未在VCF文件中找到样本信息 | No sample information found in VCF file")
            return False
        
        self.config.sample_names = sample_names
        
        # 🔄 处理VCF数据行 | Process VCF data lines
        with self.vcf_reader.open_file() as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                # 🔍 解析数据行 | Parse data line
                if self._process_variant_line(line.strip(), sample_names):
                    self.processed_snps += 1
                
                self.total_snps += 1
                
                # 📊 进度报告 | Progress report
                if self.total_snps % 10000 == 0:
                    self.logger.info(f"📊 已处理 {self.total_snps} 个变异位点 | Processed {self.total_snps} variants")
        
        self.config.total_snps = self.total_snps
        self.logger.info(f"✅ VCF处理完成 | VCF processing completed")
        self.logger.info(f"📊 总变异位点数 | Total variants: {self.total_snps}")
        self.logger.info(f"🎯 通过过滤的位点数 | Variants passed filters: {self.processed_snps}")
        
        return True
    
    def _process_variant_line(self, line: str, sample_names: list) -> bool:
        """🔍 处理单个变异位点行 | Process single variant line"""
        fields = line.split('\t')
        
        if len(fields) < 9 + len(sample_names):
            return False
        
        # 📋 提取变异信息 | Extract variant information
        chrom, pos, var_id, ref, alt, qual, filter_col, info = fields[:8]
        format_field = fields[8]
        genotype_fields = fields[9:9+len(sample_names)]
        
        # 📏 应用质量过滤 | Apply quality filtering
        if self.config.min_qual > 0:
            try:
                if float(qual) < self.config.min_qual:
                    return False
            except (ValueError, TypeError):
                if qual == '.':
                    return False
        
        # 🧬 处理每个样本的基因型 | Process genotype for each sample
        for i, (sample_name, gt_field) in enumerate(zip(sample_names, genotype_fields)):
            self._process_sample_genotype(sample_name, gt_field, format_field)
        
        return True
    
    def _process_sample_genotype(self, sample_name: str, gt_field: str, format_field: str):
        """🧬 处理单个样本的基因型 | Process genotype for single sample"""
        if not gt_field or gt_field == '.':
            if not self.config.exclude_missing:
                self.sample_stats[sample_name]['missing'] += 1
                self.detailed_stats[sample_name]['./.'] += 1
            return
        
        # 🔍 解析格式字段 | Parse format field
        format_keys = format_field.split(':')
        gt_values = gt_field.split(':')
        
        if len(gt_values) != len(format_keys):
            if not self.config.exclude_missing:
                self.sample_stats[sample_name]['missing'] += 1
                self.detailed_stats[sample_name]['./.'] += 1
            return
        
        # 🧬 提取基因型 | Extract genotype
        if 'GT' not in format_keys:
            return
        
        gt_index = format_keys.index('GT')
        gt_str = gt_values[gt_index]
        
        # 📏 应用深度过滤 | Apply depth filtering
        if self.config.min_depth > 0 and 'DP' in format_keys:
            dp_index = format_keys.index('DP')
            try:
                depth = int(gt_values[dp_index])
                if depth < self.config.min_depth:
                    return
            except (ValueError, IndexError):
                pass
        
        # 🔍 解析并统计基因型 | Parse and count genotype
        gt_type = self.genotype_parser.parse_genotype(gt_str)
        gt_category = self.genotype_parser.get_genotype_category(gt_str)
        
        if gt_type == 'missing' and self.config.exclude_missing:
            return
        
        # 📊 更新统计 | Update statistics
        self.sample_stats[sample_name][gt_type] += 1
        self.detailed_stats[sample_name][gt_category] += 1

class StatisticsCalculator:
    """📈 统计计算器 | Statistics Calculator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def calculate_rates(self, sample_stats: dict) -> dict:
        """📊 计算各种比率 | Calculate various rates"""
        results = {}
        
        for sample_name, stats in sample_stats.items():
            total_calls = sum(stats.values())
            total_valid = total_calls - stats.get('missing', 0)
            
            if total_calls == 0:
                continue
            
            # 📋 计算基本统计 | Calculate basic statistics
            sample_result = {
                'total_sites': total_calls,
                'valid_calls': total_valid,
                'missing_calls': stats.get('missing', 0),
                'hom_ref_count': stats.get('hom_ref', 0),
                'het_count': stats.get('het', 0),
                'hom_alt_count': stats.get('hom_alt', 0),
            }
            
            # 📈 计算比率 | Calculate rates
            if total_calls > 0:
                sample_result.update({
                    'missing_rate': stats.get('missing', 0) / total_calls,
                    'call_rate': total_valid / total_calls,
                })
            
            if total_valid > 0:
                sample_result.update({
                    'hom_ref_rate': stats.get('hom_ref', 0) / total_valid,
                    'het_rate': stats.get('het', 0) / total_valid,
                    'hom_alt_rate': stats.get('hom_alt', 0) / total_valid,
                    'heterozygosity_rate': stats.get('het', 0) / total_valid,
                    'homozygosity_rate': (stats.get('hom_ref', 0) + stats.get('hom_alt', 0)) / total_valid,
                })
            else:
                sample_result.update({
                    'hom_ref_rate': 0.0, 'het_rate': 0.0, 'hom_alt_rate': 0.0,
                    'heterozygosity_rate': 0.0, 'homozygosity_rate': 0.0,
                })
            
            results[sample_name] = sample_result
        
        return results