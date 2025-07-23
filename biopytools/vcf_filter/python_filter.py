"""
Python原生筛选模块 | Python Native Filtering Module
"""

import os
from .position_filter import PositionFilter
from .quality_filter import QualityFilter
from .sample_filter import SampleFilter
from .variant_filter import VariantFilter

class PythonFilter:
    """Python原生筛选器 | Python Native Filter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        
        # 初始化各个筛选器 | Initialize filters
        self.position_filter = PositionFilter(config, logger)
        self.quality_filter = QualityFilter(config, logger)
        self.sample_filter = SampleFilter(config, logger)
        self.variant_filter = VariantFilter(config, logger)
    
    def filter_vcf(self, vcf_reader) -> str:
        """使用Python进行VCF筛选 | Filter VCF using Python"""
        
        if self.config.verbose:
            self.logger.info("使用Python进行VCF筛选 | Filtering VCF using Python")
        
        # 读取头部信息 | Read header
        header_lines, samples = vcf_reader.read_vcf_header()
        
        # 筛选样本 | Filter samples
        filtered_header, filtered_samples, keep_indices = \
            self.sample_filter.filter_samples(header_lines, samples)
        
        # 统计变异位点 | Count variants
        total_variants = 0
        kept_variants = 0
        
        # 写入筛选后的VCF文件 | Write filtered VCF file
        with open(self.config.output_file, 'w') as out_f:
            # 写入头部 | Write header
            for line in filtered_header:
                out_f.write(line + '\n')
            
            # 处理变异位点 | Process variants
            for variant in vcf_reader.iterate_variants():
                total_variants += 1
                
                # 应用各种筛选条件 | Apply various filtering conditions
                if not self.position_filter.should_keep_variant(variant):
                    continue
                
                if not self.quality_filter.should_keep_variant(variant):
                    continue
                
                if not self.variant_filter.should_keep_variant(variant):
                    continue
                
                # 筛选样本数据 | Filter sample data
                if keep_indices:
                    variant = self.sample_filter.filter_variant_samples(variant, keep_indices)
                
                # 写入变异位点 | Write variant
                out_f.write(variant['raw_line'] + '\n')
                kept_variants += 1
                
                if self.config.verbose and total_variants % 10000 == 0:
                    self.logger.info(f"已处理 {total_variants} 个变异位点 | "
                                   f"Processed {total_variants} variants")
        
        if self.config.verbose:
            self.logger.info(f"筛选完成: {total_variants} -> {kept_variants} 个变异位点 | "
                           f"Filtering completed: {total_variants} -> {kept_variants} variants")
        
        return self.config.output_file
