"""
样本筛选模块 | Sample Filtering Module
"""

import os
from typing import List, Dict

class SampleFilter:
    """样本筛选器 | Sample Filter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def filter_samples(self, header_lines: List[str], samples: List[str]) -> tuple:
        """筛选样本 | Filter samples"""
        
        if (self.config.keep_samples is None and 
            self.config.remove_samples is None):
            return header_lines, samples, list(range(len(samples)))
        
        # 确定要保留的样本索引 | Determine sample indices to keep
        keep_indices = []
        filtered_samples = []
        
        for i, sample in enumerate(samples):
            keep = True
            
            # 如果指定了要保留的样本 | If keep_samples is specified
            if self.config.keep_samples is not None:
                keep = sample in self.config.keep_samples
            
            # 如果指定了要移除的样本 | If remove_samples is specified
            if self.config.remove_samples is not None:
                if sample in self.config.remove_samples:
                    keep = False
            
            if keep:
                keep_indices.append(i)
                filtered_samples.append(sample)
        
        # 更新头部信息 | Update header
        new_header_lines = []
        for line in header_lines:
            if line.startswith('#CHROM'):
                # 重建样本列 | Rebuild sample columns
                fields = line.split('\t')
                new_fields = fields[:9] + filtered_samples
                new_header_lines.append('\t'.join(new_fields))
            else:
                new_header_lines.append(line)
        
        if self.config.verbose:
            self.logger.info(f"样本筛选: {len(samples)} -> {len(filtered_samples)} | "
                           f"Sample filtering: {len(samples)} -> {len(filtered_samples)}")
        
        return new_header_lines, filtered_samples, keep_indices
    
    def filter_variant_samples(self, variant: dict, keep_indices: List[int]) -> dict:
        """筛选变异位点的样本数据 | Filter variant sample data"""
        if not keep_indices or len(variant['samples']) == 0:
            return variant
        
        # 过滤样本基因型数据 | Filter sample genotype data
        filtered_samples = [variant['samples'][i] for i in keep_indices 
                          if i < len(variant['samples'])]
        
        # 重建行数据 | Rebuild line data
        fields = variant['raw_line'].split('\t')
        new_fields = fields[:9] + filtered_samples
        
        variant['samples'] = filtered_samples
        variant['raw_line'] = '\t'.join(new_fields)
        
        return variant
