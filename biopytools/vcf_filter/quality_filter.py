"""
质量筛选模块 | Quality Filtering Module
"""

import re
import os
from typing import Dict, Any

class QualityFilter:
    """质量筛选器 | Quality Filter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def should_keep_variant(self, variant: dict) -> bool:
        """判断是否保留变异位点基于质量 | Determine if variant should be kept based on quality"""
        
        # 质量值筛选 | Quality score filtering
        if self.config.quality_threshold is not None:
            qual_str = variant['qual']
            if qual_str == '.' or qual_str == 'PASS':
                pass  # 保留 | Keep
            else:
                try:
                    qual_value = float(qual_str)
                    if qual_value < self.config.quality_threshold:
                        return False
                except ValueError:
                    # 无法解析质量值，跳过 | Cannot parse quality value, skip
                    pass
        
        # 深度筛选 | Depth filtering
        if self.config.min_depth is not None or self.config.max_depth is not None:
            dp_value = self._extract_depth(variant)
            if dp_value is not None:
                if self.config.min_depth is not None and dp_value < self.config.min_depth:
                    return False
                if self.config.max_depth is not None and dp_value > self.config.max_depth:
                    return False
        
        # 双等位基因位点筛选 | Biallelic sites filtering
        if self.config.biallelic_only:
            alt_alleles = variant['alt'].split(',')
            if len(alt_alleles) > 1:
                return False
        
        # 移除插入缺失 | Remove indels
        if self.config.remove_indels:
            ref = variant['ref']
            alt = variant['alt']
            # 简单判断是否为indel | Simple indel detection
            if len(ref) > 1 or any(len(a) > 1 for a in alt.split(',')):
                return False
        
        return True
    
    def _extract_depth(self, variant: dict) -> int:
        """提取深度信息 | Extract depth information"""
        info = variant['info']
        
        # 从INFO字段提取DP | Extract DP from INFO field
        dp_match = re.search(r'DP=(\d+)', info)
        if dp_match:
            return int(dp_match.group(1))
        
        return None
