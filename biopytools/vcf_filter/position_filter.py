"""
位置筛选模块 | Position Filtering Module
"""

import os
from typing import Union, List

class PositionFilter:
    """位置筛选器 | Position Filter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def should_keep_variant(self, variant: dict) -> bool:
        """判断是否保留变异位点 | Determine if variant should be kept"""
        
        # 染色体筛选 | Chromosome filtering
        if self.config.chr_name is not None:
            if isinstance(self.config.chr_name, list):
                if variant['chrom'] not in self.config.chr_name:
                    return False
            else:
                if variant['chrom'] != str(self.config.chr_name):
                    return False
        
        # 位置筛选 | Position filtering
        pos = variant['pos']
        
        if self.config.start is not None and pos < self.config.start:
            return False
        
        if self.config.end is not None and pos > self.config.end:
            return False
        
        return True
