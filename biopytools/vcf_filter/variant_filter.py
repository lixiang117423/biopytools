"""
变异位点筛选模块 | Variant Site Filtering Module
"""

class VariantFilter:
    """变异位点筛选器 | Variant Site Filter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def should_keep_variant(self, variant: dict) -> bool:
        """判断是否保留变异位点 | Determine if variant should be kept"""
        
        # ID筛选 | ID filtering
        if self.config.keep_ids is not None:
            if variant['id'] not in self.config.keep_ids:
                return False
        
        if self.config.remove_ids is not None:
            if variant['id'] in self.config.remove_ids:
                return False
        
        return True
