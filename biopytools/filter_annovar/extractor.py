"""
变异提取模块 | Variant Extraction Module
"""

import os
from typing import List, Dict
from .parser import VariantParser
from .coordinate import CoordinateConverter

class VariantExtractor:
    """变异提取器 | Variant Extractor"""
    
    def __init__(self, logger):
        self.logger = logger
        self.parser = VariantParser(logger)
        self.converter = CoordinateConverter()
    
    def extract_variants_in_region(self, variant_file: str, chrom: str, 
                                   region_start: int, region_end: int,
                                   gene_start: int, gene_end: int, 
                                   strand: str, is_exonic: bool = False) -> List[Dict]:
        """
        提取指定区域内的变异
        Extract variants in specified region
        """
        self.logger.info(f"📄 正在处理变异文件 | Processing variant file: {os.path.basename(variant_file)}")
        variants = []
        
        if not os.path.exists(variant_file):
            self.logger.warning(f"⚠️  文件不存在 | File not found: {variant_file}")
            return variants
        
        with open(variant_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                
                variant_data = self.parser.parse_variant_line(line, is_exonic)
                if variant_data is None:
                    continue
                
                # 检查是否在目标区域 | Check if in target region
                if variant_data['chrom'] == chrom and region_start <= variant_data['position'] <= region_end:
                    # 计算基因坐标 | Calculate gene coordinate
                    gene_coord, location_type = self.converter.calculate_gene_coordinate(
                        variant_data['position'], gene_start, gene_end, strand
                    )
                    variant_data['gene_coordinate'] = gene_coord
                    variant_data['location_type'] = location_type
                    variants.append(variant_data)
        
        self.logger.info(f"✅ 找到 {len(variants)} 个变异 | Found {len(variants)} variants")
        return variants
