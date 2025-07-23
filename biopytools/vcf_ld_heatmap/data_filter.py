"""
数据过滤模块 | Data Filtering Module
"""

import numpy as np
from typing import Tuple

try:
    import allel
    ALLEL_AVAILABLE = True
except ImportError:
    ALLEL_AVAILABLE = False
    allel = None

class DataFilter:
    """数据过滤器 | Data Filter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        
        # 检查 scikit-allel 是否可用
        if not ALLEL_AVAILABLE:
            error_msg = (
                "错误: 需要安装scikit-allel库 | Error: scikit-allel library required\n"
                "请运行 | Please run: pip install scikit-allel"
            )
            self.logger.error(error_msg)
            raise ImportError(error_msg)
    
    def filter_variants(self, gt, positions: np.ndarray) -> Tuple:
        """过滤变异位点 | Filter variants"""
        n_variants_original = gt.shape[0]
        
        # 计算等位基因频率 | Calculate allele frequencies
        ac = gt.count_alleles()
        af = ac.to_frequencies()
        
        # MAF过滤 | MAF filtering
        maf = np.minimum(af[:, 0], af[:, 1])  # 取较小的等位基因频率 | Take minor allele frequency
        maf_mask = maf >= self.config.maf
        
        # 移除单态位点 | Remove monomorphic sites
        polymorphic_mask = ac.is_segregating()
        
        # 组合过滤条件 | Combine filtering conditions
        combined_mask = maf_mask & polymorphic_mask
        
        # 应用过滤 | Apply filtering
        gt_filtered = gt.compress(combined_mask, axis=0)
        positions_filtered = positions[combined_mask]
        
        # 如果SNP数量超过限制，随机抽样 | Random sampling if SNP count exceeds limit
        if gt_filtered.shape[0] > self.config.max_snps:
            if self.config.verbose:
                self.logger.info(f"SNP数量 ({gt_filtered.shape[0]}) 超过限制 ({self.config.max_snps})，进行随机抽样 | "
                               f"SNP count ({gt_filtered.shape[0]}) exceeds limit ({self.config.max_snps}), random sampling")
            
            indices = np.random.choice(gt_filtered.shape[0], self.config.max_snps, replace=False)
            indices = np.sort(indices)  # 保持位置顺序 | Maintain position order
            gt_filtered = gt_filtered.take(indices, axis=0)
            positions_filtered = positions_filtered[indices]
        
        if self.config.verbose:
            self.logger.info(f"过滤后剩余 {gt_filtered.shape[0]} 个变异位点 (原始: {n_variants_original}) | "
                           f"Remaining {gt_filtered.shape[0]} variants after filtering (original: {n_variants_original})")
        
        return gt_filtered, positions_filtered
