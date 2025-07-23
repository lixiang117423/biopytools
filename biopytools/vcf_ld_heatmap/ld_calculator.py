"""
连锁不平衡计算模块 | Linkage Disequilibrium Calculation Module
"""

import numpy as np
from scipy.stats import pearsonr

try:
    import allel
    ALLEL_AVAILABLE = True
except ImportError:
    ALLEL_AVAILABLE = False
    allel = None

class LDCalculator:
    """连锁不平衡计算器 | Linkage Disequilibrium Calculator"""
    
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
    
    def calculate_ld_matrix(self, gt) -> np.ndarray:
        """计算连锁不平衡矩阵 | Calculate linkage disequilibrium matrix"""
        if self.config.verbose:
            self.logger.info("正在计算连锁不平衡矩阵 | Calculating linkage disequilibrium matrix...")
        
        # 转换为数值编码 (计算次要等位基因数量) | Convert to numeric encoding (count minor alleles)
        gt_numeric = gt.to_n_alt()
        
        n_variants = gt_numeric.shape[0]
        ld_matrix = np.full((n_variants, n_variants), np.nan)
        
        # 计算r²值 | Calculate r² values
        for i in range(n_variants):
            for j in range(i, n_variants):
                if i == j:
                    ld_matrix[i, j] = 1.0
                else:
                    # 移除缺失值 | Remove missing values
                    mask = ~(np.isnan(gt_numeric[i]) | np.isnan(gt_numeric[j]))
                    if np.sum(mask) > 10:  # 至少需要10个有效样本 | At least 10 valid samples needed
                        try:
                            r, _ = pearsonr(gt_numeric[i][mask], gt_numeric[j][mask])
                            ld_matrix[i, j] = r**2 if not np.isnan(r) else 0.0
                            ld_matrix[j, i] = ld_matrix[i, j]  # 对称矩阵 | Symmetric matrix
                        except:
                            ld_matrix[i, j] = 0.0
                            ld_matrix[j, i] = 0.0
                    else:
                        ld_matrix[i, j] = 0.0
                        ld_matrix[j, i] = 0.0
            
            if self.config.verbose and (i + 1) % 50 == 0:
                self.logger.info(f"已完成 {i + 1}/{n_variants} 个位点的LD计算 | "
                               f"Completed LD calculation for {i + 1}/{n_variants} loci")
        
        return ld_matrix
