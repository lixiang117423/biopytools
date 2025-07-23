"""
矩阵导出模块 | Matrix Export Module
"""

import pandas as pd
import numpy as np

class MatrixExporter:
    """矩阵导出器 | Matrix Exporter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def save_ld_matrix(self, ld_matrix: np.ndarray, positions: np.ndarray):
        """保存LD矩阵到CSV文件 | Save LD matrix to CSV file"""
        if not self.config.save_matrix:
            return
        
        if self.config.verbose:
            self.logger.info(f"正在保存LD矩阵到 | Saving LD matrix to: {self.config.save_matrix}")
        
        # 创建DataFrame | Create DataFrame
        df = pd.DataFrame(ld_matrix, 
                         index=positions, 
                         columns=positions)
        
        # 保存到CSV | Save to CSV
        df.to_csv(self.config.save_matrix)
        
        if self.config.verbose:
            self.logger.info("LD矩阵保存完成 | LD matrix save completed")
