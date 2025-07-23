"""
Augustus基因预测流水线工具包 | Augustus Gene Prediction Pipeline Toolkit
功能: 完整的Augustus基因预测训练和评估流水线 | 
Features: Complete Augustus gene prediction training and evaluation pipeline
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-23

使用示例 | Usage Examples:
    from biopytools.augustus_pipeline import AugustusPipeline, PipelineConfig
    
    # 创建流水线 | Create pipeline
    pipeline = AugustusPipeline(
        species_name="Rice_NLR_Model",
        genome_file="genome.fa",
        gff_file="annotations.gff3",
        output_dir="augustus_results"
    )
    
    # 运行完整流水线 | Run complete pipeline
    pipeline.run_complete_pipeline()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import AugustusPipeline
from .config import PipelineConfig

__all__ = ['AugustusPipeline', 'PipelineConfig']
