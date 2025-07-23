"""
多转录组Augustus基因预测工具包 | Multiple RNA-seq Augustus Gene Prediction Toolkit
功能: 结合多个转录组数据进行Augustus基因预测的完整流程 | 
Features: Complete pipeline for Augustus gene prediction with multiple RNA-seq datasets
作者 | Author: Xiang LI  
版本 | Version: v1.0 - 修复版本 | Fixed version
日期 | Date: 2025-07-17

使用示例 | Usage Examples:
    # 方法1: 自动发现样本 (推荐) | Method 1: Auto-discover samples (Recommended)
    python run_augustus_multi_rnaseq.py -g genome.fasta -i /path/to/fastq -s model_name
    
    # 方法2: 使用配置文件 | Method 2: Using config file
    python run_augustus_multi_rnaseq.py -g genome.fasta -c samples.txt -s model_name
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import AugustusMultiRNASeq
from .config import AugustusConfig

__all__ = ['AugustusMultiRNASeq', 'AugustusConfig']
