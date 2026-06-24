"""
QIIME2微生物组多样性分析工具|QIIME2 Microbiome Diversity Analysis Tool
功能: 双端扩增子数据自动识别并完成QIIME2全流程(ASV/OTU、丰度表、分类注释、多样性、抽平)|
Features: Auto-detect paired-end amplicon data and run full QIIME2 pipeline
"""

from .pipeline import Qiime2Pipeline
from .config import Qiime2Config

__version__ = "1.0.0"

__all__ = [
    'Qiime2Pipeline',
    'Qiime2Config'
]
