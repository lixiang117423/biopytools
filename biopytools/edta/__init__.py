"""
EDTA转座子注释工具包|EDTA TE Annotation Toolkit
功能: 基于EDTA的全基因组转座子鉴定、分类和注释|Features: Whole-genome transposon identification, classification, and annotation based on EDTA

使用示例|Usage Examples:
    from biopytools.edta import EDTARunner, EDTAConfig

    # 创建EDTA配置|Create EDTA config
    config = EDTAConfig(
        genome="plant.fa",
        cds="cds.fa",
        threads=24
    )

    # 运行EDTA注释|Run EDTA annotation
    runner = EDTARunner(config)
    runner.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import main
from .config import EDTAConfig, PanEDTAConfig
from .edta_runner import EDTARunner
from .panedta_runner import PanEDTARunner

__all__ = [
    'main',
    'EDTAConfig',
    'PanEDTAConfig',
    'EDTARunner',
    'PanEDTARunner'
]
