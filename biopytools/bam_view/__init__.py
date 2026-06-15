"""
BAM比对可视化工具包|BAM Alignment Visualization Toolkit
功能: 从BAM文件生成交互式比对可视化|
Features: Generate interactive alignment visualization from BAM files
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-08

使用示例|Usage Examples:
    from biopytools.bam_view import BamViewer, BamViewConfig

    # 创建可视化器|Create viewer
    viewer = BamViewer(
        bam_file="alignments.bam",
        reference="reference.fa",
        region="chr1:1000-2000",
        output_format="html"
    )

    # 运行可视化|Run visualization
    viewer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import BamViewer
from .config import BamViewConfig
from .calculator import BamViewCalculator

__all__ = ['BamViewer', 'BamViewConfig', 'BamViewCalculator']
