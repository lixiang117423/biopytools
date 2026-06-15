"""
OrthoFinder泛基因组分析工具包|OrthoFinder Pangenome Analysis Toolkit
功能: 基于OrthoFinder的泛基因组分析完整流程，支持核心基因组、附属基因组和特异基因组分类|
Features: Complete pangenome analysis pipeline based on OrthoFinder, supporting core, accessory and specific genome classification
版本|Version: v1.0.0

使用示例|Usage Examples:
    from biopytools.orthofinder import OrthoFinderPangenomeAnalyzer, PangenomeConfig

    config = PangenomeConfig(
        input_dir="protein_sequences/",
        output_dir="pangenome_results",
        softcore_missing_threshold=1,
        search_program="blast"
    )
    analyzer = OrthoFinderPangenomeAnalyzer(**config.__dict__)
    analyzer.run_analysis()
"""

__version__ = "1.0.0"

from .main import OrthoFinderPangenomeAnalyzer
from .config import PangenomeConfig

__all__ = ['OrthoFinderPangenomeAnalyzer', 'PangenomeConfig']
