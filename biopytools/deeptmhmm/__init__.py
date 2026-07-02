"""
DeepTMHMM 1.0跨膜螺旋/信号肽预测模块|DeepTMHMM 1.0 TM Helix & Signal Peptide Prediction Module

功能|Function: 使用DeepTMHMM预测蛋白质跨膜螺旋、信号肽和拓扑结构
使用示例|Usage Examples:
    from biopytools.deeptmhmm import DeeptmhmmPredictor

    predictor = DeeptmhmmPredictor(
        input_file="proteins.fa",
        output_dir="output"
    )
    predictor.run()
"""

__version__ = "1.0.0"

from .main import DeeptmhmmPredictor, main
from .config import DeeptmhmmConfig

__all__ = ['DeeptmhmmPredictor', 'DeeptmhmmConfig', 'main']
