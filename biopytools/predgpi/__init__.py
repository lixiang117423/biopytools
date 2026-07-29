"""
PredGPI GPI锚定蛋白预测模块|PredGPI GPI-anchored protein prediction module
"""

# 仅导出类,不导出main函数,避免与子模块main.py同名属性冲突
from .main import PredGPIPredictor
from .config import PredGPIConfig

__version__ = "1.0.0"

__all__ = ["PredGPIConfig", "PredGPIPredictor"]
