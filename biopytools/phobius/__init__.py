"""
Phobius跨膜拓扑+信号肽预测模块|Phobius TM topology & signal peptide prediction module
"""

# 仅导出类(不导出main函数,避免与子模块main.py同名属性冲突,保证 patch("...phobius.main.x") 可用)
# Export classes only (not the main fn) to avoid shadowing the main submodule,
# keeping dotted patches like patch("...phobius.main.x") resolvable.
from .main import PhobiusPredictor
from .config import PhobiusConfig

__version__ = "1.0.0"

__all__ = ["PhobiusConfig", "PhobiusPredictor"]
