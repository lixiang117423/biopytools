"""选择性扫荡检测模块|Selective Sweep Detection Module

输入VCF+群体分组信息,计算pi/Tajima's D/Fst/RAiSD μ统计量,
合并为composite_score并输出候选选择性扫荡区域。
|Input VCF + pop info; compute pi/TajimaD/Fst/RAiSD mu statistics,
merge into composite_score, output candidate selective sweep regions.
"""

from .config import SweepMergeConfig
from .merge_stats import SweepMerger
from .pipeline import SweepPipeline

__version__ = '1.0.0'

__all__ = ['SweepMergeConfig', 'SweepMerger', 'SweepPipeline']
