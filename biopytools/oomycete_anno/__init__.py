"""疫霉菌基因组注释模块|Oomycete Genome Annotation Module

复刻 T2T 文章的"证据驱动 Augustus"手工流程:
RepeatMasker -> HISAT2(二代RNA-seq) -> bam2hints -> GeneMark训练 -> Augustus预测
|Replicates the T2T evidence-driven Augustus pipeline:
RepeatMasker -> HISAT2(short RNA-seq) -> bam2hints -> GeneMark training -> Augustus prediction

按手头证据 graceful degradation: 缺哪种证据就跳过对应步骤。
|Graceful degradation: skip steps whose evidence is absent.
"""

from .config import OomyceteAnnoConfig
from .pipeline import OomyceteAnnoRunner

__version__ = "1.0.0"

__all__ = [
    "OomyceteAnnoConfig",
    "OomyceteAnnoRunner",
]
