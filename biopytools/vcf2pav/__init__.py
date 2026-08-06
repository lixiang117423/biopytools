"""VCF转PAV矩阵|VCF to PAV matrix converter"""

from .config import Vcf2PavConfig

# converter.py 由后续任务创建，此处守卫导入保证骨架阶段包可导入、
# 且 converter 落盘后无需改动本文件即可自动导出
# |converter.py lands in a later task; guarding the import keeps the package
# importable during the skeleton phase, and auto-exports once it exists
try:
    from .converter import Vcf2PavConverter
except ImportError:
    Vcf2PavConverter = None

__version__ = "1.0.0"

__all__ = ["Vcf2PavConfig", "Vcf2PavConverter"]
