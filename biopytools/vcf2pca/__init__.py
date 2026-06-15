"""
VCF2PCA分析工具包|VCF2PCA Analysis Toolkit
功能: VCF文件主成分分析，支持VCF2PCACluster和PLINK两种后端|
Features: VCF principal component analysis supporting VCF2PCACluster and PLINK backends
作者|Author: Xiang LI
版本|Version: v2.0 - 双backend版本|Dual-backend version
日期|Date: 2026-03-23

使用示例|Usage Examples:
    from biopytools.vcf2pca import VCF2PCARunner, VCF2PCAConfig

    # 使用VCF2PCACluster后端（默认）|Use VCF2PCACluster backend (default)
    runner = VCF2PCARunner(
        vcf_file="variants.vcf",
        output_dir="pca_results",
        backend="v2p",
        components=10,
        cluster=True,
        cluster_method="kmeans",
        cluster_k=3
    )

    # 使用PLINK后端|Use PLINK backend
    runner = VCF2PCARunner(
        vcf_file="variants.vcf",
        output_dir="pca_results",
        backend="plink",
        components=10,
        maf=0.05,
        plot=True
    )

    # 运行分析|Run analysis
    runner.run_analysis()
"""

__version__ = "2.0.0"
__author__ = "Xiang LI"

from .main import VCF2PCARunner
from .config import VCF2PCAConfig
from .v2p_analyzer import V2PAnalyzer

# 兼容性保留|Legacy compatibility
from .main import main
from .utils import VCF2PCALogger, PCALogger, CommandRunner, check_dependencies

__all__ = [
    'VCF2PCARunner',
    'VCF2PCAConfig',
    'V2PAnalyzer',
    'VCF2PCALogger',
    'PCALogger',
    'CommandRunner',
    'check_dependencies',
    'main'
]
