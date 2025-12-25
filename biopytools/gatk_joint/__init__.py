"""
GATK Joint Genotyping 工具包 | GATK Joint Genotyping Toolkit
功能: VCF/GVCF文件联合分型的完整流程，支持自动识别文件类型、质控过滤和结果整合 | 
Features: Complete pipeline for VCF/GVCF joint genotyping, supporting auto file type detection, QC filtering and result integration
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-10-11

使用示例 | Usage Examples:
    from biopytools.gatk_joint import GATKJointGenotyper, JointConfig
    
    # 创建分析器 | Create genotyper
    genotyper = GATKJointGenotyper(
        input_dir="gvcf_folder",
        reference="ref.fasta",
        output_dir="joint_results",
        threads=88
    )
    
    # 运行分析 | Run analysis
    genotyper.run_pipeline()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import GATKJointGenotyper
from .config import JointConfig

__all__ = ['GATKJointGenotyper', 'JointConfig']
