"""
🌳 RAxML系统发育分析工具包 | RAxML Phylogenetic Analysis Toolkit
功能: RAxML最大似然系统发育树构建的完整流程 | 
Features: Complete pipeline for RAxML maximum likelihood phylogenetic tree construction
作者 | Author: Bioinformatics Pipeline Team  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-08-24

使用示例 | Usage Examples:
    from biopytools.raxml import RAxMLAnalyzer, RAxMLConfig
    
    # 创建分析器 | Create analyzer
    analyzer = RAxMLAnalyzer(
        sequence_file="alignment.phy",
        output_name="phylo_tree",
        model="GTRGAMMA",
        threads=88
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Bioinformatics Pipeline Team"

from .main import RAxMLAnalyzer
from .config import RAxMLConfig

__all__ = ['RAxMLAnalyzer', 'RAxMLConfig']
