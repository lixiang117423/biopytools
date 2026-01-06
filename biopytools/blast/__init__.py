"""
BLAST序列比对分析工具包|BLAST Alignment Analysis Toolkit
功能: 多序列文件与目标序列的BLAST比对分析完整流程|Complete pipeline for BLAST alignment analysis
作者|Author: Xiang LI
版本|Version: 2.0.0 - 标准化版本|Standardized version
日期|Date: 2025-12-19

使用示例|Usage Examples:
    from biopytools.blast import BLASTAnalyzer, BLASTConfig

    # 创建分析器|Create analyzer
    config = BLASTConfig(
        input="/path/to/sequences.fa",
        target="/path/to/target.fa",
        output_dir="./blast_results"
    )
    analyzer = BLASTAnalyzer(config)

    # 运行分析|Run analysis
    analyzer.run_pipeline()
"""

__version__ = "2.0.0"
__author__ = "Xiang LI"

from .main import BLASTAnalyzer
from .config import BLASTConfig

__all__ = ['BLASTAnalyzer', 'BLASTConfig']