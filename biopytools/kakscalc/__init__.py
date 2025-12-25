"""
🧬 Ka/Ks Calculator 2.0 Python封装
功能: 基于KaKs_Calculator2.0的同义和非同义替换率分析工具 | Ka/Ks substitution rate analysis toolkit
支持: 多种计算方法和批量分析 | Multiple calculation methods and batch analysis
作者 | Author: Bioinformatics Pipeline Team
版本 | Version: v1.0.0 - Python模块化版本 | Python modular version
日期 | Date: 2025-08-15

主要功能 | Main Features:
- 🧮 Ka/Ks比率计算 | Ka/Ks ratio calculation
- 🔬 17种计算方法支持 | 17 calculation methods supported
- 📊 批量序列对分析 | Batch sequence pair analysis
- ✅ 智能输入验证 | Intelligent input validation
- 📈 结果统计分析 | Statistical result analysis
- 📋 多格式输出 | Multiple output formats

使用示例 | Usage Examples:
    from biopytools.kakscalc import KaKsAnalyzer, KaKsConfig
    
    # 创建分析器 | Create analyzer
    analyzer = KaKsAnalyzer(
        fasta1="species1.fasta",
        fasta2="species2.fasta", 
        pairs="pairs.txt",
        output_dir="results/"
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()
"""

__version__ = "1.0.0"
__author__ = "Bioinformatics Pipeline Team"

from .main import KaKsAnalyzer
from .config import KaKsConfig
from .logger import Logger

__all__ = ['KaKsAnalyzer', 'KaKsConfig', 'Logger']
