"""
ADMIXTURE群体结构分析工具包 | ADMIXTURE Population Structure Analysis Toolkit
功能: VCF到ADMIXTURE分析的完整流程，支持群体结构分析和协变量生成 | 
Features: Complete pipeline from VCF to ADMIXTURE analysis, supporting population structure analysis and covariate generation
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化重构版 | Modular refactored version
日期 | Date: 2025-07-17

使用示例 | Usage Examples:
    from biopytools.admixture_toolkit import AdmixtureAnalyzer, AdmixtureConfig
    
    # 创建分析器 | Create analyzer
    analyzer = AdmixtureAnalyzer(
        vcf_file="data.vcf.gz",
        output_dir="admixture_results",
        min_k=2,
        max_k=10
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()

命令行参数对照表 | Command Line Parameter Reference:
    短参数 | Short    长参数 | Long              默认值 | Default    说明 | Description
    -v               --vcf                      -                  输入VCF文件 | Input VCF file
    -o               --output                   admixture_results  输出目录 | Output directory
    -k               --min-k                    2                  最小K值 | Minimum K value
    -K               --max-k                    10                 最大K值 | Maximum K value
    -c               --cv-folds                 5                  交叉验证折数 | CV folds
    -t               --threads                  4                  线程数 | Threads
    -m               --maf                      0.01               MAF阈值 | MAF threshold
    -M               --missing                  0.1                缺失率阈值 | Missing rate threshold
    -H               --hwe                      1e-6               HWE p值阈值 | HWE p-value threshold
    -s               --skip-preprocessing       False              跳过预处理 | Skip preprocessing
    -i               --keep-intermediate        False              保留中间文件 | Keep intermediate files
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import AdmixtureAnalyzer
from .config import AdmixtureConfig

__all__ = ['AdmixtureAnalyzer', 'AdmixtureConfig']
