"""
群体遗传分析工具包 | Population Genetics Analysis Toolkit
功能: 基于VCF文件的完整群体遗传参数分析套件 |
Features: Complete population genetics parameter analysis suite based on VCF files
作者 | Author: Xiang LI  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-25

使用示例 | Usage Examples:
    from biopytools.popgen_toolkit import PopulationGeneticsAnalyzer, PopGenConfig
    
    # 创建分析器 | Create analyzer
    analyzer = PopulationGeneticsAnalyzer(
        vcf_file="variants.vcf.gz",
        output_dir="popgen_results",
        group_file="groups.txt",
        calculate_all=True
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()

命令行参数对照表 | Command Line Parameter Reference:
    短参数 | Short    长参数 | Long              默认值 | Default       说明 | Description
    -v               --vcf                      -                     输入VCF文件 | Input VCF file
    -o               --output                   popgen_output         输出目录 | Output directory
    -g               --groups                   None                  分组文件 | Group file
    -w               --windows                  [10k,100k,500k]       窗口大小 | Window sizes
    -m               --maf                      0.01                  MAF阈值 | MAF threshold
    -M               --missing                  0.1                   缺失率阈值 | Missing rate threshold
    -H               --hwe                      1e-6                  HWE阈值 | HWE threshold
    -f               --format                   txt                   输出格式 | Output format
    -t               --threads                  4                     线程数 | Threads
    
    --all            --all                      True                  计算所有参数 | Calculate all
    --fst            --fst                      False                 计算Fst | Calculate Fst
    --pi             --pi                       False                 计算π | Calculate π
    --theta-w        --theta-w                  False                 计算θw | Calculate θw
    --tajima-d       --tajima-d                 False                 计算Tajima's D | Calculate Tajima's D
    --ibd            --ibd                      False                 计算IBD | Calculate IBD
    --ld             --ld                       False                 计算LD | Calculate LD
    --ne             --ne                       False                 计算Ne | Calculate Ne
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import PopulationGeneticsAnalyzer
from .config import PopGenConfig

__all__ = ['PopulationGeneticsAnalyzer', 'PopGenConfig']
