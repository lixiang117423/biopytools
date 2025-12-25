"""
基因组组装工具包 | Genome Assembly Toolkit 🧬
功能: 多数据类型整合的染色体级基因组组装完整流程 | 
Features: Complete pipeline for chromosome-level genome assembly with multi-data integration
作者 | Author: GenomeAssembly Team  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-08-29

使用示例 | Usage Examples:
    from biopytools.genomeasm import GenomeAssembler, AssemblyConfig
    
    # 创建组装器 | Create assembler
    assembler = GenomeAssembler(
        input_dir="raw_data/",
        output_dir="assembly_results/",
        project_name="genome_project",
        hic_strategy="complete_juicer"
    )
    
    # 运行分析 | Run analysis
    assembler.run_assembly()
"""

__version__ = "1.0.0"
__author__ = "GenomeAssembly Team"

# from .main import GenomeAssembler
# from .config import AssemblyConfig

# __all__ = ['GenomeAssembler', 'AssemblyConfig']

from .main import GenomeAssembler, main  # 添加main函数导入
from .config import AssemblyConfig

__all__ = ['GenomeAssembler', 'AssemblyConfig', 'main']  # 添加main到导出列表
