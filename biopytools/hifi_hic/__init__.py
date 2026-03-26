"""
基因组组装工具包 |Genome Assembly Toolkit
功能: 使用 hifiasm 进行 HiFi + Hi-C 基因组组装|
Features: HiFi + Hi-C genome assembly using hifiasm
作者|Author: Claude  
版本|Version: v1.0 - 基础版本|Basic version
日期|Date: 2025-12-17

使用示例|Usage Examples:
    from biopytools.genome_assembler import GenomeAssembler, AssemblyConfig
    
    # 创建组装器|Create assembler
    assembler = GenomeAssembler(
        hifi_data="hifi.fq",
        hic_r1="hic_R1.fastq.gz",
        hic_r2="hic_R2.fastq.gz",
        prefix="sample1",
        threads=12
    )
    
    # 运行组装|Run assembly
    assembler.run_assembly()
"""

__version__ = "1.0.0"
__author__ = "Claude"

import sys
import os
from pathlib import Path

# 添加当前目录到sys.path，解决相对导入问题
package_dir = Path(__file__).parent
if str(package_dir) not in sys.path:
    sys.path.insert(0, str(package_dir))

from config import AssemblyConfig
from main import GenomeAssembler

__all__ = ['GenomeAssembler', 'AssemblyConfig']
