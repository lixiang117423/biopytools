"""
ALLHiC Pipeline工具包 🧬 | ALLHiC Pipeline Toolkit
功能: Hi-C基因组支架构建和染色体级别组装 | 
Features: Hi-C genome scaffolding and chromosome-level assembly
作者 | Author: Claude  
版本 | Version: v5.4 - Asmkit Edition
日期 | Date: 2025-12-17

使用示例 | Usage Examples:
    from biopytools.allhic import ALLHiCPipeline, PipelineConfig
    
    # 创建流水线 | Create pipeline
    pipeline = ALLHiCPipeline(
        reference="draft.asm.fasta",
        read1="reads_R1.fastq.gz",
        read2="reads_R2.fastq.gz",
        chr_num=12,
        threads=88
    )
    
    # 运行流水线 | Run pipeline
    pipeline.run_pipeline()
"""

__version__ = "5.4.0"
__author__ = "Claude"

import sys
import os
from pathlib import Path

# 添加当前目录到sys.path，解决相对导入问题
package_dir = Path(__file__).parent
if str(package_dir) not in sys.path:
    sys.path.insert(0, str(package_dir))

from config import PipelineConfig
from main import ALLHiCPipeline

__all__ = ['ALLHiCPipeline', 'PipelineConfig']
