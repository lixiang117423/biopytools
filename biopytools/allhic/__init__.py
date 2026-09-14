"""
ALLHiC Pipeline工具包 |ALLHiC Pipeline Toolkit
功能: Hi-C基因组支架构建和染色体级别组装|
Features: Hi-C genome scaffolding and chromosome-level assembly
作者|Author: Claude  
版本|Version: v5.4 - Asmkit Edition
日期|Date: 2025-12-17

使用示例|Usage Examples:
    from biopytools.allhic import ALLHiCPipeline, PipelineConfig
    
    # 创建流水线|Create pipeline
    pipeline = ALLHiCPipeline(
        reference="draft.asm.fasta",
        read1="reads_R1.fastq.gz",
        read2="reads_R2.fastq.gz",
        chr_num=12,
        threads=12
    )
    
    # 运行流水线|Run pipeline
    pipeline.run_pipeline()
"""

__version__ = "5.4.0"
__author__ = "Claude"

# 包内统一相对导入(§一): 不得用 sys.path hack + 顶层绝对导入,
# 否则 config.py 的 `..common.paths` 相对导入会失败(import biopytools.allhic 即报
# "attempted relative import with no known parent package")
# |Intra-package relative imports only; the old sys.path hack made
# config.py's `..common.paths` fail at import time
from .config import PipelineConfig
from .main import ALLHiCPipeline

__all__ = ['ALLHiCPipeline', 'PipelineConfig']
