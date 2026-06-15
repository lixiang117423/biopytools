"""
CPhasing基因组分相和挂载工具包|CPhasing Genome Phasing and Scaffolding Toolkit

功能|Features:
- 基于Hi-C/Pore-C数据的多倍体基因组分相|Polyploid genome phasing based on Hi-C/Pore-C data
- 染色体级别挂载|Chromosome-level scaffolding
- 支持所有CPhasing子命令|Supports all CPhasing subcommands
- 完整的conda环境管理和日志记录|Complete conda env management and logging

使用示例|Usage Examples:
    >>> from biopytools.cphasing import CPhasingConfig, CPhasingRunner

    >>> # pipeline模式（默认）|pipeline mode (default)
    >>> config = CPhasingConfig(
    ...     fasta="genome.fa",
    ...     hic1="R1.fq.gz",
    ...     hic2="R2.fq.gz",
    ...     groups="8:4",
    ...     threads=10
    ... )
    >>> runner = CPhasingRunner(config)
    >>> runner.run()

    >>> # 其他子命令|Other subcommands
    >>> config = CPhasingConfig(
    ...     subcommand="alleles",
    ...     fasta="genome.fa",
    ...     extra_args=["-t", "12"]
    ... )
    >>> runner = CPhasingRunner(config)
    >>> runner.run()
"""

__version__ = "1.0.0"

from .config import CPhasingConfig
from .runner import CPhasingRunner, CPhasingPipeline
from .main import main

__all__ = [
    'CPhasingConfig',
    'CPhasingRunner',
    'CPhasingPipeline',
    'main',
]
