"""
RagTag基因组scaffolding工具|RagTag Genome Scaffolding Toolkit
功能: 基于参考基因组进行序列scaffolding，支持序列ID重命名和分类输出|
Features: Reference-based genome scaffolding with sequence ID renaming and categorized output
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-01-15

使用示例|Usage Examples:
    from biopytools.ragtag import RagTagScaffolder, RagTagConfig

    # 创建scaffolder|Create scaffolder
    scaffolder = RagTagScaffolder(
        reference="ref.fa",
        query="query.fa",
        sample_name="Sample1",
        threads=12
    )

    # 运行scaffolding|Run scaffolding
    scaffolder.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import RagTagScaffolder
from .config import RagTagConfig

__all__ = ['RagTagScaffolder', 'RagTagConfig']
