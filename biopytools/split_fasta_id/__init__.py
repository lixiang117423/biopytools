"""
🧬 FASTA序列ID分割工具包 | FASTA Sequence ID Splitting Toolkit
功能: 按空格或制表符分割FASTA序列名称行，提取指定位置的元素
Features: Split FASTA sequence name lines by space or tab, extract element at specified position
"""

from .main import FastaIDSplitter
from .config import SplitConfig

__all__ = ['FastaIDSplitter', 'SplitConfig']
