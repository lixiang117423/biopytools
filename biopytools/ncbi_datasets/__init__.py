"""NCBI datasets 基因组批量下载工具包|NCBI datasets genome batch download toolkit
功能: 输入 taxon 编号, 用官方 datasets CLI 批量下载该 taxon 下所有基因组的序列/注释
|Features: download genome sequences/annotations of all assemblies under a taxon
via the official NCBI datasets CLI
"""

__version__ = "1.0.0"

from .config import NCBIDatasetsConfig
from .downloader import NCBIDatasetsDownloader

__all__ = ['NCBIDatasetsConfig', 'NCBIDatasetsDownloader']
