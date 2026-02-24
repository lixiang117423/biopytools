"""
公共测序数据下载工具包|Public Sequencing Data Download Toolkit
功能: 基于iSeq工具从GSA/SRA/ENA/DDBJ数据库下载测序数据和元数据|
Features: Download sequencing data and metadata from GSA/SRA/ENA/DDBJ databases based on iSeq tool
作者|Author: Xiang LI
版本|Version: 1.0.0
日期|Date: 2026-02-24

使用示例|Usage Examples:
    from biopytools.iseq import ISeqDownloader, ISeqConfig

    # 创建下载器|Create downloader
    downloader = ISeqDownloader(
        accession="PRJNA1014406",
        output_dir="./data"
    )

    # 运行下载|Run download
    downloader.run()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import ISeqDownloader
from .config import ISeqConfig

__all__ = ['ISeqDownloader', 'ISeqConfig']
