"""
ENA数据下载工具包 | ENA Data Downloader Toolkit
功能: 从ENA (European Nucleotide Archive) 下载元数据和FASTQ文件的完整工具 | 
Features: Complete toolkit for downloading metadata and FASTQ files from ENA (European Nucleotide Archive)
作者 | Author: Claude  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-23

使用示例 | Usage Examples:
    from biopytools.ena_downloader import ENADownloader, DownloadConfig
    
    # 只下载元数据 | Download metadata only
    downloader = ENADownloader(
        accession="PRJNA661210",
        output_dir="ena_results"
    )
    downloader.download_metadata()
    
    # 下载元数据并生成下载脚本 | Download metadata and generate download script
    downloader = ENADownloader(
        accession="PRJNA661210",
        protocol="ftp",
        method="save",
        output_dir="ena_results"
    )
    downloader.run_full_pipeline()
"""

__version__ = "1.0.0"
__author__ = "Claude"

from .main import ENADownloader
from .config import DownloadConfig

__all__ = ['ENADownloader', 'DownloadConfig']
