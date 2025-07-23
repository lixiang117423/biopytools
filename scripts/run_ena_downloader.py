#!/usr/bin/env python3
"""
ENA数据下载运行脚本 | ENA Data Download Runner Script
这是一个简化的入口脚本，用于运行ENA数据下载 | Simple entry script for running ENA data download

用法 | Usage:
    python run_ena_downloader.py --accession PRJNA661210 --metadata-only
    
示例 | Examples:
    # 仅下载元数据 | Download metadata only
    python run_ena_downloader.py -a PRJNA661210 -M
    
    # 下载元数据并生成FTP下载脚本 | Download metadata and generate FTP download script
    python run_ena_downloader.py -a PRJNA661210 -p ftp -m save
    
    # 下载元数据并生成Aspera下载脚本 | Download metadata and generate Aspera download script
    python run_ena_downloader.py -a PRJNA661210 -p aspera -k ~/.aspera/connect/etc/asperaweb_id_dsa.openssh
    
    # 直接执行下载 | Execute download directly
    python run_ena_downloader.py -a PRJNA661210 -p ftp -m run
    
    # 输出Excel格式的元数据 | Output Excel format metadata
    python run_ena_downloader.py -a PRJNA661210 -f xlsx -o results
    
    # 自定义元数据字段 | Custom metadata fields
    python run_ena_downloader.py -a PRJNA661210 -F fastq_ftp study_title sample_accession
"""

from biopytools.ena_downloader.main import main

if __name__ == "__main__":
    main()
