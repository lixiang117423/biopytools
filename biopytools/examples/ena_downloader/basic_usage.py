#!/usr/bin/env python3
"""
ENA下载工具基本使用示例 | ENA Downloader Basic Usage Examples
"""

from ena_downloader import ENADownloader

def example_metadata_only():
    """示例1: 仅下载元数据 | Example 1: Download metadata only"""
    print("=" * 50)
    print("示例1: 仅下载元数据 | Example 1: Metadata only")
    print("=" * 50)
    
    downloader = ENADownloader(
        accession="PRJNA661210",
        output_dir="example_metadata_only",
        metadata_format="tsv"
    )
    
    result = downloader.download_metadata()
    if result:
        print(f"元数据文件: {result}")
    else:
        print("下载失败")

def example_ftp_script():
    """示例2: 生成FTP下载脚本 | Example 2: Generate FTP download script"""
    print("=" * 50)
    print("示例2: 生成FTP下载脚本 | Example 2: FTP download script")
    print("=" * 50)
    
    downloader = ENADownloader(
        accession="PRJNA661210",
        protocol="ftp",
        method="save",
        output_dir="example_ftp_script"
    )
    
    success = downloader.run_full_pipeline()
    print(f"脚本生成{'成功' if success else '失败'}")

def example_custom_fields():
    """示例3: 自定义字段 | Example 3: Custom fields"""
    print("=" * 50)
    print("示例3: 自定义字段 | Example 3: Custom fields")
    print("=" * 50)
    
    downloader = ENADownloader(
        accession="PRJNA661210",
        fields=["fastq_ftp", "study_title", "sample_accession", "scientific_name"],
        metadata_format="xlsx",
        output_dir="example_custom_fields"
    )
    
    result = downloader.download_metadata()
    if result:
        print(f"自定义字段元数据文件: {result}")

def example_aspera_script():
    """示例4: 生成Aspera下载脚本 | Example 4: Generate Aspera download script"""
    print("=" * 50)
    print("示例4: 生成Aspera下载脚本 | Example 4: Aspera download script")
    print("=" * 50)
    
    # 注意: 需要有效的Aspera密钥路径 | Note: Requires valid Aspera key path
    aspera_key = "~/.aspera/connect/etc/asperaweb_id_dsa.openssh"
    
    try:
        downloader = ENADownloader(
            accession="PRJNA661210",
            protocol="aspera",
            method="save",
            aspera_key=aspera_key,
            output_dir="example_aspera_script"
        )
        
        success = downloader.run_full_pipeline()
        print(f"Aspera脚本生成{'成功' if success else '失败'}")
        
    except Exception as e:
        print(f"Aspera配置错误: {e}")

if __name__ == "__main__":
    print("ENA下载工具使用示例 | ENA Downloader Usage Examples")
    print()
    
    # 运行示例 | Run examples
    example_metadata_only()
    print()
    
    example_ftp_script()
    print()
    
    example_custom_fields()
    print()
    
    example_aspera_script()
    print()
    
    print("所有示例执行完成 | All examples completed")
