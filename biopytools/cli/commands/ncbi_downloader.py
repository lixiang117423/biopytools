"""
NcbiDownloaderCLI命令 | NcbiDownloader CLI Command
"""

import click
from ...ncbi_downloader.main import main as ncbi_downloader_main

@click.command()
@click.option('-i', '--input-dir', required=True, help='输入目录 | Input directory')
@click.option('-o', '--output-dir', required=True, help='输出目录 | Output directory')
@click.option('-t', '--threads', default=4, type=int, help='线程数 | Number of threads')
# TODO: 根据实际模块添加更多选项 | Add more options based on actual module
@click.pass_context
def ncbi_downloader(ctx, input_dir, output_dir, threads):
    """NcbiDownloader工具 | NcbiDownloader Tool"""
    import sys
    
    # 构建参数列表 | Build argument list
    args = ['biopytools-ncbi_downloader', '-i', input_dir, '-o', output_dir, '-t', str(threads)]
    
    # 设置sys.argv并调用主函数 | Set sys.argv and call main function
    original_argv = sys.argv
    sys.argv = args
    
    try:
        ncbi_downloader_main()
    finally:
        sys.argv = original_argv
