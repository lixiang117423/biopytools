"""
SraConverterCLI命令 | SraConverter CLI Command
"""

import click
from ...sra_converter.main import main as sra_converter_main

@click.command()
@click.option('-i', '--input-dir', required=True, help='输入目录 | Input directory')
@click.option('-o', '--output-dir', required=True, help='输出目录 | Output directory')
@click.option('-t', '--threads', default=4, type=int, help='线程数 | Number of threads')
# TODO: 根据实际模块添加更多选项 | Add more options based on actual module
@click.pass_context
def sra_converter(ctx, input_dir, output_dir, threads):
    """SraConverter工具 | SraConverter Tool"""
    import sys
    
    # 构建参数列表 | Build argument list
    args = ['biopytools-sra_converter', '-i', input_dir, '-o', output_dir, '-t', str(threads)]
    
    # 设置sys.argv并调用主函数 | Set sys.argv and call main function
    original_argv = sys.argv
    sys.argv = args
    
    try:
        sra_converter_main()
    finally:
        sys.argv = original_argv
