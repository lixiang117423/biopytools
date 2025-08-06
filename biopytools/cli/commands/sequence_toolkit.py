"""
SequenceToolkitCLI命令 | SequenceToolkit CLI Command
"""

import click
from ...sequence_toolkit.main import main as sequence_toolkit_main

@click.command()
@click.option('-i', '--input-dir', required=True, help='输入目录 | Input directory')
@click.option('-o', '--output-dir', required=True, help='输出目录 | Output directory')
@click.option('-t', '--threads', default=4, type=int, help='线程数 | Number of threads')
# TODO: 根据实际模块添加更多选项 | Add more options based on actual module
@click.pass_context
def sequence_toolkit(ctx, input_dir, output_dir, threads):
    """SequenceToolkit工具 | SequenceToolkit Tool"""
    import sys
    
    # 构建参数列表 | Build argument list
    args = ['biopytools-sequence_toolkit', '-i', input_dir, '-o', output_dir, '-t', str(threads)]
    
    # 设置sys.argv并调用主函数 | Set sys.argv and call main function
    original_argv = sys.argv
    sys.argv = args
    
    try:
        sequence_toolkit_main()
    finally:
        sys.argv = original_argv
