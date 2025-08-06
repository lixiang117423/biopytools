"""
BismarkPipelineeCLI命令 | BismarkPipelinee CLI Command
"""

import click
from ...bismark_pipelinee.main import main as bismark_pipelinee_main

@click.command()
@click.option('-i', '--input-dir', required=True, help='输入目录 | Input directory')
@click.option('-o', '--output-dir', required=True, help='输出目录 | Output directory')
@click.option('-t', '--threads', default=4, type=int, help='线程数 | Number of threads')
# TODO: 根据实际模块添加更多选项 | Add more options based on actual module
@click.pass_context
def bismark_pipelinee(ctx, input_dir, output_dir, threads):
    """BismarkPipelinee工具 | BismarkPipelinee Tool"""
    import sys
    
    # 构建参数列表 | Build argument list
    args = ['biopytools-bismark_pipelinee', '-i', input_dir, '-o', output_dir, '-t', str(threads)]
    
    # 设置sys.argv并调用主函数 | Set sys.argv and call main function
    original_argv = sys.argv
    sys.argv = args
    
    try:
        bismark_pipelinee_main()
    finally:
        sys.argv = original_argv
