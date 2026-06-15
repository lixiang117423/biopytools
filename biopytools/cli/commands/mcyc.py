"""
甲烷循环基因丰度分析|Methane Cycle Gene Abundance Analysis Command
"""

import click
import sys
import os


def _lazy_import_mcyc_main():
    """延迟加载mcyc主函数|Lazy load mcyc main function"""
    try:
        from ...mcyc.main import main as mcyc_main
        return mcyc_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


def _validate_threads(threads):
    """验证线程数|Validate thread count"""
    if threads < 1:
        raise click.BadParameter("线程数必须大于0|Thread count must be greater than 0")
    return threads


@click.command(
    short_help='甲烷循环基因丰度分析|Methane cycle gene abundance analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.argument('input-list',
                callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
                required=True)
@click.option('--output', '-o',
              default=None,
              show_default=True,
              type=click.Path(),
              help='输出目录路径|Output directory path')
@click.option('--threads', '-t',
              default=12,
              show_default=True,
              type=int,
              callback=lambda ctx, param, value: _validate_threads(value) if value is not None else None,
              help='线程数|Thread count')
@click.option('--mcyc-base',
              default=None,
              type=click.Path(),
              help='MCycDB基础路径|MCycDB base path')
@click.option('--skip-diamond',
              is_flag=True,
              help='跳过Diamond比对|Skip Diamond alignment')
@click.option('--keep-temp',
              is_flag=True,
              help='保留临时文件|Keep temporary files')
def mcyc(input_list, output, threads, mcyc_base, skip_diamond, keep_temp):
    """
    甲烷循环基因丰度分析工具|Methane Cycle Gene Abundance Analysis Tool

    基于MCycDB数据库计算甲烷循环基因丰度(TPM,CLR)
    Calculate methane cycle gene abundance (TPM, CLR) based on MCycDB database

    示例|Examples: biopytools mcyc samples.txt
    """

    # 延迟加载|Lazy load
    mcyc_main = _lazy_import_mcyc_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['mcyc.py']

    # 必需参数|Required parameter
    args.extend([input_list])

    # 可选参数|Optional parameters
    if output is not None:
        args.extend(['--output', output])

    if threads != 4:
        args.extend(['--threads', str(threads)])

    if mcyc_base is not None:
        args.extend(['--mcyc-base', mcyc_base])

    if skip_diamond:
        args.append('--skip-diamond')

    if keep_temp:
        args.append('--keep-temp')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        mcyc_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
