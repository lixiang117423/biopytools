"""
FOF文件生成命令|FOF File Generation Command
"""

import click
import sys
import os


def _lazy_import_fof_main():
    """延迟加载FOF生成主函数|Lazy load FOF generation main function"""
    try:
        from ...fof.main import main as fof_main
        return fof_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式，支持软链接)|Validate path exists (only in non-help mode, supports symlinks)"""
    if not _is_help_request() and path and not os.path.lexists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='生成样品名到文件路径的FOF映射表|Generate sample-to-file-path FOF mapping table',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入文件或目录|Input file or directory')
@click.option('-o', '--output',
              required=True,
              type=click.Path(),
              help='输出FOF文件路径|Output FOF file path')
@click.option('-s', '--suffix',
              multiple=True,
              help='文件后缀过滤（可多次指定，如 -s .fastq.gz -s .fq.gz）|'
                   'File suffix filter (can be specified multiple times)')
@click.option('-r', '--recursive',
              is_flag=True,
              default=False,
              help='递归扫描子目录|Recursively scan subdirectories')
def fof(input, output, suffix, recursive):
    """
    FOF文件生成工具|FOF File Generation Tool

    生成样品名到文件绝对路径的tab分割映射表|Generate tab-separated sample-to-file-path mapping table

    示例|Examples: biopytools fof -i ./data/ -o samples.fof
    """

    # 延迟加载|Lazy loading
    fof_main = _lazy_import_fof_main()

    # 构建参数列表|Build argument list
    args = ['fof.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数|Optional parameters
    for s in suffix:
        args.extend(['-s', s])

    if recursive:
        args.extend(['-r'])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        fof_main()
    except SystemExit as e:
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断操作|Analysis interrupted by user", err=True)
        sys.exit(130)
    except Exception as e:
        click.echo(f"分析执行失败|Analysis execution failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
