"""
染色体重命名命令|Chromosome Rename Command
"""

import click
import sys
import os


def _lazy_import_rename_main():
    """延迟加载染色体重命名主函数|Lazy load chromosome rename main function"""
    try:
        from ...rename_chromosomes.main import main as rename_main
        return rename_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if _is_help_request():
        return file_path
    if not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='染色体重命名工具|Chromosome rename tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件路径|Input FASTA file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出FASTA文件路径|Output FASTA file path')
@click.option('--number', '-n',
              required=True,
              type=int,
              help='染色体数量|Number of chromosomes')
@click.option('--keep-all',
              is_flag=True,
              default=False,
              help='输出所有序列（默认只输出前n个序列）|Output all sequences (default: only first n)')
def rename_chromosomes(input, output, number, keep_all):
    """
    染色体重命名工具|Chromosome Rename Tool

    将FASTA文件中的序列名重命名为标准染色体格式|Rename sequences in FASTA file to standard chromosome format

    示例|Examples: biopytools rename-chromosomes -i genome.fa -o renamed.fa -n 20
    """

    # 延迟加载|Lazy load
    rename_main = _lazy_import_rename_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['rename_chromosomes.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])
    args.extend(['-n', str(number)])

    # 可选参数|Optional parameters
    if keep_all:
        args.append('--keep-all')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        rename_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
