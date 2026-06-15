"""
基因组挂载率统计命令|Genome Mount Rate Statistics Command

"""

import click
import sys
import os


def _lazy_import_genome_mount_rate_main():
    """延迟加载genome_mount_rate主函数|Lazy load genome_mount_rate main function"""
    try:
        from ...genome_mount_rate.main import main as genome_mount_rate_main
        return genome_mount_rate_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='基因组挂载率统计工具|Genome mount rate statistics tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件|Input FASTA file path')
@click.option('--number', '-n',
              required=True,
              type=int,
              help='序列数量|Number of sequences')
@click.option('--sort',
              is_flag=True,
              help='按长度从大到小排序(计算最长N条)|Sort by length descending (calculate longest N)')
def genome_mount_rate(input, number, sort):
    """
    基因组挂载率统计工具|Genome Mount Rate Statistics Tool

    计算FASTA文件中前N条（或最长N条）序列占总基因组长度的百分比|Calculate percentage of top N (or longest N) sequences

    示例|Example: biopytools genome-mount-rate -i genome.fa -n 10 --sort
    """

    # 延迟加载|Lazy loading
    genome_mount_rate_main = _lazy_import_genome_mount_rate_main()

    # 构建参数列表|Build argument list
    args = ['genome_mount_rate.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-n', str(number)])

    # 可选参数|Optional parameters
    if sort:
        args.append('--sort')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        genome_mount_rate_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
