"""
FASTQ文件过滤命令|FASTQ File Filtering Command
"""

import click
import sys
import os


def _lazy_import_fastq_gc_filter_main():
    """延迟加载fastq_gc_filter主函数|Lazy load fastq_gc_filter main function"""
    try:
        from ...fastq_gc_filter.main import main as fastq_gc_filter_main
        return fastq_gc_filter_main
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
    short_help='FASTQ文件GC含量和序列长度过滤|FASTQ file GC content and sequence length filtering',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTQ文件|Input FASTQ file (支持.gz压缩|supports .gz compression)')
@click.option('--output', '-o',
              required=True,
              help='输出FASTQ文件|Output FASTQ file (支持.gz压缩|supports .gz compression)')
@click.option('--min-gc',
              type=float,
              default=25.0,
              show_default=True,
              help='最小GC含量百分比|Minimum GC content percentage')
@click.option('--max-gc',
              type=float,
              default=100.0,
              show_default=True,
              help='最大GC含量百分比|Maximum GC content percentage')
@click.option('--min-length',
              type=int,
              default=50,
              show_default=True,
              help='最短序列长度|Minimum sequence length')
@click.option('--max-length',
              type=int,
              default=None,
              help='最长序列长度|Maximum sequence length')
def fastq_gc_filter(input, output, min_gc, max_gc, min_length, max_length):
    """
    FASTQ文件GC含量和序列长度过滤工具|FASTQ File GC Content and Sequence Length Filtering Tool

    根据GC含量和序列长度筛选FASTQ文件中的reads|Filter FASTQ reads by GC content and sequence length

    示例|Examples: biopytools fastq-gc-filter -i input.fastq -o filtered.fastq --min-gc 25 --min-length 50
    """

    # 延迟加载|Lazy loading
    fastq_gc_filter_main = _lazy_import_fastq_gc_filter_main()

    # 构建参数列表|Build argument list
    args = ['fastq_gc_filter.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # GC含量参数|GC content parameters
    if min_gc != 25.0:
        args.extend(['--min-gc', str(min_gc)])

    if max_gc != 100.0:
        args.extend(['--max-gc', str(max_gc)])

    # 序列长度参数|Sequence length parameters
    if min_length != 50:
        args.extend(['--min-length', str(min_length)])

    if max_length is not None:
        args.extend(['--max-length', str(max_length)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        fastq_gc_filter_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
