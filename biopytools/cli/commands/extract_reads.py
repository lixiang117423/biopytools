"""
基于contig-reads对应关系提取fastq reads|Extract fastq reads by contig-reads mapping
"""

import click
import sys
import os


def _lazy_import_extract_reads_main():
    """延迟加载extract_reads主函数|Lazy load extract_reads main function"""
    try:
        from ...extract_reads.main import main as extract_reads_main
        return extract_reads_main
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
    short_help='基于contig-reads对应关系提取fastq reads|Extract fastq reads by contig-reads mapping',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--mapping', '-m',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='contig-reads对应关系文件(TSV格式)|contig-reads mapping file (TSV format)')
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTQ文件|Input FASTQ file')
@click.option('--output', '-o',
              required=True,
              help='输出文件|Output file')
@click.option('--no-compress',
              is_flag=True,
              help='不压缩输出文件|Do not compress output files')
def extract_reads(mapping, input, output, no_compress):
    """
    基于contig-reads对应关系从fastq提取指定reads|Extract specified reads from fastq by contig-reads mapping

    示例|Example: biopytools extract-reads -m contig_reads.tsv -i input.fq.gz -o output.fq.gz
    """

    # 延迟加载|Lazy loading
    extract_reads_main = _lazy_import_extract_reads_main()

    # 构建参数列表|Build argument list
    args = ['extract_reads.py']

    # 必需参数|Required parameters
    args.extend(['-m', mapping])
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数|Optional parameters
    if no_compress:
        args.append('--no-compress')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        extract_reads_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
