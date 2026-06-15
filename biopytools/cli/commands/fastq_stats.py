"""
FASTQ文件统计命令|FASTQ File Statistics Command
"""

import click
import sys
import os


def _lazy_import_fastq_stats_main():
    """延迟加载fastq_stats主函数|Lazy load fastq_stats main function"""
    try:
        from ...fastq_stats.main import main as fastq_stats_main
        return fastq_stats_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if path and not _is_help_request() and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='FASTQ文件统计工具|FASTQ file statistics tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入FASTQ文件或目录|Input FASTQ file or directory')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出文件路径(.csv或.xlsx)|Output file path (.csv or .xlsx)')
@click.option('--pattern', '-p',
              default=None,
              help='FASTQ文件匹配模式，如"*_1.clean.fq.gz"|FASTQ file matching pattern, e.g., "*_1.clean.fq.gz"')
@click.option('--threads', '-t',
              default=12,
              type=int,
              show_default=True,
              help='线程数|Number of threads')
def fastq_stats(input, output, pattern, threads):
    """
    FASTQ文件统计工具|FASTQ File Statistics Tool

    基于seqkit的高性能FASTQ文件统计工具，支持自动配对双末端文件，输出CSV/Excel格式
    High-performance FASTQ file statistics tool based on seqkit, supports automatic paired-end file matching, outputs CSV/Excel format

    示例|Examples: biopytools fq-stats -i /data/fastq/ -o results.csv -p "*_1.clean.fq.gz"
    """

    # 延迟加载|Lazy loading
    fastq_stats_main = _lazy_import_fastq_stats_main()

    # 构建参数列表|Build argument list
    args = ['fastq_stats.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数|Optional parameters
    if pattern:
        args.extend(['-p', pattern])

    if threads != 12:
        args.extend(['-t', str(threads)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        fastq_stats_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
