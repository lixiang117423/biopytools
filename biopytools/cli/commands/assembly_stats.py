"""
基因组装配统计命令|Genome Assembly Statistics Command
"""

import click
import sys
import os


def _lazy_import_assembly_stats_main():
    """延迟加载assembly_stats主函数|Lazy load assembly_stats main function"""
    try:
        from ...assembly_stats.main import main as assembly_stats_main
        return assembly_stats_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if not _is_help_request():
        if not os.path.exists(path):
            raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='基因组装配统计|Genome Assembly Statistics',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入文件或文件夹|Input file or directory path')
@click.option('-l', '--min-length',
              type=int,
              default=1,
              show_default=True,
              help='最小序列长度过滤|Minimum sequence length cutoff')
@click.option('-s',
              is_flag=True,
              help='Grep友好输出格式|Print grep-friendly output')
@click.option('-t',
              is_flag=True,
              help='Tab分隔输出|Print tab-delimited output')
@click.option('-u',
              is_flag=True,
              help='Tab分隔输出且无header|Print tab-delimited output without header')
@click.option('-o', '--output-dir',
              default='./assembly_stats_output',
              show_default=True,
              help='输出目录|Output directory')
def assembly_stats(input, min_length, s, t, u, output_dir):
    """
    基因组装配序列长度统计工具|Genome Assembly Sequence Length Statistics Tool

    报告FASTA和FASTQ文件的序列长度统计信息|Reports sequence length statistics from fasta and/or fastq files

    示例|Examples: biopytools assembly-stats -i genome.fa
    """

    # 延迟加载|Lazy loading
    assembly_stats_main = _lazy_import_assembly_stats_main()

    # 构建参数列表|Build argument list
    args = ['assembly_stats.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])

    # 可选参数|Optional parameters
    if min_length != 1:
        args.extend(['-l', str(min_length)])

    if output_dir != './assembly_stats_output':
        args.extend(['-o', output_dir])

    # 布尔选项|Boolean options
    if s:
        args.append('-s')

    if t:
        args.append('-t')

    if u:
        args.append('-u')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        assembly_stats_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
