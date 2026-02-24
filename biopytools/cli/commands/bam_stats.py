"""
BAM文件批量统计分析命令|BAM File Batch Statistics Analysis Command
"""

import click
import sys
import os


def _lazy_import_bam_stats_main():
    """延迟加载bam_stats主函数|Lazy load bam_stats main function"""
    try:
        from ...bam_stats.main import main as bam_stats_main
        return bam_stats_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input_exists(input_path):
    """验证输入路径存在(仅在非帮助模式)|Validate input path exists (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(input_path):
        raise click.BadParameter(f"输入路径不存在|Input path does not exist: {input_path}")
    return input_path


def _validate_parent_dir_writable(file_path):
    """验证父目录可写(仅在非帮助模式)|Validate parent directory is writable (only in non-help mode)"""
    if not _is_help_request():
        parent_dir = os.path.dirname(os.path.abspath(file_path))
        if not os.path.exists(parent_dir):
            try:
                os.makedirs(parent_dir, exist_ok=True)
            except OSError:
                raise click.BadParameter(f"无法创建输出目录|Cannot create output directory: {parent_dir}")
        elif not os.access(parent_dir, os.W_OK):
            raise click.BadParameter(f"输出目录不可写|Output directory is not writable: {parent_dir}")
    return file_path


@click.command(
    short_help='BAM文件批量统计|BAM File Batch Statistics',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input-dir', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_input_exists(value) if value else None,
              help='BAM文件或目录|BAM file or directory containing BAM files')
@click.option('--output-file', '-o',
              required=True,
              callback=lambda ctx, param, value: _validate_parent_dir_writable(value) if value else None,
              help='输出报告文件|Output report file path (supports .xlsx and .csv formats)')
@click.option('--processes', '-p',
              type=int,
              help='并行进程数|Number of parallel processes')
@click.option('--log-dir',
              default='.',
              show_default=True,
              type=click.Path(),
              help='日志输出目录|Log output directory')
def bam_stats(input_dir, output_file, processes, log_dir):
    """
    BAM文件批量统计分析工具|BAM File Batch Statistics Analysis Tool

    批量处理BAM文件并生成统计报告|Batch process BAM files and generate statistics reports

    示例|Examples: biopytools bam-stats -i ./bam_files -o report.xlsx
                 biopytools bam-stats -i sample.bam -o report.xlsx
    """

    # 延迟加载|Lazy load: import only when actually called
    bam_stats_main = _lazy_import_bam_stats_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['bam_stats.py']

    # 必需参数|Required parameters
    args.extend(['--input-dir', input_dir])
    args.extend(['--output-file', output_file])

    # 可选参数(仅在非默认值时添加)|Optional parameters (add only when non-default)
    if processes is not None:
        args.extend(['--processes', str(processes)])

    if log_dir != '.':
        args.extend(['--log-dir', log_dir])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        bam_stats_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"运行时错误|Runtime Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
