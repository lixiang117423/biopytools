"""
FASTQ质控命令|FASTQ Quality Control Command
Optimized version: lazy loading for faster response
"""

import click
import sys
import os


def _lazy_import_fastp_main():
    """延迟加载fastp主函数|Lazy load fastp main function"""
    try:
        from ...fastp.main import main as fastp_main
        return fastp_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_input(input_path):
    """验证输入路径(目录或文件)|Validate input path (directory or file)"""
    if not _is_help_request():
        if not os.path.exists(input_path):
            raise click.BadParameter(f"输入路径不存在|Input path does not exist: {input_path}")
        if not (os.path.isdir(input_path) or os.path.isfile(input_path)):
            raise click.BadParameter(f"输入路径必须是目录或文件|Input path must be a directory or file: {input_path}")
    return input_path


def _validate_output_dir(dir_path):
    """验证输出目录路径(仅在非帮助模式)|Validate output directory path (only in non-help mode)"""
    if not _is_help_request():
        # 输出目录可由程序创建|Output directory can be created by program
        parent_dir = os.path.dirname(os.path.abspath(dir_path))
        if parent_dir and not os.path.exists(parent_dir):
            raise click.BadParameter(f"输出目录的父目录不存在|Parent directory of output does not exist: {parent_dir}")
    return dir_path


@click.command(
    short_help='FASTQ质控批处理工具|FASTQ Quality Control Batch Processing Tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_input(value) if value else None,
              help='输入原始FASTQ数据目录或文件|Input raw FASTQ data directory or file')
@click.option('--output-dir', '-o',
              required=True,
              callback=lambda ctx, param, value: _validate_output_dir(value) if value else None,
              help='输出清洁FASTQ数据目录|Output clean FASTQ data directory')
@click.option('--fastp-path',
              default='fastp',
              show_default=True,
              help='fastp可执行文件路径|fastp executable path')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--quality-threshold', '-q',
              type=int,
              default=30,
              show_default=True,
              help='质量阈值|Quality threshold')
@click.option('--min-length', '-l',
              type=int,
              default=50,
              show_default=True,
              help='最小长度|Minimum length')
@click.option('--unqualified-percent', '-u',
              type=int,
              default=40,
              show_default=True,
              help='不合格碱基百分比阈值|Unqualified base percentage threshold')
@click.option('--n-base-limit', '-n',
              type=int,
              default=10,
              show_default=True,
              help='N碱基数量限制|N base count limit')
@click.option('--read1-suffix',
              default='_1.fq.gz',
              show_default=True,
              help='Read1文件后缀（单末端模式也使用此参数）|Read1 file suffix (also used for single-end mode)')
@click.option('--read2-suffix',
              default='_2.fq.gz',
              show_default=True,
              help='Read2文件后缀|Read2 file suffix')
@click.option('--single-end',
              is_flag=True,
              help='单末端模式|Single-end mode')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet',
              is_flag=True,
              help='静默模式(仅输出ERROR)|Quiet mode (ERROR only)')
@click.option('--log-level',
              help='日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)|Log level')
@click.option('--log-file',
              help='日志文件路径|Log file path')
@click.option('--force', '-f',
              is_flag=True,
              help='强制覆盖已存在文件|Force overwrite existing files')
@click.option('--dry-run',
              is_flag=True,
              help='模拟运行(不实际执行)|Dry run without execution')
def fastp(input, output_dir, fastp_path, threads, quality_threshold,
          min_length, unqualified_percent, n_base_limit, read1_suffix, read2_suffix, single_end,
          verbose, quiet, log_level, log_file, force, dry_run):
    """
    FASTQ质控批处理工具|FASTQ Quality Control Batch Processing Tool

    使用fastp批量处理FASTQ文件质控|Batch quality control FASTQ files using fastp

    示例|Examples:
      # 处理整个目录|Process entire directory:
      biopytools fastp -i raw_data/ -o clean_data/

      # 处理单个文件|Process single file:
      biopytools fastp -i sample_1.fq.gz -o clean_data/
    """

    # 延迟加载|Lazy load: import only when actually called
    fastp_main = _lazy_import_fastp_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['fastp.py']  # 模拟脚本名|Simulate script name

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters (仅在非默认值时添加|add only when non-default)
    if fastp_path != 'fastp':
        args.extend(['--fastp-path', fastp_path])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if quality_threshold != 30:
        args.extend(['-q', str(quality_threshold)])

    if min_length != 50:
        args.extend(['-l', str(min_length)])

    if unqualified_percent != 40:
        args.extend(['-u', str(unqualified_percent)])

    if n_base_limit != 10:
        args.extend(['-n', str(n_base_limit)])

    if read1_suffix != '_1.fq.gz':
        args.extend(['--read1-suffix', read1_suffix])

    if read2_suffix != '_2.fq.gz':
        args.extend(['--read2-suffix', read2_suffix])

    if single_end:
        args.extend(['--single-end'])

    # 日志和执行控制参数|Logging and execution control parameters
    if verbose:
        args.extend(['-' + 'v' * verbose])

    if quiet:
        args.extend(['--quiet'])

    if log_level:
        args.extend(['--log-level', log_level])

    if log_file:
        args.extend(['--log-file', log_file])

    if force:
        args.extend(['--force'])

    if dry_run:
        args.extend(['--dry-run'])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        fastp_main()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
