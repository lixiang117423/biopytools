"""
BAM to FASTQ转换命令|BAM to FASTQ Conversion Command
"""

import click
import sys
import os


def _lazy_import_bam2fastq_main():
    """延迟加载BAM to FASTQ转换主函数|Lazy load BAM to FASTQ conversion main function"""
    try:
        from ...bam2fastq.main import main as bam2fastq_main
        return bam2fastq_main
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
    short_help='BAM to FASTQ批量转换工具|BAM to FASTQ batch conversion tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入BAM文件或文件夹路径|Input BAM file or directory path')
@click.option('-o', '--output-dir',
              required=True,
              type=click.Path(),
              help='输出文件夹路径|Output directory path')
@click.option('-t', '--threads',
              default=12,
              show_default=True,
              type=int,
              help='每个BAM文件转换使用的线程数|Threads per BAM file conversion')
@click.option('-j', '--jobs',
              default=1,
              show_default=True,
              type=int,
              help='并行处理的BAM文件数量|Number of parallel BAM file processing')
@click.option('--bam2fastq-path',
              default='bam2fastq',
              show_default=True,
              type=str,
              help='bam2fastq可执行文件路径|bam2fastq executable path')
def bam2fastq(input, output_dir, threads, jobs, bam2fastq_path):
    """
    BAM to FASTQ批量转换工具|BAM to FASTQ Batch Conversion Tool

    使用bam2fastq将BAM文件批量转换为FASTQ格式|Batch convert BAM files to FASTQ format using bam2fastq

    示例|Examples: biopytools bam2fastq -i ./sample.bam -o ./fastq_dir
    """

    # 延迟加载|Lazy loading
    bam2fastq_main = _lazy_import_bam2fastq_main()

    # 构建参数列表|Build argument list
    args = ['bam2fastq.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters
    if threads != 64:
        args.extend(['-t', str(threads)])

    if jobs != 1:
        args.extend(['-j', str(jobs)])

    if bam2fastq_path != 'bam2fastq':
        args.extend(['--bam2fastq-path', bam2fastq_path])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        bam2fastq_main()
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
