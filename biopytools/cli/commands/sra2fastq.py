"""
SRA转FASTQ命令|SRA to FASTQ Conversion Command
"""

import click
import sys
import os


def _lazy_import_sra_main():
    """延迟加载sra转换主函数|Lazy load SRA conversion main function"""
    try:
        from ...sra2fastq.main import main as sra_main
        return sra_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path existence (only in non-help mode)"""
    if _is_help_request():
        return path
    if path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='SRA转FASTQ工具|SRA to FASTQ conversion tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入SRA文件或文件夹路径|Input SRA file or folder path')
@click.option('--output', '-o',
              default='./fastq_output',
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--tmpdir',
              type=click.Path(),
              help='临时目录(用于加速)|Temporary directory (for acceleration)')
@click.option('--compress/--no-compress',
              default=True,
              show_default=True,
              help='压缩输出为.gz格式|Compress output to .gz format')
@click.option('--split/--no-split',
              'split_files',
              default=True,
              show_default=True,
              help='拆分双端测序文件|Split paired-end reads')
@click.option('--min-len',
              type=int,
              default=0,
              show_default=True,
              help='最小读长过滤|Minimum read length filter')
@click.option('--clip',
              is_flag=True,
              help='剪切adapters|Clip adapters')
def sra2fastq(input, output, threads, tmpdir, compress, split_files, min_len, clip):
    """
    SRA转FASTQ高速转换工具|SRA to FASTQ High-Speed Conversion Tool

    使用parallel-fastq-dump进行多线程高速转换，支持批量处理|High-speed SRA to FASTQ conversion using parallel-fastq-dump with multi-threading support and batch processing

    示例|Examples:
    biopytools sra2fastq -i sra_dir/ -o fastq_output
    """

    # 延迟加载|Lazy load
    sra_main = _lazy_import_sra_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['sra2fastq.py']
    args.extend(['-i', input])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if output != './fastq_output':
        args.extend(['-o', output])

    if threads != 88:
        args.extend(['-t', str(threads)])

    if tmpdir:
        args.extend(['--tmpdir', tmpdir])

    if not compress:
        args.append('--no-compress')

    if not split_files:
        args.append('--no-split')

    if min_len != 0:
        args.extend(['--min-len', str(min_len)])

    if clip:
        args.append('--clip')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        sra_main()
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
