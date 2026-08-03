"""
NeedleIdentity命令|Needle Identity Command
"""

import os
import sys

import click


def _lazy_import_needle_identity_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...needle_identity.main import main as needle_identity_main
        return needle_identity_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file exists (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='序列两两identity计算(EMBOSS needle)|Pairwise sequence identity (EMBOSS needle)',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120),
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA序列文件|Input FASTA file')
@click.option('--output-dir', '-o',
              default='./output',
              help='输出目录|Output directory (default: ./output)')
@click.option('--needle-path',
              default=None,
              help='needle可执行文件路径|needle executable path')
@click.option('--threads',
              type=int,
              default=12,
              help='并行线程数|Threads (default: 12)')
@click.option('--gapopen',
              type=float,
              default=10.0,
              help='gap开放罚分|Gap open penalty (default: 10.0)')
@click.option('--gapextend',
              type=float,
              default=0.5,
              help='gap延伸罚分|Gap extend penalty (default: 0.5)')
def needle_identity(input, output_dir, needle_path, threads, gapopen, gapextend):
    """
    序列两两identity计算(EMBOSS needle)|Pairwise sequence identity (EMBOSS needle)

    使用EMBOSS needle计算输入FASTA中所有序列对的两两identity|Compute pairwise identity of all sequence pairs in input FASTA using EMBOSS needle

    示例|Examples: biopytools needle-identity -i sequences.fa -o output_dir/
    """
    # 延迟加载|Lazy loading
    needle_identity_main = _lazy_import_needle_identity_main()

    # 构建参数列表|Build argument list
    args = ['needle_identity.py']
    args.extend(['-i', input])

    if output_dir != './output':
        args.extend(['-o', output_dir])
    if needle_path is not None:
        args.extend(['--needle-path', needle_path])
    if threads != 12:
        args.extend(['--threads', str(threads)])
    if gapopen != 10.0:
        args.extend(['--gapopen', str(gapopen)])
    if gapextend != 0.5:
        args.extend(['--gapextend', str(gapextend)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        needle_identity_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
