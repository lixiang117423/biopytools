"""
Bismark甲基化分析命令|Bismark Methylation Analysis Command
"""

import click
import sys
import os


def _lazy_import_bismark_main():
    """延迟加载bismark主函数|Lazy load bismark main function"""
    try:
        from ...bismark.main import main as bismark_main
        return bismark_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _lazy_import_version():
    """延迟加载版本信息|Lazy load version info"""
    try:
        from ...bismark import __version__
        return __version__
    except ImportError:
        return "1.0.0"  # 默认版本|Default version


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help', '--version'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


def _validate_dir_exists(dir_path):
    """验证目录存在(仅在非帮助模式)|Validate directory existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.isdir(dir_path):
        raise click.BadParameter(f"目录不存在|Directory does not exist: {dir_path}")
    return dir_path


# 版本选项处理函数
def _get_version_option():
    """获取版本选项装饰器|Get version option decorator"""
    if _is_help_request():
        # 帮助模式
        return click.version_option(version="1.0.0", prog_name='Bismark Pipeline')
    else:
        # 正常模式
        version = _lazy_import_version()
        return click.version_option(version=version, prog_name='Bismark Pipeline')


@click.command(
    short_help='Bismark甲基化分析|Bismark methylation analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--version', is_flag=True, expose_value=False, is_eager=True,
              callback=lambda ctx, param, value: _handle_version_request() if value else None,
              help='显示版本信息|Show version information')
@click.option('--genome', '-g',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组FASTA文件|Genome FASTA file')
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_dir_exists(value) if value else None,
              help='原始FASTQ数据目录|Raw FASTQ data directory')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help='主输出目录|Main output directory')
@click.option('--pattern', '-p',
              default='_1_clean.fq.gz',
              show_default=True,
              help='R1文件后缀模式|Suffix pattern for R1 files')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--sort-buffer',
              type=str,
              default='400G',
              show_default=True,
              help='排序缓冲区大小|Sort buffer size')
@click.option('--include-overlap',
              is_flag=True,
              help='包含重叠读段|Include overlapping reads')
def bismark(genome, input, output_dir, pattern, threads, sort_buffer, include_overlap):
    """
    Bismark甲基化分析流程|Bismark methylation analysis pipeline

    Bismark全基因组重亚硫酸盐测序分析工具|WGBS analysis tool for DNA methylation detection

    示例|Examples: biopytools bismark -g genome.fa -i ./raw_data -o ./bismark_results
    """

    # 延迟加载: 仅在实际调用时导入
    bismark_main = _lazy_import_bismark_main()

    # 构建主函数的参数列表
    args = ['bismark.py']

    # 必需参数
    args.extend(['-g', genome])
    args.extend(['-i', input])
    args.extend(['-o', output_dir])

    # 可选参数
    if pattern != '_1_clean.fq.gz':
        args.extend(['-p', pattern])
    if threads != 88:
        args.extend(['-t', str(threads)])
    if sort_buffer != '400G':
        args.extend(['--sort-buffer', sort_buffer])

    # 处理反向布尔标志
    # argparse使用'--no-no-overlap' (action='store_false') 设置 no_overlap=True (默认排除重叠)
    # click使用--include-overlap 对应 --no-no-overlap 设置 no_overlap 为 False
    if include_overlap:
        args.append('--no-no-overlap')

    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数
        bismark_main()
    except SystemExit as e:
        # 处理正常程序退出
        if e.code != 0:
            click.secho(f"退出码|Exit code: {e.code}", fg='red', err=True)
        sys.exit(e.code)
    except Exception as e:
        click.secho(f"发生意外错误|An unexpected error occurred: {e}", fg='red', err=True)
        sys.exit(1)
    finally:
        # 无论结果如何都恢复原始sys.argv|Restore original sys.argv regardless of outcome
        sys.argv = original_argv


def _handle_version_request():
    """处理版本信息请求|Handle version information request"""
    version = _lazy_import_version()
    click.echo(f"Bismark Pipeline version {version}")
    sys.exit(0)


# 直接运行此文件时用于测试|If running this file directly for testing
if __name__ == '__main__':
    bismark()