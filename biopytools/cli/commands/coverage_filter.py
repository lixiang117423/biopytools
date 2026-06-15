"""
覆盖度过滤命令|Coverage Filter Command
"""

import click
import sys
import os


def _lazy_import_coverage_filter_main():
    """延迟加载coverage_filter主函数|Lazy load coverage_filter main function"""
    try:
        from ...coverage_filter.main import main as coverage_filter_main
        return coverage_filter_main
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


def _validate_parent_dir_writable(path):
    """验证父目录可写(仅在非帮助模式)|Validate parent directory is writable (only in non-help mode)"""
    if not _is_help_request():
        parent_dir = os.path.dirname(os.path.abspath(path))
        if not os.path.exists(parent_dir):
            try:
                os.makedirs(parent_dir, exist_ok=True)
            except OSError:
                raise click.BadParameter(f"无法创建输出目录|Cannot create output directory: {parent_dir}")
        elif not os.access(parent_dir, os.W_OK):
            raise click.BadParameter(f"输出目录不可写|Output directory is not writable: {parent_dir}")
    return path


@click.command(
    short_help='覆盖度过滤工具|Coverage Filter Tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--bam-file', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='BAM文件|BAM file path')
@click.option('--fasta-file', '-f',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组FASTA文件|Genome FASTA file')
@click.option('--output-prefix', '-o',
              required=True,
              help='输出文件前缀|Output file prefix')
@click.option('--output-dir', '-d',
              default='.',
              show_default=True,
              callback=lambda ctx, param, value: _validate_parent_dir_writable(value) if value else None,
              help='输出目录|Output directory')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--high-cov',
              type=float,
              default=90.0,
              show_default=True,
              help='高质量覆盖度阈值|High quality coverage threshold')
@click.option('--medium-cov-min',
              type=float,
              default=50.0,
              show_default=True,
              help='中等质量最小覆盖度|Medium quality minimum coverage')
def coverage_filter(bam_file, fasta_file, output_prefix, output_dir, threads, high_cov,
                    medium_cov_min):
    """
    覆盖度过滤工具|Coverage Filter Tool

    基于BAM覆盖度对序列进行质量分级和过滤|Sequence quality classification and filtering based on BAM coverage

    示例|Example: biopytools coverage-filter -i sample.bam -f genome.fa -o filtered
    """

    # 延迟加载|Lazy loading
    coverage_filter_main = _lazy_import_coverage_filter_main()

    # 构建参数列表|Build argument list
    args = ['coverage_filter.py']

    # 必需参数|Required parameters
    args.extend(['-i', bam_file])
    args.extend(['-f', fasta_file])
    args.extend(['-o', output_prefix])

    # 可选参数|Optional parameters (仅在非默认值时添加)
    if output_dir != '.':
        args.extend(['-d', output_dir])

    if threads != 12:
        args.extend(['-t', str(threads)])

    if high_cov != 90.0:
        args.extend(['--high-cov', str(high_cov)])

    if medium_cov_min != 50.0:
        args.extend(['--medium-cov-min', str(medium_cov_min)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        coverage_filter_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"运行时错误|Runtime Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
