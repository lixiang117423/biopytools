"""
BAM覆盖度统计命令|BAM Coverage Statistics Command
"""

import click
import sys
import os


def _lazy_import_bam_cov_main():
    """延迟加载BAM覆盖度统计主函数|Lazy load BAM coverage statistics main function"""
    try:
        from ...bam_cov.main import main as bam_cov_main
        return bam_cov_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='BAM覆盖度统计|BAM coverage statistics',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入路径（BAM文件或包含BAM的目录）|Input path (BAM file or directory containing BAM files)')
@click.option('--chromosome', '-c',
              required=True,
              help='染色体名称|Chromosome name (e.g., chr1, Chr12)')
@click.option('--start', '-s',
              required=True,
              type=int,
              help='起始位置|Start position (1-based)')
@click.option('--end', '-e',
              type=int,
              help='终止位置|End position')
@click.option('--output-dir', '-o',
              default='./bam_coverage_stats_output',
              show_default=True,
              type=str,
              help='输出目录|Output directory')
@click.option('--output-prefix', '-p',
              default='coverage',
              show_default=True,
              type=str,
              help='输出文件前缀|Output file prefix')
@click.option('--min-mapq',
              type=int,
              default=0,
              show_default=True,
              help='最小mapping质量|Minimum mapping quality')
@click.option('--min-baseq',
              type=int,
              default=0,
              show_default=True,
              help='最小碱基质量|Minimum base quality')
@click.option('--no-merge',
              is_flag=True,
              help='不合并样本输出|Do not merge sample outputs')
@click.option('--no-summary',
              is_flag=True,
              help='不生成统计摘要|Do not generate summary statistics')
@click.option('--verbose', '-v',
              count=True,
              help='增加输出详细程度|Increase output verbosity')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode')
@click.option('--log-file',
              type=str,
              help='日志文件路径|Log file path')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
@click.option('--dry-run',
              is_flag=True,
              help='试运行模式|Dry run mode')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
def bam_cov(input, chromosome, start, end, output_dir, output_prefix,
            min_mapq, min_baseq, no_merge, no_summary, verbose, quiet,
            log_file, log_level, dry_run, threads):
    """
    BAM覆盖度统计工具|BAM Coverage Statistics Tool

    计算BAM文件中指定区域的reads覆盖度|Calculate per-base coverage for specified regions in BAM files

    示例|Examples: biopytools bam-cov -i sample.bam -c chr1 -s 1000000 -e 2000000
    """

    # 延迟加载|Lazy loading
    bam_cov_main = _lazy_import_bam_cov_main()

    # 构建参数列表|Build argument list
    args = ['bam_cov.py']

    # 必需参数|Required arguments
    args.extend(['-c', chromosome])
    args.extend(['-s', str(start)])
    args.extend(['-i', input])

    # 位置参数|Position parameters
    if end is not None:
        args.extend(['-e', str(end)])

    # 输出配置|Output configuration
    if output_dir != './bam_coverage_stats_output':
        args.extend(['-o', output_dir])

    if output_prefix != 'coverage':
        args.extend(['-p', output_prefix])

    # 过滤参数|Filtering parameters
    if min_mapq != 0:
        args.extend(['--min-mapq', str(min_mapq)])

    if min_baseq != 0:
        args.extend(['--min-baseq', str(min_baseq)])

    # 输出选项|Output options
    if no_merge:
        args.append('--no-merge')

    if no_summary:
        args.append('--no-summary')

    # 日志选项|Logging options
    if verbose > 0:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    if log_file:
        args.extend(['--log-file', log_file])

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    # 高级选项|Advanced options
    if dry_run:
        args.append('--dry-run')

    if threads != 1:
        args.extend(['-t', str(threads)])

    # 临时修改sys.argv|Temporarily modify sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        bam_cov_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断操作|Operation interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
