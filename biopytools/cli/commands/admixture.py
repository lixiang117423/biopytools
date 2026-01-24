"""
ADMIXTURE群体结构分析命令|ADMIXTURE Population Structure Analysis Command
懒加载版本优化|Optimized version: lazy loading for faster response
"""

import click
import sys
import os


def _lazy_import_admixture_main():
    """延迟加载ADMIXTURE主函数|Lazy load ADMIXTURE main function"""
    try:
        from ...admixture.main import main as admixture_main
        return admixture_main
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
    short_help='ADMIXTURE群体结构分析工具|ADMIXTURE population structure analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--vcf', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入VCF文件路径|Input VCF file path')
@click.option('--output', '-o',
              default='admixture_results',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--min-k', '-k',
              type=int,
              default=2,
              show_default=True,
              help='最小K值|Minimum K value')
@click.option('--max-k', '-K',
              type=int,
              default=10,
              show_default=True,
              help='最大K值|Maximum K value')
@click.option('--cv-folds', '-c',
              type=int,
              default=5,
              show_default=True,
              help='交叉验证折数|Cross-validation folds')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--maf', '-m',
              type=float,
              default=0.05,
              show_default=True,
              help='最小等位基因频率阈值|MAF threshold')
@click.option('--missing', '-M',
              type=float,
              default=0.1,
              show_default=True,
              help='缺失率阈值|Missing rate threshold')
@click.option('--hwe', '-H',
              type=float,
              default=1e-6,
              show_default=True,
              help='HWE p值阈值|HWE p-value threshold')
@click.option('--skip-preprocessing', '-s',
              is_flag=True,
              help='跳过VCF预处理|Skip VCF preprocessing')
@click.option('--keep-intermediate',
              is_flag=True,
              help='保留中间文件|Keep intermediate files')
@click.option('--verbose', '-v',
              count=True,
              help='详细输出模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)')
@click.option('--quiet',
              is_flag=True,
              help='静默模式(仅ERROR)|Quiet mode (ERROR only)')
@click.option('--log-level',
              help='日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)|Log level')
@click.option('--log-file',
              help='日志文件路径|Log file path')
@click.option('--force', '-f',
              is_flag=True,
              help='强制覆盖已存在的文件|Force overwrite existing files')
@click.option('--dry-run',
              is_flag=True,
              help='试运行模式(不实际执行)|Dry run without execution')
def admixture(vcf, output, min_k, max_k, cv_folds, threads, maf, missing,
              hwe, skip_preprocessing, keep_intermediate,
              verbose, quiet, log_level, log_file, force, dry_run):
    """
    ADMIXTURE群体结构分析工具|ADMIXTURE Population Structure Analysis Tool

    使用ADMIXTURE进行群体结构分析，支持自动VCF预处理、质量控制、K值交叉验证|Use ADMIXTURE for population structure analysis with automatic VCF preprocessing, quality control, K-value cross-validation

    示例|Examples: biopytools admixture -i input.vcf -o results
    """

    # 延迟加载|Lazy loading
    admixture_main = _lazy_import_admixture_main()

    # 构建参数列表|Build argument list
    args = ['admixture.py']

    # 必需参数|Required parameters
    args.extend(['-i', vcf])

    # 可选参数|Optional parameters
    if output != 'admixture_results':
        args.extend(['-o', output])

    if min_k != 2:
        args.extend(['-k', str(min_k)])

    if max_k != 10:
        args.extend(['-K', str(max_k)])

    if cv_folds != 5:
        args.extend(['-c', str(cv_folds)])

    if threads != 64:
        args.extend(['-t', str(threads)])

    # 质控参数|Quality control parameters
    if maf != 0.05:
        args.extend(['-m', str(maf)])

    if missing != 0.1:
        args.extend(['-M', str(missing)])

    if hwe != 1e-6:
        args.extend(['-H', str(hwe)])

    # 布尔选项|Boolean options
    if skip_preprocessing:
        args.append('--skip-preprocessing')

    if keep_intermediate:
        args.append('--keep-intermediate')

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

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        admixture_main()
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
