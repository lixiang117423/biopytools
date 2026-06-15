"""
VCF基因型统计命令|VCF Genotype Statistics Analysis Command
"""

import click
import sys


def _lazy_import_vcf_stats_main():
    """延迟加载vcf_stats_sample主函数|Lazy load vcf_stats_sample main function"""
    try:
        from ...vcf_sample_hete.main import main as vcf_stats_main
        return vcf_stats_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help='VCF基因型统计|VCF genotype statistics analysis',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              type=click.Path(exists=True),
              help='输入VCF文件路径|Input VCF file path')
@click.option('--output', '-o',
              default='vcf_stats_output',
              type=str,
              show_default=True,
              help='输出目录|Output directory')
@click.option('--min-depth', '-d',
              default=0,
              type=int,
              show_default=True,
              help='最小深度过滤阈值|Minimum depth filter threshold')
@click.option('--min-qual', '-q',
              default=0.0,
              type=float,
              show_default=True,
              help='最小质量分数过滤阈值|Minimum quality score filter threshold')
@click.option('--exclude-missing', '-e',
              is_flag=True,
              help='排除缺失基因型|Exclude missing genotypes')
@click.option('--no-detailed', '-D',
              is_flag=True,
              help='不输出详细统计结果|Do not output detailed statistics')
@click.option('--no-summary', '-S',
              is_flag=True,
              help='不输出汇总统计结果|Do not output summary statistics')
@click.option('--verbose',
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
def vcf_sample_hete(input, output, min_depth, min_qual, exclude_missing, no_detailed, no_summary,
                    verbose, quiet, log_file, log_level, dry_run):
    """
    VCF基因型统计分析工具|VCF Genotype Statistics Analysis Tool

    对VCF文件中的基因型进行统计分析|Calculate genotype statistics from VCF files

    示例|Examples: biopytools vcf-sample-hete -i variants.vcf -o vcf_stats_output
    """

    # 延迟加载|Lazy loading
    vcf_stats_main = _lazy_import_vcf_stats_main()

    # 构建参数列表|Build argument list
    args = ['vcf_sample_hete.py']

    # 必需参数|Required arguments
    args.extend(['-i', input])

    # 输出配置|Output configuration
    if output != 'vcf_stats_output':
        args.extend(['-o', output])

    # 过滤参数|Filtering parameters
    if min_depth != 0:
        args.extend(['-d', str(min_depth)])

    if min_qual != 0.0:
        args.extend(['-q', str(min_qual)])

    if exclude_missing:
        args.append('-e')

    # 输出控制|Output control
    if no_detailed:
        args.append('-D')

    if no_summary:
        args.append('-S')

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

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        vcf_stats_main()
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
