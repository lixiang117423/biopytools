"""
VCF基因型提取命令|VCF Genotype Extraction Command
"""

import click
import sys


def _lazy_import_vcf2genotype_main():
    """延迟加载vcf2genotype主函数|Lazy load vcf2genotype main function"""
    try:
        from ...vcf2genotype.main import main as vcf2genotype_main
        return vcf2genotype_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help='VCF基因型提取工具|VCF genotype extraction tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入VCF文件路径|Input VCF file path')
@click.option('--output', '-o',
              default='vcf2genotype',
              type=str,
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--samples', '-s',
              default='all',
              type=str,
              show_default=True,
              help='样本选择：all(所有样本)或逗号分隔的样本名称|Sample selection: all or comma-separated sample names')
@click.option('--biallelic-only',
              is_flag=True,
              help='只保留双等位位点|Keep only biallelic sites')
@click.option('--min-length',
              type=int,
              default=None,
              help='最小变异长度（包含），默认不限制|Minimum variant length (inclusive), default no limit')
@click.option('--max-length',
              type=int,
              default=None,
              help='最大变异长度（包含），默认不限制|Maximum variant length (inclusive), default no limit')
@click.option('--each', '-e',
              default='n',
              type=click.Choice(['yes', 'y', 'no', 'n']),
              show_default=True,
              help='按染色体拆分输出文件：yes/y或no/n|Split output files by chromosome: yes/y or no/n')
@click.option('--output-type', '-t',
              default='txt',
              type=click.Choice(['txt', 'csv']),
              show_default=True,
              help='输出文件格式|Output file format')
@click.option('--output-dir',
              default='./',
              type=click.Path(),
              show_default=True,
              help='输出目录|Output directory')
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
def vcf2genotype(input, output, samples, biallelic_only, min_length, max_length, each, output_type, output_dir,
                 verbose, quiet, log_file, log_level, dry_run):
    """
    VCF基因型提取工具|VCF Genotype Extraction Tool

    从VCF文件中提取基因型信息，支持多种输出格式|Extract genotype information from VCF files with multiple output formats

    示例|Examples: biopytools vcf2genotype -i variants.vcf -o genotypes
    """

    # 延迟加载|Lazy load
    vcf2genotype_main = _lazy_import_vcf2genotype_main()

    # 构建主函数参数列表|Build argument list for main function
    args = ['vcf2genotype.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])

    # 可选参数(仅在非默认时添加)|Optional parameters (add only when non-default)
    if output != 'vcf2genotype':
        args.extend(['-o', output])

    if samples != 'all':
        args.extend(['-s', samples])

    if output_type != 'txt':
        args.extend(['-t', output_type])

    if output_dir != './':
        args.extend(['--output-dir', output_dir])

    if biallelic_only:
        args.append('--biallelic-only')

    if min_length is not None:
        args.extend(['--min-length', str(min_length)])

    if max_length is not None:
        args.extend(['--max-length', str(max_length)])

    if each != 'n':
        args.extend(['-e', each])

    if dry_run:
        args.append('--dry-run')

    # 日志选项|Logging options
    if verbose > 0:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    if log_file:
        args.extend(['--log-file', log_file])

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        vcf2genotype_main()
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
