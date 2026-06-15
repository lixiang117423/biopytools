"""
GFF3基因转录本信息提取命令|GFF3 Gene Transcript Information Extraction Command
"""

import click
import sys
import os


def _lazy_import_geneinfo_main():
    """延迟加载GFF3基因信息主函数|Lazy load GFF3 gene info main function"""
    try:
        from ...geneinfo.main import main as geneinfo_main
        return geneinfo_main
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


@click.command(short_help='GFF3基因转录本信息提取|GFF3 gene transcript info extraction',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入GFF3文件或目录路径|Input GFF3 file or directory path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出TSV文件路径|Output TSV file path')
@click.option('--gene-type',
              default='gene',
              show_default=True,
              help='基因特征类型|Gene feature type')
@click.option('--transcript-types',
              multiple=True,
              default=['mRNA', 'transcript'],
              show_default=True,
              help='转录本特征类型|Transcript feature types')
@click.option('--verbose', '-v',
              count=True,
              help='增加输出详细程度|Increase output verbosity')
@click.option('--quiet',
              is_flag=True,
              help='静默模式，仅输出错误信息|Quiet mode, only output errors')
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
              help='试运行模式，不实际执行|Dry run mode, no actual execution')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
def geneinfo(input, output, gene_type, transcript_types, verbose, quiet,
             log_file, log_level, dry_run, threads):
    """
    GFF3基因转录本信息提取工具|GFF3 Gene Transcript Information Extraction Tool

    从GFF3文件中为每个转录本提取整合的基因和转录本信息|Extract integrated gene and transcript information for each transcript from GFF3 files

    示例|Examples: biopytools geneinfo -i genome.gff3 -o gene_info.tsv
    """

    # 延迟加载|Lazy loading
    geneinfo_main = _lazy_import_geneinfo_main()

    # 构建参数列表|Build argument list
    args = ['geneinfo.py']

    # 必需参数|Required arguments
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数|Optional arguments
    if gene_type != 'gene':
        args.extend(['--gene-type', gene_type])

    # transcript-types - 检查是否为默认值
    # transcript-types - check if using default values
    default_transcript_types = ('mRNA', 'transcript')
    if transcript_types and tuple(transcript_types) != default_transcript_types:
        for tt in transcript_types:
            args.extend(['--transcript-types', tt])

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
        geneinfo_main()
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
