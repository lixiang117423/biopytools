"""基因信息+序列合并表 CLI 包装|gene-table CLI wrapper"""

import click
import sys
import os


def _lazy_import_gene_table_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...gene_table.main import main as gene_table_main
        return gene_table_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (non-help mode only)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='基因信息+序列合并表|Gene info + sequence merged table',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--genome', '-g', required=True,
              callback=lambda c, p, v: _validate_file_exists(v),
              help='基因组 FASTA|Genome FASTA')
@click.option('--gff', '-f', required=True,
              callback=lambda c, p, v: _validate_file_exists(v),
              help='GFF3 注释(支持 .gz)|GFF3 annotation (gz supported)')
@click.option('--output', '-o', required=True,
              help='输出表路径(或目录)|Output table path (or directory)')
@click.option('--prefix', default=None,
              help='输出前缀 + Sample 列(默认取 GFF 文件名)|Output prefix + Sample column')
@click.option('--longest-only', 'longest_only', is_flag=True,
              help='每基因仅保留最长转录本(默认全部)|Keep only longest transcript per gene')
@click.option('--min-length', type=int, default=0, show_default=True,
              help='基因 DNA 最小长度过滤|Min gene-DNA length filter')
@click.option('--gffread', default=None,
              help='gffread 路径(默认自动检测)|gffread path (auto-detected)')
@click.option('--log-file', default=None, help='日志文件|Log file')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO', show_default=True, help='日志级别|Log level')
@click.option('--verbose', '-v', is_flag=True, help='详细日志|Verbose')
def gene_table(genome, gff, output, prefix, longest_only, min_length,
               gffread, log_file, log_level, verbose):
    """
    基因信息+序列合并表|Gene info + sequence merged table

    输出一张每行一个转录本的 TSV(Gene_ID/Transcript_ID/坐标/Gene_DNA/CDS/Protein),
    并同时输出 gene.fa / cds.fa / pep.fa。|Outputs one row-per-transcript TSV plus gene.fa/cds.fa/pep.fa.

    示例|Examples: biopytools gene-table -g genome.fa -f input.gff -o out.tsv
    """

    gene_table_main = _lazy_import_gene_table_main()

    args = ['gene_table.py', '-g', genome, '-f', gff, '-o', output]
    if prefix:
        args.extend(['--prefix', prefix])
    if longest_only:
        args.append('--longest-only')
    if min_length != 0:
        args.extend(['--min-length', str(min_length)])
    if gffread:
        args.extend(['--gffread', gffread])
    if log_file:
        args.extend(['--log-file', log_file])
    if log_level != 'INFO':
        args.extend(['--log-level', log_level])
    if verbose:
        args.append('--verbose')

    original_argv = sys.argv
    sys.argv = args
    try:
        gene_table_main()
    except SystemExit as e:
        if e.code not in (0, None):
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|Interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
