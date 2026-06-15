"""
NCBI分类学注释CLI命令|NCBI Taxonomy Annotation CLI Command
"""

import click
import sys
import os


def _lazy_import_ncbi_taxo_main():
    """延迟加载ncbi_taxo主函数|Lazy load ncbi_taxo main function"""
    try:
        from ...ncbi_taxo.main import main as ncbi_taxo_main
        return ncbi_taxo_main
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


@click.command(
    short_help='NCBI分类学注释工具|NCBI Taxonomy Annotation Tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入文件（BLAST结果或accession列表）|Input file (BLAST results or accession list)')
@click.option('--output-prefix', '-o',
              required=True,
              help='输出文件前缀|Output file prefix')
@click.option('--input-type',
              type=click.Choice(['auto', 'blast', 'accession']),
              default='auto',
              show_default=True,
              help='输入文件类型|Input file type')
@click.option('--blast-column',
              type=int,
              default=2,
              show_default=True,
              help='BLAST结果中accession所在的列（从1开始）|Column containing accession in BLAST (1-based)')
@click.option('--min-length', '-l',
              type=int,
              default=1000,
              show_default=True,
              help='最小比对长度过滤（bp）|Minimum alignment length filter (bp)')
@click.option('--fetch-titles',
              is_flag=True,
              help='获取accession的序列描述（需要edirect）|Fetch accession titles (requires edirect)')
@click.option('--taxid-db',
              default='~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz',
              show_default=True,
              help='TaxID数据库路径|TaxID database path')
@click.option('--lineage-format',
              default='{k};{p};{c};{o};{f};{g};{s}',
              show_default=True,
              help='分类层级格式|Lineage format string')
@click.option('--no-full-lineage',
              is_flag=True,
              help='不保留完整lineage|Do not keep full lineage')
@click.option('--stats-by',
              multiple=True,
              type=click.Choice(['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']),
              default=['genus', 'species'],
              help='统计层级（可多选）|Statistics levels (can be specified multiple times)')
@click.option('--stats-target',
              type=click.Choice(['blast_hits', 'unique_accessions', 'both']),
              default='both',
              show_default=True,
              help='统计对象|Statistics target')
@click.option('--stats-output',
              type=click.Choice(['txt', 'csv']),
              default='txt',
              show_default=True,
              help='统计输出格式|Statistics output format')
@click.option('--threads', '-t',
              type=int,
              default=4,
              show_default=True,
              help='线程数|Number of threads')
def ncbi_taxo(input, output_prefix, input_type, blast_column, min_length, fetch_titles, taxid_db,
              lineage_format, no_full_lineage, stats_by, stats_target,
              stats_output, threads):
    """
    NCBI分类学注释工具|NCBI Taxonomy Annotation Tool

    从BLAST结果或accession列表获取NCBI分类学注释和统计信息|Get NCBI taxonomy annotations and statistics from BLAST results or accession lists.

    示例|Example: biopytools ncbi-taxo -i blast_results.txt -o output_prefix
    """

    # 延迟加载|Lazy loading
    ncbi_taxo_main = _lazy_import_ncbi_taxo_main()

    # 构建参数列表|Build argument list
    args = ['ncbi_taxo.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output_prefix])

    # 输入类型配置|Input type configuration
    if input_type != 'auto':
        args.extend(['--input-type', input_type])

    if blast_column != 2:
        args.extend(['--blast-column', str(blast_column)])

    if min_length != 1000:
        args.extend(['--min-length', str(min_length)])

    if fetch_titles:
        args.append('--fetch-titles')

    # 数据库配置|Database configuration
    if taxid_db != '~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz':
        args.extend(['--taxid-db', taxid_db])

    # 分类学配置|Taxonomy configuration
    if lineage_format != '{k};{p};{c};{o};{f};{g};{s}':
        args.extend(['--lineage-format', lineage_format])

    if no_full_lineage:
        args.append('--no-full-lineage')

    # 统计配置|Statistics configuration
    if stats_by != ('genus', 'species'):
        args.extend(['--stats-by'] + list(stats_by))

    if stats_target != 'both':
        args.extend(['--stats-target', stats_target])

    if stats_output != 'txt':
        args.extend(['--stats-output', stats_output])

    # 性能配置|Performance configuration
    if threads != 4:
        args.extend(['-t', str(threads)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        ncbi_taxo_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
