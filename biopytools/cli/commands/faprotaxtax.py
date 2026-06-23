"""
FAPROTAX功能注释命令|FAPROTAX Functional Annotation Command
"""

import click
import sys
import os


def _lazy_import_faprotaxtax_main():
    """延迟加载faprotaxtax主函数|Lazy load faprotaxtax main function"""
    try:
        from ...faprotaxtax.main import main as faprotaxtax_main
        return faprotaxtax_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在（仅在非帮助模式）|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='FAPROTAX功能注释|FAPROTAX functional annotation',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input-table',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入OTU/ASV表（BIOM或TSV格式）|Input OTU/ASV table (BIOM or TSV format)')
@click.option('-o', '--output-dir',
              default='./faprotaxtax_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-g', '--groups-file',
              help='FAPROTAX功能组数据库文件路径|Path to FAPROTAX groups database file')
@click.option('--collapse-table-path',
              help='collapse_table.py脚本路径|Path to collapse_table.py script')
@click.option('--collapse-by-metadata',
              help='用于功能注释的BIOM元数据字段名|BIOM metadata field for functional annotation')
@click.option('--group-leftovers-as',
              help='未匹配到功能组的记录归为此组名|Group name for unmatched records')
@click.option('-n', '--normalize',
              type=click.Choice([
                  'none', 'columns_before_collapsing', 'rows_before_collapsing',
                  'columns_after_collapsing', 'rows_after_collapsing',
                  'columns_before_collapsing_excluding_unassigned',
                  'rows_before_collapsing_excluding_unassigned'
              ]),
              default='none',
              show_default=True,
              help='标准化方式|Normalization method')
@click.option('--average',
              type=click.Choice(['none', 'across_records', 'across_group_members',
                                 'across_used_group_members', 'maximum', 'minimum',
                                 'minimum_across_records']),
              default='none',
              show_default=True,
              help='组内聚合方式|Aggregation method within groups')
@click.option('--row-names-are-in-column',
              help='包含行名的列名或索引|Column name or index containing row names')
@click.option('--output-format',
              type=click.Choice(['auto', 'BIOM', 'classical']),
              default='auto',
              show_default=True,
              help='输出格式|Output format')
@click.option('-t', '--threads',
              type=int,
              default=1,
              show_default=True,
              help='线程数|Number of threads')
@click.option('-f', '--force',
              is_flag=True,
              help='覆盖已存在的输出文件|Overwrite existing output files')
@click.option('-v', '--verbose',
              is_flag=True,
              help='详细输出|Verbose output')
def faprotaxtax(input_table, output_dir, groups_file, collapse_table_path,
                collapse_by_metadata, group_leftovers_as, normalize, average,
                row_names_are_in_column, output_format, threads, force, verbose):
    """
    FAPROTAX微生物群落功能注释|FAPROTAX Microbial Community Functional Annotation

    基于分类学注释将微生物群落OTU/ASV表转换为功能丰度表|Convert taxonomic profiles to functional abundance profiles

    示例|Example: biopytools faprotaxtax -i otu_table.biom -o faprotaxtax_output
    """

    faprotaxtax_main = _lazy_import_faprotaxtax_main()

    args = ['faprotaxtax.py']
    args.extend(['-i', input_table])

    if output_dir != './faprotaxtax_output':
        args.extend(['-o', output_dir])

    if groups_file:
        args.extend(['-g', groups_file])

    if collapse_table_path:
        args.extend(['--collapse-table-path', collapse_table_path])

    if collapse_by_metadata:
        args.extend(['--collapse-by-metadata', collapse_by_metadata])

    if group_leftovers_as:
        args.extend(['--group-leftovers-as', group_leftovers_as])

    if normalize != 'none':
        args.extend(['-n', normalize])

    if average != 'none':
        args.extend(['--average', average])

    if row_names_are_in_column:
        args.extend(['--row-names-are-in-column', row_names_are_in_column])

    if output_format != 'auto':
        args.extend(['--output-format', output_format])

    if threads != 1:
        args.extend(['-t', str(threads)])

    if force:
        args.append('-f')

    if verbose:
        args.append('-v')

    original_argv = sys.argv
    sys.argv = args

    try:
        faprotaxtax_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
