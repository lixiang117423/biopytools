"""
AGP转表格命令|AGP to Table Converter Command
"""

import click
import sys
import os


def _lazy_import_agp2table_main():
    """延迟加载agp2table主函数|Lazy load agp2table main function"""
    try:
        from ...agp2table.main import main as agp2table_main
        return agp2table_main
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
    short_help='AGP转表格工具|AGP to table converter',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='AGP文件路径|AGP file path')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出表格文件路径|Output table file path')
@click.option('--format', '-f',
              default='txt',
              show_default=True,
              type=click.Choice(['txt', 'tsv', 'csv', 'xlsx']),
              help='输出格式|Output format')
@click.option('--statistics',
              is_flag=True,
              help='添加统计信息|Add statistics')
@click.option('--no-headers',
              is_flag=True,
              help='不添加表头|Do not add headers')
@click.option('--no-grouping',
              is_flag=True,
              help='不按scaffold分组|Do not group by scaffold')
def agp2table(input, output, format, statistics, no_headers, no_grouping):
    """
    AGP转表格工具|AGP to Table Converter

    将AGP格式文件转换为易读的表格格式，支持统计信息|Convert AGP format files to readable table format with statistics

    示例|Examples: biopytools agp2table -i assembly.agp -o assembly_table.txt
    """

    # 延迟加载|Lazy loading
    agp2table_main = _lazy_import_agp2table_main()

    # 构建参数列表|Build argument list
    args = ['agp2table.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 可选参数|Optional parameters
    if format != 'txt':
        args.extend(['-f', format])

    # 布尔选项|Boolean options
    if statistics:
        args.append('--statistics')

    if no_headers:
        args.append('--no-headers')

    if no_grouping:
        args.append('--no-grouping')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        agp2table_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
