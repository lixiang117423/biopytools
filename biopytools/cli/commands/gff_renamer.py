"""
GFF重命名命令|GFF Renamer Command
"""

import click
import sys
import os


def _lazy_import_gff_renamer_main():
    """延迟加载gff_renamer主函数|Lazy load gff_renamer main function"""
    try:
        from ...gff_renamer.main import main as gff_renamer_main
        return gff_renamer_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(仅在非帮助模式)|Validate path exists (only in non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(short_help="GFF文件ID规范化|GFF File ID Standardization",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              help='输入GFF文件|Input GFF file')
@click.option('-o', '--output',
              required=True,
              help='输出GFF文件|Output GFF file')
@click.option('-p', '--prefix',
              required=True,
              help='ID前缀|ID prefix (e.g., CDRT, AGIS)')
@click.option('-s', '--species',
              required=True,
              help='物种缩写|Species abbreviation (e.g., Ov, Os)')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='并行线程数|Number of parallel threads')
@click.option('--output-mrna-mapping',
              is_flag=True,
              default=False,
              help='输出mRNA映射文件|Output mRNA mapping file')
@click.option('--mrna-mapping-file',
              default=None,
              help='mRNA映射文件路径（可选）|mRNA mapping file path (optional)')
@click.option('--chr-mapping',
              default=None,
              help='染色体映射文件路径|Chromosome mapping file path')
@click.option('--naming-format',
              type=click.Choice(['standard', 'simple', 'compact']),
              default='standard',
              show_default=True,
              help='命名格式|Naming format (standard/simple/compact)')
@click.option('--include-utr',
              is_flag=True,
              default=False,
              help='包含UTR特征|Include UTR features (five_prime_UTR, three_prime_UTR)')
@click.option('--skip-clean',
              is_flag=True,
              default=False,
              help='跳过AGAT清洗步骤|Skip AGAT GFF cleaning step')
def gff_renamer(input, output, prefix, species, threads, output_mrna_mapping, mrna_mapping_file,
                chr_mapping, naming_format, include_utr, skip_clean):
    """GFF文件ID规范化工具|GFF File ID Standardization Tool

    将GFF文件中的gene、mRNA、exon、CDS等特征的ID规范化为标准格式
    Standardize gene, mRNA, exon, CDS feature IDs in GFF file

    示例|Examples: biopytools gff-renamer -i input.gff -o output.gff -p CDRT -s Ov
    """

    # 延迟加载|Lazy load: import only when actually called
    gff_renamer_main_func = _lazy_import_gff_renamer_main()

    # 构建主函数的参数列表|Build argument list for main function
    args = ['gff_renamer.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])
    args.extend(['-p', prefix])
    args.extend(['-s', species])

    # 可选参数|Optional parameters
    args.extend(['-t', str(threads)])

    # mRNA映射文件参数|mRNA mapping file parameters
    if output_mrna_mapping:
        args.append('--output-mrna-mapping')
    if mrna_mapping_file:
        args.extend(['--mrna-mapping-file', mrna_mapping_file])

    # 新增参数|New parameters
    if chr_mapping:
        args.extend(['--chr-mapping', chr_mapping])
    args.extend(['--naming-format', naming_format])
    if include_utr:
        args.append('--include-utr')
    if skip_clean:
        args.append('--skip-clean')

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始主函数|Call original main function
        gff_renamer_main_func()
    except SystemExit as e:
        # 处理正常程序退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
