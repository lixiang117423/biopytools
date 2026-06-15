"""
Protein Stats分析命令|Protein Stats Analysis Command
"""

import click
import sys
import os


def _lazy_import_protein_stats_main():
    """延迟加载protein_stats主函数|Lazy load protein_stats main function"""
    try:
        from ...protein_stats.main import main as protein_stats_main
        return protein_stats_main
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
    short_help='Protein Stats分析工具|Protein Stats analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--protein-fasta',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='蛋白序列FASTA文件|Protein sequence FASTA file')
@click.option('-o', '--output-file',
              default='protein_stats.tsv',
              show_default=True,
              help='输出文件路径|Output file path')
@click.option('--output-format',
              type=click.Choice(['tsv', 'csv', 'excel']),
              default='tsv',
              show_default=True,
              help='输出文件格式|Output file format')
@click.option('--no-length',
              is_flag=True,
              help='不计算序列长度|Do not calculate sequence length')
@click.option('--no-mw',
              is_flag=True,
              help='不计算分子量|Do not calculate molecular weight')
@click.option('--no-pi',
              is_flag=True,
              help='不计算等电点|Do not calculate isoelectric point')
@click.option('--aa-composition',
              is_flag=True,
              help='计算氨基酸组成|Calculate amino acid composition')
@click.option('--instability-index',
              is_flag=True,
              help='计算不稳定指数|Calculate instability index')
@click.option('--gravy',
              is_flag=True,
              help='计算脂肪指数(疏水性)|Calculate gravy (hydropathy)')
@click.option('--aromaticity',
              is_flag=True,
              help='计算芳香性|Calculate aromaticity')
def protein_stats(protein_fasta, output_file, output_format,
                   no_length, no_mw, no_pi, aa_composition,
                   instability_index, gravy, aromaticity):
    """
    Protein Stats分析工具：计算蛋白质序列理化性质|Protein Stats Analysis Tool: Calculate protein sequence properties

    示例|Examples: biopytools protein-stats -i proteins.fa
    """

    # 延迟加载|Lazy loading
    protein_stats_main = _lazy_import_protein_stats_main()

    # 构建参数列表|Build argument list
    args = ['protein_stats.py']

    # 必需参数|Required parameters
    args.extend(['-i', protein_fasta])

    # 输出配置|Output configuration
    if output_file != 'protein_stats.tsv':
        args.extend(['-o', output_file])

    if output_format != 'tsv':
        args.extend(['--output-format', output_format])

    # 布尔选项|Boolean options
    if no_length:
        args.append('--no-length')

    if no_mw:
        args.append('--no-mw')

    if no_pi:
        args.append('--no-pi')

    if aa_composition:
        args.append('--aa-composition')

    if instability_index:
        args.append('--instability-index')

    if gravy:
        args.append('--gravy')

    if aromaticity:
        args.append('--aromaticity')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        protein_stats_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
