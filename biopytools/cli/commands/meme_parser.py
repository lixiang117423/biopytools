"""
MEME Parser分析命令|MEME Parser Analysis Command
"""

import click
import sys
import os


def _lazy_import_meme_parser_main():
    """延迟加载meme_parser主函数|Lazy load meme_parser main function"""
    try:
        from ...meme_parser.main import main as meme_parser_main
        return meme_parser_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_or_dir_exists(path):
    """验证文件或目录存在(仅在非帮助模式)|Validate file or directory exists (only in non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"文件或目录不存在|File or directory does not exist: {path}")
    return path


@click.command(
    short_help='MEME Parser分析工具|MEME Parser analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input-file',
              required=True,
              callback=lambda ctx, param, value: _validate_file_or_dir_exists(value) if value else None,
              help='输入FASTA文件或MEME输出文件(xml/txt)|Input FASTA file or MEME output file (xml/txt)')
@click.option('--parse-only',
              is_flag=True,
              help='仅解析模式，不运行MEME|Parse-only mode, do not run MEME')
@click.option('-o', '--output-prefix',
              default='meme_results',
              show_default=True,
              help='输出文件前缀|Output file prefix')
@click.option('--output-dir',
              default='.',
              show_default=True,
              help='输出目录|Output directory')
@click.option('--no-tsv',
              is_flag=True,
              help='不输出TSV文件|Do not output TSV file')
@click.option('--no-csv',
              is_flag=True,
              help='不输出CSV文件|Do not output CSV file')
@click.option('--no-excel',
              is_flag=True,
              help='不输出Excel文件|Do not output Excel file')
@click.option('--meme-path',
              default='~/miniforge3/envs/meme_v.5.5.9/bin/meme',
              show_default=True,
              help='MEME可执行文件路径|MEME executable path')
@click.option('-protein',
              is_flag=True,
              default=True,
              help='输入序列为蛋白质(默认)|Input sequences are protein (default)')
@click.option('-dna',
              is_flag=True,
              help='输入序列为DNA|Input sequences are DNA')
@click.option('-mod',
              default='zoops',
              show_default=True,
              type=click.Choice(['zoops', 'anr', 'oor']),
              help='Motif分布模式|Motif distribution mode')
@click.option('-nmotifs',
              default=10,
              show_default=True,
              type=int,
              help='Motif数量|Number of motifs')
@click.option('-minw',
              default=6,
              show_default=True,
              type=int,
              help='最小motif宽度|Minimum motif width')
@click.option('-maxw',
              default=50,
              show_default=True,
              type=int,
              help='最大motif宽度|Maximum motif width')
@click.option('-objfun',
              default='classic',
              show_default=True,
              type=click.Choice(['classic', 'de', 'ce', 'cd']),
              help='目标函数|Objective function')
@click.option('-markov-order',
              default=0,
              show_default=True,
              type=int,
              help='Markov链阶数|Markov chain order')
@click.option('--no-extract-motifs',
              is_flag=True,
              help='不提取motif序列|Do not extract motif sequences (default: extract)')
@click.option('--input-fasta',
              help='原始FASTA文件路径（解析模式时用于提取motif序列）|Original FASTA file path (for motif extraction in parse-only mode)')
def meme_parser(input_file, parse_only, output_prefix, output_dir,
                no_tsv, no_csv, no_excel, meme_path,
                protein, dna, mod, nmotifs, minw, maxw, objfun, markov_order,
                no_extract_motifs, input_fasta):
    """
    MEME Parser分析工具：运行MEME并解析结果|MEME Parser Analysis Tool: Run MEME and parse results

    示例|Examples: biopytools meme-parser -i proteins.fa -protein

    # 仅解析已有MEME输出
    biopytools meme-parser -i meme_out/meme.xml --parse-only
    """

    # 延迟加载|Lazy loading
    meme_parser_main = _lazy_import_meme_parser_main()

    # 处理序列类型参数互斥|Handle mutually exclusive sequence type parameters
    if dna:
        protein = False

    # 处理motif提取参数|Handle motif extraction parameter
    extract_motif_seqs = not no_extract_motifs

    # 构建参数列表|Build argument list
    args = ['meme_parser.py']

    # 必需参数|Required parameters
    args.extend(['-i', input_file])

    # 解析模式|Parse-only mode
    if parse_only:
        args.append('--parse-only')

    # 输出配置|Output configuration
    if output_prefix != 'meme_results':
        args.extend(['-o', output_prefix])

    if output_dir != '.':
        args.extend(['--output-dir', output_dir])

    # 布尔选项|Boolean options
    if no_tsv:
        args.append('--no-tsv')

    if no_csv:
        args.append('--no-csv')

    if no_excel:
        args.append('--no-excel')

    # MEME软件路径|MEME software path
    default_meme_path = '~/miniforge3/envs/meme_v.5.5.9/bin/meme'
    if meme_path != default_meme_path:
        args.extend(['--meme-path', meme_path])

    # MEME参数|MEME parameters
    # protein默认为True，所以总是添加（除非dna为True）
    if protein and not dna:
        args.append('-protein')

    if dna:
        args.append('-dna')

    if mod != 'zoops':
        args.extend(['-mod', mod])

    if nmotifs != 10:
        args.extend(['-nmotifs', str(nmotifs)])

    if minw != 6:
        args.extend(['-minw', str(minw)])

    if maxw != 50:
        args.extend(['-maxw', str(maxw)])

    if objfun != 'classic':
        args.extend(['-objfun', objfun])

    if markov_order != 0:
        args.extend(['-markov-order', str(markov_order)])

    # 序列提取选项|Sequence extraction options
    if not extract_motif_seqs:
        args.append('--no-extract-motifs')

    if input_fasta:
        args.extend(['--input-fasta', input_fasta])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        exit_code = meme_parser_main()
        sys.exit(exit_code)
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
