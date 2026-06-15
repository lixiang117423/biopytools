"""
Resistify分析命令|Resistify Analysis Command
"""

import click
import sys
import os


def _lazy_import_resistify_main():
    """延迟加载resistify主函数|Lazy load resistify main function"""
    try:
        from ...resistify.main import main as resistify_main
        return resistify_main
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


@click.command(
    short_help='Resistify NLR分析工具|Resistify NLR analysis tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入蛋白质FASTA文件或目录|Input protein FASTA file or directory')
@click.option('-o', '--output-dir',
              default='./resistify_output',
              show_default=True,
              help='输出目录|Output directory')
@click.option('-t', '--threads',
              type=int,
              default=12,
              show_default=True,
              help='并行线程数|Number of parallel threads')
@click.option('--resistify-path',
              help='Resistify可执行文件路径|Resistify executable path')
@click.option('--skip-resistify',
              is_flag=True,
              help='跳过Resistify运行，仅解析已有结果|Skip Resistify run, only parse existing results')
@click.option('--output-prefix',
              default='resistify_results',
              show_default=True,
              help='解析结果文件前缀|Parsed results file prefix')
@click.option('--no-tsv',
              is_flag=True,
              help='不输出TSV文件|Do not output TSV file')
@click.option('--no-csv',
              is_flag=True,
              help='不输出CSV文件|Do not output CSV file')
@click.option('--no-excel',
              is_flag=True,
              help='不输出Excel文件|Do not output Excel file')
@click.option('--extract-nlr',
              is_flag=True,
              help='提取NLR序列|Extract NLR sequences')
@click.option('--extract-nbarc',
              is_flag=True,
              help='提取NB-ARC序列|Extract NB-ARC sequences')
@click.option('--filter-classification',
              help='按分类筛选(如TN, CNL, NL等)|Filter by classification (e.g., TN, CNL, NL)')
@click.option('--min-length',
              type=int,
              help='最小序列长度|Minimum sequence length')
@click.option('--max-length',
              type=int,
              help='最大序列长度|Maximum sequence length')
@click.option('--min-lrr-length',
              type=int,
              help='最小LRR长度|Minimum LRR length')
@click.option('--include-motifs',
              is_flag=True,
              help='包含motifs详情|Include motifs details')
def resistify(input, output_dir, threads, resistify_path, skip_resistify,
              output_prefix, no_tsv, no_csv, no_excel,
              extract_nlr, extract_nbarc, filter_classification,
              min_length, max_length, min_lrr_length, include_motifs):
    """
    Resistify NLR分析工具：运行Resistify并解析结果|Resistify NLR Analysis Tool: Run Resistify and parse results

    示例|Examples: biopytools resistify -i proteins.pep.fa -o ./resistify_output
    """

    # 延迟加载|Lazy loading
    resistify_main = _lazy_import_resistify_main()

    # 构建参数列表|Build argument list
    args = ['resistify.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output_dir])
    args.extend(['-t', str(threads)])

    # Resistify运行参数|Resistify run parameters
    if resistify_path:
        args.extend(['--resistify-path', resistify_path])
    if skip_resistify:
        args.append('--skip-resistify')

    # 输出配置|Output configuration
    if output_prefix != 'resistify_results':
        args.extend(['--output-prefix', output_prefix])

    # 布尔选项|Boolean options
    if no_tsv:
        args.append('--no-tsv')
    if no_csv:
        args.append('--no-csv')
    if no_excel:
        args.append('--no-excel')
    if extract_nlr:
        args.append('--extract-nlr')
    if extract_nbarc:
        args.append('--extract-nbarc')
    if include_motifs:
        args.append('--include-motifs')

    # 筛选选项|Filtering options
    if filter_classification:
        args.extend(['--filter-classification', filter_classification])
    if min_length is not None:
        args.extend(['--min-length', str(min_length)])
    if max_length is not None:
        args.extend(['--max-length', str(max_length)])
    if min_lrr_length is not None:
        args.extend(['--min-lrr-length', str(min_lrr_length)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        exit_code = resistify_main()
        sys.exit(exit_code)
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
