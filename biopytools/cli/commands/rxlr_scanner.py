"""
RxLR效应蛋白扫描命令|RxLR Effector Protein Scanner Command
"""

import click
import sys
import os


def _lazy_import_rxlr_main():
    """延迟加载rxlr_scanner主函数|Lazy load rxlr_scanner main function"""
    try:
        from ...rxlr_scanner.main import main as rxlr_main
        return rxlr_main
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
    short_help='RxLR效应蛋白扫描工具|RxLR effector protein scanner',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('-i', '--input',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件|Input FASTA file path')
@click.option('-o', '--output-prefix',
              required=True,
              help='输出文件前缀|Output file prefix')
@click.option('--window-start',
              type=int,
              default=20,
              show_default=True,
              help='窗口起始位置(20对应第21位)|Window start position (20 for position 21)')
@click.option('--window-end',
              type=int,
              default=120,
              show_default=True,
              help='窗口结束位置|Window end position')
@click.option('--min-length',
              type=int,
              default=120,
              show_default=True,
              help='最小序列长度|Minimum sequence length')
@click.option('--output-dir',
              default='./rxlr_scanner_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--no-excel',
              is_flag=True,
              help='不生成Excel输出|Do not generate Excel output')
@click.option('--no-tsv',
              is_flag=True,
              help='不生成TSV输出|Do not generate TSV output')
@click.option('-v', '--verbose',
              is_flag=True,
              help='详细输出模式|Verbose output mode')
@click.option('--log-file',
              is_flag=True,
              help='生成日志文件|Generate log file')
def rxlr_scanner(input, output_prefix, window_start, window_end, min_length,
                output_dir, no_excel, no_tsv, verbose, log_file):
    """
    RxLR效应蛋白扫描工具|RxLR Effector Protein Scanner

    在蛋白质序列中搜索RxLR和EER基序，用于鉴定卵菌效应蛋白|Search for RxLR and EER motifs in protein sequences for oomycete effector identification

    示例|Example: biopytools rxlr-scanner -i proteins.fa -o rxlr_results
    """

    # 延迟加载|Lazy loading
    rxlr_main = _lazy_import_rxlr_main()

    # 构建参数列表|Build argument list
    args = ['rxlr_scanner.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output_prefix])

    # 窗口参数|Window parameters
    if window_start != 20:
        args.extend(['--window-start', str(window_start)])

    if window_end != 120:
        args.extend(['--window-end', str(window_end)])

    if min_length != 120:
        args.extend(['--min-length', str(min_length)])

    # 输出选项|Output options
    if output_dir != './rxlr_scanner_output':
        args.extend(['--output-dir', output_dir])

    if no_excel:
        args.append('--no-excel')

    if no_tsv:
        args.append('--no-tsv')

    # 日志选项|Logging options
    if verbose:
        args.append('--verbose')

    if log_file:
        args.append('--log-file')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        rxlr_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
