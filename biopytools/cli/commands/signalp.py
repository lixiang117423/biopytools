"""
SignalP 6.0信号肽预测命令|SignalP 6.0 Signal Peptide Prediction Command
"""

import click
import sys
import os


def _lazy_import_signalp_main():
    """延迟加载signalp主函数|Lazy load signalp main function"""
    try:
        from ...signalp.main import main as signalp_main
        return signalp_main
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
    short_help='SignalP 6.0信号肽预测工具|SignalP 6.0 signal peptide prediction tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件(氨基酸序列)|Input FASTA file (amino acid sequences)')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--organism', '-org',
              default='eukarya',
              type=click.Choice(['eukarya', 'other', 'euk']),
              help='生物类型|Organism type (default: eukarya)')
@click.option('--format', '-fmt',
              default='txt',
              type=click.Choice(['txt', 'png', 'eps', 'all', 'none']),
              help='输出格式|Output format (default: txt)')
@click.option('--mode', '-m',
              default='fast',
              type=click.Choice(['fast', 'slow', 'slow-sequential']),
              help='预测模式|Prediction mode (default: fast)')
@click.option('--bsize', '-bs',
              type=int,
              default=12,
              show_default=True,
              help='批处理大小|Batch size')
@click.option('--write-procs', '-wp',
              type=int,
              default=12,
              show_default=True,
              help='写入进程数|Number of write processes')
@click.option('--torch-num-threads', '-tt',
              type=int,
              default=12,
              show_default=True,
              help='PyTorch线程数|PyTorch threads')
@click.option('--signalp-path',
              default='~/miniforge3/envs/signalp6/bin/signalp6',
              show_default=True,
              help='SignalP程序路径|SignalP program path')
@click.option('--model-dir', '-md',
              default=None,
              help='模型权重目录|Model weights directory')
@click.option('--skip-resolve',
              is_flag=True,
              help='跳过冲突解析|Skip conflict resolution')
@click.option('--cleanup-plots',
              is_flag=True,
              default=True,
              show_default=True,
              help='自动删除plot文件|Auto-delete plot files')
@click.option('--keep-plots',
              is_flag=True,
              help='保留plot文件|Keep plot files (overrides --cleanup-plots)')
def signalp(input, output_dir, organism, format, mode, bsize, write_procs,
            torch_num_threads, signalp_path, model_dir, skip_resolve, cleanup_plots, keep_plots):
    """
    SignalP 6.0信号肽预测工具|SignalP 6.0 Signal Peptide Prediction Tool

    使用深度学习模型预测蛋白质信号肽|Predict signal peptides using deep learning models

    示例|Examples: biopytools signalp -i proteins.faa -o output_dir
    """

    # 延迟加载|Lazy loading
    signalp_main = _lazy_import_signalp_main()

    # 构建参数列表|Build argument list
    args = ['signalp.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output_dir])

    # 可选参数|Optional parameters
    if organism != 'eukarya':
        args.extend(['--organism', organism])

    if format != 'txt':
        args.extend(['--format', format])

    if mode != 'fast':
        args.extend(['--mode', mode])

    if bsize != 12:
        args.extend(['--bsize', str(bsize)])

    if write_procs != 12:
        args.extend(['--write-procs', str(write_procs)])

    if torch_num_threads != 12:
        args.extend(['--torch-num-threads', str(torch_num_threads)])

    if signalp_path != '~/miniforge3/envs/signalp6/bin/signalp6':
        args.extend(['--signalp-path', signalp_path])

    if model_dir:
        args.extend(['--model-dir', model_dir])

    if skip_resolve:
        args.append('--skip-resolve')

    # 处理plot清理参数|Handle plot cleanup parameters
    if keep_plots:
        args.append('--keep-plots')
    elif not cleanup_plots:
        args.append('--cleanup-plots')

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        signalp_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
