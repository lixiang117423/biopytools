"""
MSA可视化命令|MSA Visualization Command

"""

import click
import sys
import os


def _lazy_import_msaviz_main():
    """延迟加载msaviz主函数|Lazy load msaviz main function"""
    try:
        from ...msaviz.main import main as msaviz_main
        return msaviz_main
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
    short_help='MSA可视化工具（自动比对+可视化）|MSA Visualization Tool (Auto-align + Visualize)',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--infile', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入序列文件或MSA文件|Input sequences or MSA file')
@click.option('--outfile', '-o',
              required=True,
              help='输出可视化文件|Output visualization file (*.png|*.jpg|*.svg|*.pdf)')
# 比对参数|Alignment parameters
@click.option('--skip-align',
              is_flag=True,
              help='跳过MAFFT比对（输入已是比对结果）|Skip MAFFT alignment (input is already aligned)')
@click.option('--mafft-path',
              default='mafft',
              help='MAFFT可执行文件路径|MAFFT executable path (default: mafft)')
@click.option('--mafft-params',
              default='--auto',
              help='MAFFT参数|MAFFT parameters (default: --auto)')
@click.option('--threads',
              type=int,
              default=12,
              help='MAFFT线程数|MAFFT threads (default: 4)')
@click.option('--keep-alignment/--no-keep-alignment',
              default=True,
              help='保留比对结果文件|Keep alignment result file (default: True)')
# 格式参数|Format parameters
@click.option('--format',
              default='fasta',
              help='MSA文件格式|MSA file format (default: fasta)')
@click.option('--color-scheme',
              default='Zappo',
              help='颜色方案|Color scheme (default: Zappo)')
# 区域参数|Region parameters
@click.option('--start',
              type=int,
              default=1,
              help='起始位置(1-based)|Start position (1-based, default: 1)')
@click.option('--end',
              type=int,
              help='结束位置(1-based)|End position (1-based, default: MSA length)')
# 显示参数|Display parameters
@click.option('--wrap-length',
              type=int,
              default=60,
              help='换行长度|Wrap length (default: 60)')
@click.option('--wrap-space-size',
              type=float,
              default=3.0,
              help='换行间距|Wrap space size (default: 3.0)')
@click.option('--label-type',
              type=click.Choice(['id', 'description']),
              default='id',
              help='标签类型|Label type (default: id)')
@click.option('--show-label/--no-show-label',
              default=True,
              help='显示序列标签|Show sequence labels (default: True)')
@click.option('--show-grid',
              is_flag=True,
              help='显示网格|Show grid')
@click.option('--show-count',
              is_flag=True,
              help='显示字符统计|Show sequence character count')
@click.option('--show-consensus',
              is_flag=True,
              help='显示consensus序列|Show consensus sequence')
@click.option('--consensus-color',
              default='#1f77b4',
              help='Consensus颜色|Consensus color (default: #1f77b4)')
@click.option('--consensus-size',
              type=float,
              default=2.0,
              help='Consensus大小|Consensus size (default: 2.0)')
@click.option('--sort',
              is_flag=True,
              help='按NJ树排序|Sort by NJ tree')
@click.option('--dpi',
              type=int,
              default=300,
              help='图像DPI|Figure DPI (default: 300)')
def msaviz(infile, outfile, skip_align, mafft_path, mafft_params, threads, keep_alignment,
           format, color_scheme, start, end, wrap_length, wrap_space_size,
           show_label, label_type, show_grid, show_count, show_consensus,
           consensus_color, consensus_size, sort, dpi):
    """
    MSA可视化工具（自动比对+可视化）|MSA Visualization Tool (Auto-align + Visualize)

    多序列比对可视化工具，默认自动运行MAFFT比对|Multiple sequence alignment visualization tool with automatic MAFFT alignment

    示例|Examples: biopytools msaviz -i sequences.fa -o output.png
    """

    # 延迟加载|Lazy loading
    msaviz_main = _lazy_import_msaviz_main()

    # 构建参数列表|Build argument list
    args = ['msaviz.py']

    # 必需参数|Required parameters
    args.extend(['-i', infile])
    args.extend(['-o', outfile])

    # 比对参数|Alignment parameters
    if skip_align:
        args.append('--skip-align')

    if mafft_path != 'mafft':
        args.extend(['--mafft-path', mafft_path])

    if mafft_params != '--auto':
        args.extend(['--mafft-params', mafft_params])

    if threads != 4:
        args.extend(['--threads', str(threads)])

    if not keep_alignment:
        args.append('--no-keep-alignment')

    # 格式参数|Format parameters
    if format != 'fasta':
        args.extend(['--format', format])

    if color_scheme != 'Zappo':
        args.extend(['--color-scheme', color_scheme])

    # 区域参数|Region parameters
    if start != 1:
        args.extend(['--start', str(start)])

    if end is not None:
        args.extend(['--end', str(end)])

    # 显示参数|Display parameters
    if wrap_length != 60:
        args.extend(['--wrap-length', str(wrap_length)])

    if wrap_space_size != 3.0:
        args.extend(['--wrap-space-size', str(wrap_space_size)])

    if label_type != 'id':
        args.extend(['--label-type', label_type])

    if not show_label:
        args.append('--no-show-label')

    # 布尔选项|Boolean options
    if show_grid:
        args.append('--show-grid')

    if show_count:
        args.append('--show-count')

    if show_consensus:
        args.append('--show-consensus')

    if consensus_color != '#1f77b4':
        args.extend(['--consensus-color', consensus_color])

    if consensus_size != 2.0:
        args.extend(['--consensus-size', str(consensus_size)])

    if sort:
        args.append('--sort')

    # 输出参数|Output parameters
    if dpi != 300:
        args.extend(['--dpi', str(dpi)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        msaviz_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
