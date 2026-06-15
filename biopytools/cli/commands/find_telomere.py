"""
端粒识别命令|Telomere Finder Command
"""

import click
import sys
import os


def _lazy_import_telomere_finder_main():
    """延迟加载telomere_finder主函数|Lazy load telomere finder main function"""
    try:
        from ...find_telomere.main import main as telomere_finder_main
        return telomere_finder_main
    except ImportError as e:
        click.echo(f"Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在性(仅在非帮助模式下)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and file_path and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(short_help="端粒识别工具|Telomere Finder",
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--genome', '-g',
              required=False,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='基因组文件|Genome sequence file (FASTA format)')
@click.option('--mode', '-m',
              type=click.Choice(['explore', 'find', 'search', 'plot', 'pipeline']),
              default='pipeline',
              show_default=True,
              help='分析模式|Analysis mode')
@click.option('--output-dir', '-o',
              default='./telomere_output',
              show_default=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--prefix', '-p',
              default='telomere',
              show_default=True,
              help='输出前缀|Output prefix')
@click.option('--clade', '-c',
              help='分类群名称|Clade name')
@click.option('--window', '-w',
              type=int,
              default=10000,
              show_default=True,
              help='窗口大小|Window size')
@click.option('--search-string', '-s',
              help='搜索字符串|Search string')
@click.option('--format',
              type=click.Choice(['tsv', 'bedgraph']),
              default='tsv',
              show_default=True,
              help='输出格式|Output format')
@click.option('--tsv', '-t',
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='TSV文件|TSV file')
@click.option('--plot-height',
              type=int,
              default=200,
              show_default=True,
              help='图像高度|Plot height')
@click.option('--plot-width',
              type=int,
              default=1000,
              show_default=True,
              help='图像宽度|Plot width')
@click.option('--tidk-path',
              default='~/miniforge3/envs/tidk_v.0.2.65/bin/tidk',
              show_default=True,
              help='tidk软件路径|tidk software path')
@click.option('--explore-min',
              type=int,
              default=5,
              show_default=True,
              help='最小重复长度|Min repeat length')
@click.option('--explore-max',
              type=int,
              default=12,
              show_default=True,
              help='最大重复长度|Max repeat length')
@click.option('--explore-threshold',
              type=int,
              default=100,
              show_default=True,
              help='重复阈值|Repeat threshold')
@click.option('--explore-distance',
              type=float,
              default=0.01,
              show_default=True,
              help='距离比例|Distance ratio')
@click.option('--verbose', '-v',
              is_flag=True,
              help='详细输出模式|Verbose output mode')
@click.option('--log-file',
              help='日志文件路径|Log file path')
@click.option('--print-clades',
              is_flag=True,
              help='打印支持的分类群列表|Print supported clades list')
def find_telomere(genome, mode, output_dir, prefix, clade, window, search_string,
                  format, tsv, plot_height, plot_width, tidk_path,
                  explore_min, explore_max, explore_threshold, explore_distance,
                  verbose, log_file, print_clades):
    """
    端粒识别分析工具|Telomere Finder Analysis Tool

    示例|Examples: biopytools find-telomere -g genome.fa -o results
    """

    # 打印分类群列表|Print clades list
    if print_clades:
        from ...find_telomere.utils import CladeDatabase
        CladeDatabase.print_clade_table()
        sys.exit(0)

    # 检查必需参数|Check required parameters
    if not genome:
        raise click.BadParameter('缺少必需选项"--genome"/"-g"|Missing required option "--genome"/"-g"')

    # 延迟加载|Lazy loading: import only when actually called
    telomere_finder_main = _lazy_import_telomere_finder_main()

    # 构建参数列表|Build argument list
    args = ['telomere_finder.py']

    # Required parameters
    args.extend(['-g', genome])
    args.extend(['--mode', mode])

    # Output parameters
    if output_dir != './telomere_output':
        args.extend(['--output-dir', output_dir])

    if prefix != 'telomere':
        args.extend(['--prefix', prefix])

    # 模式特定参数|Mode-specific parameters
    if mode == 'find' and clade:
        args.extend(['--clade', clade])

    if window != 10000:
        args.extend(['--window', str(window)])

    if mode == 'search' and search_string:
        args.extend(['--search-string', search_string])

    if mode == 'search' and format != 'tsv':
        args.extend(['--format', format])

    if mode == 'plot' and tsv:
        args.extend(['--tsv', tsv])

    if mode == 'plot' and plot_height != 200:
        args.extend(['--plot-height', str(plot_height)])

    if mode == 'plot' and plot_width != 1000:
        args.extend(['--plot-width', str(plot_width)])

    if mode == 'explore':
        if explore_min != 5:
            args.extend(['--explore-min', str(explore_min)])
        if explore_max != 12:
            args.extend(['--explore-max', str(explore_max)])
        if explore_threshold != 100:
            args.extend(['--explore-threshold', str(explore_threshold)])
        if explore_distance != 0.01:
            args.extend(['--explore-distance', str(explore_distance)])

    # 软件配置|Software configuration
    if tidk_path != '~/miniforge3/envs/tidk_v.0.2.65/bin/tidk':
        args.extend(['--tidk-path', tidk_path])

    # 日志参数|Logging parameters
    if verbose:
        args.append('--verbose')

    if log_file:
        args.extend(['--log-file', log_file])

    # 保存和恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用原始main函数|Call original main function
        telomere_finder_main()
    except SystemExit as e:
        # 处理程序正常退出|Handle normal program exit
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
