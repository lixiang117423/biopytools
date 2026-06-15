"""
PanDepth覆盖度计算命令|PanDepth Coverage Calculation Command
"""

import click
import sys
import os


def _lazy_import_pandepth_main():
    """延迟加载pandepth主函数|Lazy load pandepth main function"""
    try:
        from ...pandepth.main import main as pandepth_main
        return pandepth_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='PanDepth覆盖度计算工具|PanDepth coverage calculation tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入BAM文件或BAM文件目录|Input BAM file or BAM file directory')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--gff', '-g',
              type=str,
              help='GFF/GTF文件用于基因覆盖度|GFF/GTF file for gene coverage')
@click.option('--feature', '-f',
              type=click.Choice(['CDS', 'exon']),
              default='CDS',
              show_default=True,
              help='GFF/GTF特征类型|GFF/GTF feature type')
@click.option('--bed', '-b',
              type=str,
              help='BED文件用于特定区域覆盖度|BED file for specific region coverage')
@click.option('--window', '-w',
              type=int,
              help='滑动窗口大小(bp)|Sliding window size in bp')
@click.option('--min-mapq', '-q',
              type=int,
              default=0,
              show_default=True,
              help='最小比对质量|Minimum mapping quality')
@click.option('--min-depth', '-d',
              type=int,
              default=1,
              show_default=True,
              help='最小深度用于统计|Minimum depth for statistics')
@click.option('--exclude-flag', '-x',
              type=int,
              default=1796,
              show_default=True,
              help='排除reads的FLAG标志|FLAG bits to exclude reads')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='线程数|Number of threads')
@click.option('--reference', '-r',
              type=str,
              help='参考基因组文件|Reference genome file')
@click.option('--enable-gc', '-c',
              is_flag=True,
              help='启用GC含量计算|Enable GC content calculation')
@click.option('--all-sites', '-a',
              is_flag=True,
              help='输出所有位点深度|Output all site depths')
@click.option('--pandepth-path',
              type=str,
              default='~/software/PanDepth-2.26-Linux-x86_64/pandepth',
              show_default=True,
              help='PanDepth程序路径|PanDepth program path')
@click.option('--verbose', '-v',
              count=True,
              help='增加输出详细程度|Increase output verbosity')
@click.option('--quiet',
              is_flag=True,
              help='静默模式|Quiet mode')
@click.option('--log-file',
              type=str,
              help='日志文件路径|Log file path')
@click.option('--log-level',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
              default='INFO',
              show_default=True,
              help='日志级别|Log level')
def pandepth(input, output, gff, feature, bed, window, min_mapq, min_depth,
             exclude_flag, threads, reference, enable_gc, all_sites,
             pandepth_path, verbose, quiet, log_file, log_level):
    """
    PanDepth覆盖度计算工具|PanDepth Coverage Calculation Tool

    超快速基因组覆盖度计算工具，支持染色体、基因、特定区域和滑动窗口覆盖度分析|Ultra-fast genome coverage calculation tool, supporting chromosome, gene, specific region and sliding window coverage analysis

    示例|Examples: biopytools pandepth -i sample.bam -o coverage_results
    """

    # 延迟加载|Lazy load
    pandepth_main = _lazy_import_pandepth_main()

    # 构建参数列表|Build argument list
    args = ['pandepth.py']

    # 必需参数|Required parameters
    args.extend(['-i', input])
    args.extend(['-o', output])

    # 目标区域选项|Target region options
    if gff:
        args.extend(['-g', gff])

    if feature != 'CDS':
        args.extend(['-f', feature])

    if bed:
        args.extend(['-b', bed])

    if window:
        args.extend(['-w', str(window)])

    # 过滤选项|Filter options
    if min_mapq != 0:
        args.extend(['-q', str(min_mapq)])

    if min_depth != 1:
        args.extend(['-d', str(min_depth)])

    if exclude_flag != 1796:
        args.extend(['-x', str(exclude_flag)])

    # 其他选项|Other options
    if threads != 12:
        args.extend(['-t', str(threads)])

    if reference:
        args.extend(['-r', reference])

    if enable_gc:
        args.append('-c')

    if all_sites:
        args.append('-a')

    if pandepth_path != '~/software/PanDepth-2.26-Linux-x86_64/pandepth':
        args.extend(['--pandepth-path', pandepth_path])

    # 日志选项|Logging options
    if verbose > 0:
        args.extend(['-v'] * verbose)

    if quiet:
        args.append('--quiet')

    if log_file:
        args.extend(['--log-file', log_file])

    if log_level != 'INFO':
        args.extend(['--log-level', log_level])

    # 保存并恢复sys.argv|Save and restore sys.argv
    original_argv = sys.argv
    sys.argv = args

    try:
        # 调用主函数|Call main function
        pandepth_main()
    except SystemExit as e:
        # 处理系统退出|Handle system exit
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n用户中断|User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"运行时错误|Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
