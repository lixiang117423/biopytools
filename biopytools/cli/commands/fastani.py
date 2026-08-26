"""fastANI 全基因组ANI计算命令|fastANI whole-genome ANI command"""

import os
import sys

import click


def _lazy_import_fastani_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...fastani.main import main as fastani_main
        return fastani_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    return any(arg in {'-h', '--help'} for arg in sys.argv)


def _validate_path_exists(path):
    """验证路径存在(非帮助模式)|Validate path exists (non-help mode)"""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(
    short_help='全基因组ANI计算(FastANI)|Whole-genome ANI (FastANI)',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i', default=None,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              type=click.Path(),
              help='all-vs-all输入(目录/列表/单fasta;与-q/-r互斥)|'
                   'all-vs-all input (dir/list/single fasta)')
@click.option('--query', '-q', default=None,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              type=click.Path(),
              help='query侧输入(定向模式)|Query side (directional)')
@click.option('--reference', '-r', default=None,
              callback=lambda ctx, param, value: _validate_path_exists(value),
              type=click.Path(),
              help='reference侧输入(定向模式)|Reference side (directional)')
@click.option('--output-dir', '-o', default='./fastani_output', show_default=True,
              type=click.Path(), help='输出目录|Output directory')
@click.option('--threads', '-t', default=12, show_default=True, type=int,
              help='线程数|Thread count')
@click.option('--kmer', '-k', default=16, show_default=True, type=int,
              help='k-mer大小(<=16)|K-mer size (<=16)')
@click.option('--frag-len', default=3000, show_default=True, type=int,
              help='片段长度|Fragment length')
@click.option('--min-fraction', default=0.2, show_default=True, type=float,
              help='信任ANI的最小共享比例|Min shared fraction to trust ANI')
@click.option('--log-level', default='INFO', show_default=True,
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              help='日志级别|Log level')
def fastani(input, query, reference, output_dir, threads, kmer,
            frag_len, min_fraction, log_level):
    """fastANI全基因组ANI计算,输出矩阵与最近邻表
    |fastANI whole-genome ANI, producing matrix and nearest-neighbor table

    示例|Examples: biopytools fastani -i genome_dir/ -o output_dir/
    """

    fastani_main = _lazy_import_fastani_main()

    # 构造参数列表(默认值显式透传)|Build argv (defaults always forwarded)
    args = ['fastani']
    if input:
        args.extend(['-i', input])
    if query:
        args.extend(['-q', query])
    if reference:
        args.extend(['-r', reference])
    args.extend(['-o', output_dir, '-t', str(threads), '-k', str(kmer),
                 '--frag-len', str(frag_len),
                 '--min-fraction', str(min_fraction),
                 '--log-level', log_level])

    original_argv = sys.argv
    sys.argv = args
    try:
        fastani_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
