"""
Primer3引物设计|Primer3 Primer Design Command

透传约定|Forwarding convention: default 值也显式透传——模块侧 argparse/config
改默认值后 Click 层若不同步会静默失效, 总透传使两层默认值永远一致
|Defaults are always forwarded explicitly: module-side default changes can
otherwise be silently shadowed by stale Click defaults.
"""

import click
import sys
import os


def _lazy_import_primer3_main():
    """延迟加载primer3主函数|Lazy load primer3 main function"""
    try:
        from ...primer3.main import main as primer3_main
        return primer3_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if this is a help request"""
    help_flags = {'-h', '--help'}
    return any(arg in help_flags for arg in sys.argv)


def _validate_file_exists(file_path):
    """验证文件存在(仅在非帮助模式)|Validate file existence (only in non-help mode)"""
    if not _is_help_request() and not os.path.exists(file_path):
        raise click.BadParameter(f"文件不存在|File does not exist: {file_path}")
    return file_path


@click.command(
    short_help='Primer3引物设计工具|Primer3 primer design tool',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
@click.option('--input-fasta', '-i',
              required=True,
              callback=lambda ctx, param, value: _validate_file_exists(value) if value else None,
              help='输入FASTA文件|Input FASTA file path')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='输出目录|Output directory')
@click.option('--primer3-core-path',
              default='~/miniforge3/envs/misc/bin/primer3_core',
              show_default=True,
              help='Primer3核心程序路径|Primer3 core program path')
@click.option('--primer-min-size',
              type=int,
              default=20,
              show_default=True,
              help='最小引物长度|Minimum primer size')
@click.option('--primer-opt-size',
              type=int,
              default=20,
              show_default=True,
              help='最优引物长度|Optimal primer size')
@click.option('--primer-max-size',
              type=int,
              default=22,
              show_default=True,
              help='最大引物长度|Maximum primer size')
@click.option('--primer-min-tm',
              type=float,
              default=53.0,
              show_default=True,
              help='最小退火温度(°C)|Minimum annealing temperature (°C)')
@click.option('--primer-opt-tm',
              type=float,
              default=58.0,
              show_default=True,
              help='最优退火温度(°C)|Optimal annealing temperature (°C)')
@click.option('--primer-max-tm',
              type=float,
              default=63.0,
              show_default=True,
              help='最大退火温度(°C)|Maximum annealing temperature (°C)')
@click.option('--product-min-size',
              type=int,
              default=100,
              show_default=True,
              help='最小产物大小(bp)|Minimum product size (bp)')
@click.option('--product-max-size',
              type=int,
              default=300,
              show_default=True,
              help='最大产物大小(bp)|Maximum product size (bp)')
@click.option('--primer-num-return',
              type=int,
              default=5,
              show_default=True,
              help='返回引物对数量|Number of primer pairs to return')
@click.option('--primer-max-ns',
              type=int,
              default=0,
              show_default=True,
              help='允许的N碱基数量|Number of N bases accepted')
@click.option('--primer-gc-clamp',
              type=int,
              default=1,
              show_default=True,
              help='GC clamp数量|GC clamp count')
@click.option('--output-format',
              type=click.Choice(['csv', 'tsv', 'xlsx']),
              default='csv',
              show_default=True,
              help='输出文件格式|Output file format')
@click.option('--output-header-lang',
              type=click.Choice(['zh', 'en']),
              default='zh',
              show_default=True,
              help='输出表头语言(zh:中文, en:英文)|Output header language (zh: Chinese, en: English)')
@click.option('--method', '-m',
              type=click.Choice(['all', 'random']),
              default='all',
              show_default=True,
              help='引物设计策略: all=覆盖头尾, random=随机设计|Primer design strategy: all=cover ends, random=random design')
@click.option('--primer-end-margin',
              type=int,
              default=200,
              show_default=True,
              help='两端允许的引物位置范围bp(仅用于method=all)|Allowed margin at ends in bp (only for method=all)')
@click.option('--auto-product-size/--no-auto-product-size',
              default=True,
              show_default=False,
              help='自动根据序列长度设置产物大小范围(默认开启)|Auto set product size range based on sequence length (enabled by default)')
@click.option('--product-size-min-ratio',
              type=float,
              default=0.5,
              show_default=True,
              help='产物最小长度占序列长度的比例|Min product size ratio to sequence length (default: 0.5)')
@click.option('--product-size-max-ratio',
              type=float,
              default=1.0,
              show_default=True,
              help='产物最大长度占序列长度的比例|Max product size ratio to sequence length (default: 1.0)')
@click.option('--threads', '-t',
              type=int,
              default=12,
              show_default=True,
              help='并行进程数, 序列数达到阈值时生效(primer3_core单线程, 并行为多进程)|Parallel process count, active when sequence count reaches threshold (primer3_core is single-threaded; parallelism is multi-process)')
@click.option('--parallel-threshold',
              type=int,
              default=500,
              show_default=True,
              help='触发并行的序列数阈值, 低于该值保持单进程|Sequence count threshold to trigger parallel running; below it a single process is used')
def primer3(input_fasta, output_dir, primer3_core_path,
            primer_min_size, primer_opt_size, primer_max_size,
            primer_min_tm, primer_opt_tm, primer_max_tm,
            product_min_size, product_max_size, primer_num_return,
            primer_max_ns, primer_gc_clamp, output_format, output_header_lang,
            method, primer_end_margin, auto_product_size,
            product_size_min_ratio, product_size_max_ratio,
            threads, parallel_threshold):
    """
    Primer3引物设计工具|Primer3 Primer Design Tool

    批量设计PCR引物，支持自定义引物长度、退火温度等参数|Design PCR primers in batch with customizable parameters

    示例|Example: biopytools primer3 -i sequences.fasta -o primer3_output
    """

    # 延迟加载|Lazy loading
    primer3_main = _lazy_import_primer3_main()

    # 构建参数列表(default 总显式透传)|Build argument list (defaults always forwarded)
    args = ['primer3.py']

    # 必需参数|Required parameters
    args.extend(['-i', input_fasta])
    args.extend(['-o', output_dir])
    args.extend(['--primer3-core-path', primer3_core_path])

    # 引物长度参数|Primer size parameters
    args.extend(['--primer-min-size', str(primer_min_size)])
    args.extend(['--primer-opt-size', str(primer_opt_size)])
    args.extend(['--primer-max-size', str(primer_max_size)])

    # 退火温度参数|Temperature parameters
    args.extend(['--primer-min-tm', str(primer_min_tm)])
    args.extend(['--primer-opt-tm', str(primer_opt_tm)])
    args.extend(['--primer-max-tm', str(primer_max_tm)])

    # 产物大小参数|Product size parameters
    args.extend(['--product-min-size', str(product_min_size)])
    args.extend(['--product-max-size', str(product_max_size)])

    # 设计参数|Design parameters
    args.extend(['--primer-num-return', str(primer_num_return)])
    args.extend(['--primer-max-ns', str(primer_max_ns)])
    args.extend(['--primer-gc-clamp', str(primer_gc_clamp)])

    # 输出参数|Output parameters
    args.extend(['--output-format', output_format])
    args.extend(['--output-header-lang', output_header_lang])

    # 引物设计策略|Primer design strategy
    args.extend(['--method', method])
    args.extend(['--primer-end-margin', str(primer_end_margin)])

    # 自动产物大小(布尔开关)|Auto product size (boolean flag)
    if not auto_product_size:
        args.append('--no-auto-product-size')

    args.extend(['--product-size-min-ratio', str(product_size_min_ratio)])
    args.extend(['--product-size-max-ratio', str(product_size_max_ratio)])

    # 并行参数|Parallel parameters
    args.extend(['-t', str(threads)])
    args.extend(['--parallel-threshold', str(parallel_threshold)])

    # 执行主程序|Execute main program
    original_argv = sys.argv
    sys.argv = args

    try:
        primer3_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
