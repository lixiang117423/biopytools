"""
Primer3引物设计|Primer3 Primer Design Command
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
def primer3(input_fasta, output_dir, primer_min_size, primer_opt_size, primer_max_size,
            primer_min_tm, primer_opt_tm, primer_max_tm,
            product_min_size, product_max_size, primer_num_return, output_format, output_header_lang,
            method, primer_end_margin, auto_product_size,
            product_size_min_ratio, product_size_max_ratio):
    """
    Primer3引物设计工具|Primer3 Primer Design Tool

    批量设计PCR引物，支持自定义引物长度、退火温度等参数|Design PCR primers in batch with customizable parameters

    示例|Example: biopytools primer3 -i sequences.fasta -o primer3_output
    """

    # 延迟加载|Lazy loading
    primer3_main = _lazy_import_primer3_main()

    # 构建参数列表|Build argument list
    args = ['primer3.py']

    # 必需参数|Required parameters
    args.extend(['-i', input_fasta])
    args.extend(['-o', output_dir])

    # 引物长度参数|Primer size parameters
    if primer_min_size != 20:
        args.extend(['--primer-min-size', str(primer_min_size)])
    if primer_opt_size != 20:
        args.extend(['--primer-opt-size', str(primer_opt_size)])
    if primer_max_size != 22:
        args.extend(['--primer-max-size', str(primer_max_size)])

    # 退火温度参数|Temperature parameters
    if primer_min_tm != 53.0:
        args.extend(['--primer-min-tm', str(primer_min_tm)])
    if primer_opt_tm != 58.0:
        args.extend(['--primer-opt-tm', str(primer_opt_tm)])
    if primer_max_tm != 63.0:
        args.extend(['--primer-max-tm', str(primer_max_tm)])

    # 产物大小参数|Product size parameters
    if product_min_size != 100:
        args.extend(['--product-min-size', str(product_min_size)])
    if product_max_size != 300:
        args.extend(['--product-max-size', str(product_max_size)])

    # 其他参数|Other parameters
    if primer_num_return != 5:
        args.extend(['--primer-num-return', str(primer_num_return)])
    if output_format != 'csv':
        args.extend(['--output-format', output_format])
    if output_header_lang != 'zh':
        args.extend(['--output-header-lang', output_header_lang])
    # method参数，默认为all，只有设置为random时才添加参数
    if method != 'all':
        args.extend(['--method', method])
    if primer_end_margin != 200:
        args.extend(['--primer-end-margin', str(primer_end_margin)])
    # auto_product_size默认为True，只有禁用时才添加参数
    if not auto_product_size:
        args.append('--no-auto-product-size')
    if product_size_min_ratio != 0.5:
        args.extend(['--product-size-min-ratio', str(product_size_min_ratio)])
    if product_size_max_ratio != 1.0:
        args.extend(['--product-size-max-ratio', str(product_size_max_ratio)])

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
