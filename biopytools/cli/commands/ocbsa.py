"""
OcBSA - BSA分析工具套件|OcBSA - BSA Analysis Tool Suite
"""

import click
import sys
import os


def _lazy_import_main():
    """延迟加载ocbsa主函数|Lazy load ocbsa main function"""
    try:
        from ...ocbsa.main import main as ocbsa_main
        return ocbsa_main
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


@click.group(
    short_help='BSA分析工具套件|BSA Analysis Tool Suite',
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120)
)
def ocbsa():
    """
    OcBSA - BSA分析工具套件|OcBSA - BSA Analysis Tool Suite

    基于DHHP算法的F1/F2群体BSA分析、结果可视化和引物设计

    示例|Examples: biopytools ocbsa f1 -i input.vcf -p1 1 -p2 2 -b1 3 -b2 4 -o ./output
    """
    pass


@ocbsa.command(short_help='F1群体OcBSA DHHP分析|F1 population OcBSA DHHP analysis')
@click.option('-i', '--input-vcf',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='VCF文件路径|VCF file path')
@click.option('-p1', type=int, required=True, help='显性亲本列号|Dominant parent column')
@click.option('-p2', type=int, required=True, help='隐性亲本列号|Recessive parent column')
@click.option('-b1', type=int, required=True, help='显性表型混池列号|Dominant pool column')
@click.option('-b2', type=int, required=True, help='隐性表型混池列号|Recessive pool column')
@click.option('-o', '--output-dir', default='./output', show_default=True, help='输出目录|Output directory')
@click.option('-w', '--window-size', type=int, default=1000000, show_default=True, help='滑窗大小|Window size')
@click.option('-p', '--pvalue', type=float, default=99, show_default=True, help='p值阈值|P-value threshold')
@click.option('--parent-min-dep', type=int, default=10, show_default=True, help='亲本最低覆盖度|Parent min depth')
@click.option('--parent-max-dep', type=int, default=100, show_default=True, help='亲本最高覆盖度|Parent max depth')
@click.option('--pool-min-dep', type=int, default=20, show_default=True, help='混池最低覆盖度|Pool min depth')
@click.option('--pool-max-dep', type=int, default=500, show_default=True, help='混池最高覆盖度|Pool max depth')
def f1(input_vcf, p1, p2, b1, b2, output_dir, window_size, pvalue,
        parent_min_dep, parent_max_dep, pool_min_dep, pool_max_dep):
    """
    F1群体OcBSA DHHP分析|F1 population OcBSA DHHP analysis

    示例|Examples: biopytools ocbsa f1 -i input.vcf -p1 1 -p2 2 -b1 3 -b2 4 -o ./output
    """
    ocbsa_main = _lazy_import_main()
    args = ['ocbsa.py', 'f1']
    args.extend(['-i', input_vcf, '-p1', str(p1), '-p2', str(p2), '-b1', str(b1), '-b2', str(b2)])
    args.extend(['-o', output_dir])
    if window_size != 1000000:
        args.extend(['--window-size', str(window_size)])
    if pvalue != 99:
        args.extend(['-p', str(pvalue)])
    if parent_min_dep != 10:
        args.extend(['--parent-min-dep', str(parent_min_dep)])
    if parent_max_dep != 100:
        args.extend(['--parent-max-dep', str(parent_max_dep)])
    if pool_min_dep != 20:
        args.extend(['--pool-min-dep', str(pool_min_dep)])
    if pool_max_dep != 500:
        args.extend(['--pool-max-dep', str(pool_max_dep)])

    original_argv = sys.argv
    sys.argv = args
    try:
        ocbsa_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@ocbsa.command(short_help='F2/RILs群体SNP-index/ED分析|F2/RILs population SNP-index/ED analysis')
@click.option('-i', '--input-vcf',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='VCF文件路径|VCF file path')
@click.option('-p1', type=int, required=True, help='亲本1列号|Parent1 column')
@click.option('-p2', type=int, required=True, help='亲本2列号|Parent2 column')
@click.option('-b1', type=int, required=True, help='Pool1列号|Pool1 column')
@click.option('-b2', type=int, required=True, help='Pool2列号|Pool2 column')
@click.option('-o', '--output-dir', default='./output', show_default=True, help='输出目录|Output directory')
@click.option('--method', default='snpindex', type=click.Choice(['snpindex', 'ED']),
              show_default=True, help='分析方法|Analysis method')
@click.option('-w', '--window-size', type=int, default=1000000, show_default=True, help='滑窗大小|Window size')
@click.option('--parent-min-dep', type=int, default=10, show_default=True, help='亲本最低覆盖度|Parent min depth')
@click.option('--parent-max-dep', type=int, default=100, show_default=True, help='亲本最高覆盖度|Parent max depth')
@click.option('--pool-min-dep', type=int, default=20, show_default=True, help='混池最低覆盖度|Pool min depth')
@click.option('--pool-max-dep', type=int, default=500, show_default=True, help='混池最高覆盖度|Pool max depth')
def f2(input_vcf, p1, p2, b1, b2, output_dir, method, window_size,
        parent_min_dep, parent_max_dep, pool_min_dep, pool_max_dep):
    """
    F2/RILs群体SNP-index/ED分析|F2/RILs population SNP-index/ED analysis

    示例|Examples: biopytools ocbsa f2 -i input.vcf -p1 1 -p2 2 -b1 3 -b2 4 --method snpindex -o ./output
    """
    ocbsa_main = _lazy_import_main()
    args = ['ocbsa.py', 'f2']
    args.extend(['-i', input_vcf, '-p1', str(p1), '-p2', str(p2), '-b1', str(b1), '-b2', str(b2)])
    args.extend(['-o', output_dir, '--method', method])
    if window_size != 1000000:
        args.extend(['--window-size', str(window_size)])
    if parent_min_dep != 10:
        args.extend(['--parent-min-dep', str(parent_min_dep)])
    if parent_max_dep != 100:
        args.extend(['--parent-max-dep', str(parent_max_dep)])
    if pool_min_dep != 20:
        args.extend(['--pool-min-dep', str(pool_min_dep)])
    if pool_max_dep != 500:
        args.extend(['--pool-max-dep', str(pool_max_dep)])

    original_argv = sys.argv
    sys.argv = args
    try:
        ocbsa_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@ocbsa.command(short_help='BSA结果可视化绘图|BSA result visualization')
@click.option('-i', '--input-file',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='输入结果文件|Input result file')
@click.option('-o', '--output-file', required=True, help='输出图片路径(.png/.pdf)|Output figure path')
@click.option('--plot-type', default='ocvalue', type=click.Choice(['ocvalue', 'snpindex', 'ed']),
              show_default=True, help='图表类型|Plot type')
@click.option('--position', help='特定区域(chr,start,end)|Specific region')
@click.option('--color', default='plasma_r', show_default=True, help='颜色方案|Color scheme')
def fig(input_file, output_file, plot_type, position, color):
    """
    BSA结果可视化绘图|BSA result visualization

    示例|Examples: biopytools ocbsa fig -i result.txt -o output.png --plot-type ocvalue
    """
    ocbsa_main = _lazy_import_main()
    args = ['ocbsa.py', 'fig']
    args.extend(['-i', input_file, '-o', output_file, '--plot-type', plot_type])
    if position:
        args.extend(['--position', position])
    if color != 'plasma_r':
        args.extend(['--color', color])

    original_argv = sys.argv
    sys.argv = args
    try:
        ocbsa_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv


@ocbsa.command(short_help='BSA候选区域引物设计|BSA candidate region primer design')
@click.option('-g', '--genome',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='参考基因组路径|Reference genome path')
@click.option('-i', '--ocvalue-file',
              required=True,
              callback=lambda ctx, param, value: _validate_path_exists(value) if value else None,
              help='OcValue文件路径|OcValue file path')
@click.option('--region', required=True, help='目标区间(chr,start,end)|Target region')
@click.option('-o', '--output-dir', default='./output', show_default=True, help='输出目录|Output directory')
@click.option('-n', '--primer-num', type=int, default=10, show_default=True, help='引物对数量|Number of primer pairs')
@click.option('--flank-length', type=int, default=200, show_default=True, help='INDEL侧翼长度|INDEL flank length')
@click.option('--primer-min-size', type=int, default=18, show_default=True, help='最短引物|Min primer size')
@click.option('--primer-opt-size', type=int, default=20, show_default=True, help='最适引物|Opt primer size')
@click.option('--primer-max-size', type=int, default=24, show_default=True, help='最长引物|Max primer size')
@click.option('--product-min', type=int, default=70, show_default=True, help='最短产物|Min product size')
@click.option('--product-max', type=int, default=200, show_default=True, help='最长产物|Max product size')
@click.option('--min-tm', type=float, default=50.0, show_default=True, help='最低Tm|Min Tm')
@click.option('--max-tm', type=float, default=65.0, show_default=True, help='最高Tm|Max Tm')
@click.option('--min-gc', type=float, default=35.0, show_default=True, help='最低GC|Min GC%')
@click.option('--max-gc', type=float, default=65.0, show_default=True, help='最高GC|Max GC%')
@click.option('--tm-diff', type=float, default=0.5, show_default=True, help='Tm差异阈值|Tm diff threshold')
def primer(genome, ocvalue_file, region, output_dir, primer_num, flank_length,
            primer_min_size, primer_opt_size, primer_max_size, product_min, product_max,
            min_tm, max_tm, min_gc, max_gc, tm_diff):
    """
    BSA候选区域引物设计|BSA candidate region primer design

    示例|Examples: biopytools ocbsa primer -g genome.fa -i result.OcValue --region Chr01,100000,500000 -o ./output
    """
    ocbsa_main = _lazy_import_main()
    args = ['ocbsa.py', 'primer']
    args.extend(['-g', genome, '-i', ocvalue_file, '--region', region, '-o', output_dir])
    if primer_num != 10:
        args.extend(['-n', str(primer_num)])
    if flank_length != 200:
        args.extend(['--flank-length', str(flank_length)])
    if primer_min_size != 18:
        args.extend(['--primer-min-size', str(primer_min_size)])
    if primer_opt_size != 20:
        args.extend(['--primer-opt-size', str(primer_opt_size)])
    if primer_max_size != 24:
        args.extend(['--primer-max-size', str(primer_max_size)])
    if product_min != 70:
        args.extend(['--product-min', str(product_min)])
    if product_max != 200:
        args.extend(['--product-max', str(product_max)])
    if min_tm != 50.0:
        args.extend(['--min-tm', str(min_tm)])
    if max_tm != 65.0:
        args.extend(['--max-tm', str(max_tm)])
    if min_gc != 35.0:
        args.extend(['--min-gc', str(min_gc)])
    if max_gc != 65.0:
        args.extend(['--max-gc', str(max_gc)])
    if tm_diff != 0.5:
        args.extend(['--tm-diff', str(tm_diff)])

    original_argv = sys.argv
    sys.argv = args
    try:
        ocbsa_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
