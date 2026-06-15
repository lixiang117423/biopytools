"""
OcBSA - 命令行入口|OcBSA - Command Line Entry
"""

import argparse
import sys
from pathlib import Path

from .utils import OcbsaLogger


def _run_f1(args):
    """运行F1群体OcBSA分析|Run F1 population OcBSA analysis"""
    from .config import OcbsaConfig
    from .core import OcbsaCalculator

    config = OcbsaConfig(
        input_vcf=args.input_vcf,
        parent1=args.p1,
        parent2=args.p2,
        pool1=args.b1,
        pool2=args.b2,
        output_dir=args.output_dir,
        window_size=args.window_size,
        pvalue=args.pvalue,
        parent_min_dep=args.parent_min_dep,
        parent_max_dep=args.parent_max_dep,
        pool_min_dep=args.pool_min_dep,
        pool_max_dep=args.pool_max_dep,
    )
    config.validate()

    log_file = str(config.output_path / "99_logs" / "ocbsa_f1.log")
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    logger_mgr = OcbsaLogger(log_file=log_file)
    logger = logger_mgr.get_logger()

    calculator = OcbsaCalculator(config, logger)
    success = calculator.run()
    sys.exit(0 if success else 1)


def _run_f2(args):
    """运行F2群体BSA分析|Run F2 population BSA analysis"""
    from .config import OcbsaConfig
    from .f2 import F2bsaCalculator

    config = OcbsaConfig(
        input_vcf=args.input_vcf,
        parent1=args.p1,
        parent2=args.p2,
        pool1=args.b1,
        pool2=args.b2,
        output_dir=args.output_dir,
        method=args.method,
        window_size=args.window_size,
        parent_min_dep=args.parent_min_dep,
        parent_max_dep=args.parent_max_dep,
        pool_min_dep=args.pool_min_dep,
        pool_max_dep=args.pool_max_dep,
    )
    config.validate()

    log_file = str(config.output_path / "99_logs" / "ocbsa_f2.log")
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    logger_mgr = OcbsaLogger(log_file=log_file)
    logger = logger_mgr.get_logger()

    calculator = F2bsaCalculator(config, logger)
    success = calculator.run()
    sys.exit(0 if success else 1)


def _run_fig(args):
    """运行BSA绘图|Run BSA figure plotting"""
    from .config import BsaFigConfig
    from .fig import BsaFigPlotter

    config = BsaFigConfig(
        input_file=args.input_file,
        output_file=args.output_file,
        plot_type=args.plot_type,
        position=args.position,
        color=args.color,
    )
    config.validate()

    logger_mgr = OcbsaLogger()
    logger = logger_mgr.get_logger()

    plotter = BsaFigPlotter(config, logger)
    success = plotter.run()
    sys.exit(0 if success else 1)


def _run_primer(args):
    """运行引物设计|Run primer design"""
    from .config import BsaPrimerConfig
    from .primer import BsaPrimerDesigner

    config = BsaPrimerConfig(
        genome=args.genome,
        ocvalue_file=args.ocvalue_file,
        region=args.region,
        output_dir=args.output_dir,
        primer_num=args.primer_num,
        flank_length=args.flank_length,
        primer_min_size=args.primer_min_size,
        primer_opt_size=args.primer_opt_size,
        primer_max_size=args.primer_max_size,
        product_min=args.product_min,
        product_max=args.product_max,
        min_tm=args.min_tm,
        max_tm=args.max_tm,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        tm_diff=args.tm_diff,
    )
    config.validate()

    log_file = str(config.output_path / "99_logs" / "ocbsa_primer.log")
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    logger_mgr = OcbsaLogger(log_file=log_file)
    logger = logger_mgr.get_logger()

    designer = BsaPrimerDesigner(config, logger)
    success = designer.run()
    sys.exit(0 if success else 1)


def main():
    """主入口函数|Main entry function"""
    parser = argparse.ArgumentParser(
        description='OcBSA - BSA分析工具套件|OcBSA - BSA Analysis Tool Suite',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest='subcommand', help='子命令|Subcommand')

    # f1 子命令|f1 subcommand
    p_f1 = subparsers.add_parser('f1', help='F1群体OcBSA DHHP分析|F1 population OcBSA DHHP analysis')
    p_f1.add_argument('-i', '--input-vcf', required=True, help='VCF文件路径|VCF file path')
    p_f1.add_argument('-p1', '--p1', type=int, required=True, help='显性亲本列号|Dominant parent column')
    p_f1.add_argument('-p2', '--p2', type=int, required=True, help='隐性亲本列号|Recessive parent column')
    p_f1.add_argument('-b1', '--b1', type=int, required=True, help='显性表型混池列号|Dominant pool column')
    p_f1.add_argument('-b2', '--b2', type=int, required=True, help='隐性表型混池列号|Recessive pool column')
    p_f1.add_argument('-o', '--output-dir', default='./output', help='输出目录|Output directory')
    p_f1.add_argument('-w', '--window-size', type=int, default=1000000, help='滑窗大小|Window size')
    p_f1.add_argument('-p', '--pvalue', type=float, default=99, help='p值阈值|P-value threshold')
    p_f1.add_argument('--parent-min-dep', type=int, default=10, help='亲本最低覆盖度|Parent min depth')
    p_f1.add_argument('--parent-max-dep', type=int, default=100, help='亲本最高覆盖度|Parent max depth')
    p_f1.add_argument('--pool-min-dep', type=int, default=20, help='混池最低覆盖度|Pool min depth')
    p_f1.add_argument('--pool-max-dep', type=int, default=500, help='混池最高覆盖度|Pool max depth')

    # f2 子命令|f2 subcommand
    p_f2 = subparsers.add_parser('f2', help='F2/RILs群体SNP-index/ED分析|F2/RILs population SNP-index/ED analysis')
    p_f2.add_argument('-i', '--input-vcf', required=True, help='VCF文件路径|VCF file path')
    p_f2.add_argument('-p1', '--p1', type=int, required=True, help='亲本1列号|Parent1 column')
    p_f2.add_argument('-p2', '--p2', type=int, required=True, help='亲本2列号|Parent2 column')
    p_f2.add_argument('-b1', '--b1', type=int, required=True, help='Pool1列号|Pool1 column')
    p_f2.add_argument('-b2', '--b2', type=int, required=True, help='Pool2列号|Pool2 column')
    p_f2.add_argument('-o', '--output-dir', default='./output', help='输出目录|Output directory')
    p_f2.add_argument('--method', default='snpindex', choices=['snpindex', 'ED'], help='分析方法|Analysis method')
    p_f2.add_argument('-w', '--window-size', type=int, default=1000000, help='滑窗大小|Window size')
    p_f2.add_argument('--parent-min-dep', type=int, default=10, help='亲本最低覆盖度|Parent min depth')
    p_f2.add_argument('--parent-max-dep', type=int, default=100, help='亲本最高覆盖度|Parent max depth')
    p_f2.add_argument('--pool-min-dep', type=int, default=20, help='混池最低覆盖度|Pool min depth')
    p_f2.add_argument('--pool-max-dep', type=int, default=500, help='混池最高覆盖度|Pool max depth')

    # fig 子命令|fig subcommand
    p_fig = subparsers.add_parser('fig', help='BSA结果可视化绘图|BSA result visualization')
    p_fig.add_argument('-i', '--input-file', required=True, help='输入结果文件|Input result file')
    p_fig.add_argument('-o', '--output-file', required=True, help='输出图片路径|Output figure path')
    p_fig.add_argument('--plot-type', default='ocvalue', choices=['ocvalue', 'snpindex', 'ed'],
                       help='图表类型|Plot type')
    p_fig.add_argument('--position', help='特定区域(chr,start,end)|Specific region')
    p_fig.add_argument('--color', default='plasma_r', help='颜色方案|Color scheme')

    # primer 子命令|primer subcommand
    p_primer = subparsers.add_parser('primer', help='BSA候选区域引物设计|BSA candidate region primer design')
    p_primer.add_argument('-g', '--genome', required=True, help='参考基因组路径|Reference genome path')
    p_primer.add_argument('-i', '--ocvalue-file', required=True, help='OcValue文件路径|OcValue file path')
    p_primer.add_argument('--region', required=True, help='目标区间(chr,start,end)|Target region')
    p_primer.add_argument('-o', '--output-dir', default='./output', help='输出目录|Output directory')
    p_primer.add_argument('-n', '--primer-num', type=int, default=10, help='引物对数量|Number of primer pairs')
    p_primer.add_argument('--flank-length', type=int, default=200, help='INDEL侧翼长度|INDEL flank length')
    p_primer.add_argument('--primer-min-size', type=int, default=18, help='最短引物|Min primer size')
    p_primer.add_argument('--primer-opt-size', type=int, default=20, help='最适引物|Opt primer size')
    p_primer.add_argument('--primer-max-size', type=int, default=24, help='最长引物|Max primer size')
    p_primer.add_argument('--product-min', type=int, default=70, help='最短产物|Min product size')
    p_primer.add_argument('--product-max', type=int, default=200, help='最长产物|Max product size')
    p_primer.add_argument('--min-tm', type=float, default=50.0, help='最低Tm|Min Tm')
    p_primer.add_argument('--max-tm', type=float, default=65.0, help='最高Tm|Max Tm')
    p_primer.add_argument('--min-gc', type=float, default=35.0, help='最低GC|Min GC%')
    p_primer.add_argument('--max-gc', type=float, default=65.0, help='最高GC|Max GC%')
    p_primer.add_argument('--tm-diff', type=float, default=0.5, help='Tm差异阈值|Tm diff threshold')

    args = parser.parse_args()

    if args.subcommand == 'f1':
        _run_f1(args)
    elif args.subcommand == 'f2':
        _run_f2(args)
    elif args.subcommand == 'fig':
        _run_fig(args)
    elif args.subcommand == 'primer':
        _run_primer(args)
    else:
        parser.print_help()
        sys.exit(0)
