"""
Primer3引物设计主程序模块|Primer3 Primer Design Main Module
"""

import argparse
import sys
from .primer3_evaluator import Primer3Evaluator


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Primer3引物设计自动化脚本(模块化版本)|Primer3 Primer Design Automation Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input-fasta', required=True,
                       help='输入FASTA文件路径|Input FASTA file path')
    parser.add_argument('-o', '--output-dir', required=True,
                       help='输出目录路径|Output directory path')

    # Primer3路径配置|Primer3 path configuration
    parser.add_argument('--primer3-core-path',
                       default='~/miniforge3/envs/primer3_v.2.6.1/bin/primer3_core',
                       help='Primer3核心程序路径|Primer3 core program path')

    # 引物长度参数|Primer size parameters
    parser.add_argument('--primer-min-size', type=int, default=20,
                       help='最小引物长度|Minimum primer size')
    parser.add_argument('--primer-opt-size', type=int, default=20,
                       help='最优引物长度|Optimal primer size')
    parser.add_argument('--primer-max-size', type=int, default=22,
                       help='最大引物长度|Maximum primer size')

    # 退火温度参数|Annealing temperature parameters
    parser.add_argument('--primer-min-tm', type=float, default=53.0,
                       help='最小退火温度(°C)|Minimum annealing temperature (°C)')
    parser.add_argument('--primer-opt-tm', type=float, default=58.0,
                       help='最优退火温度(°C)|Optimal annealing temperature (°C)')
    parser.add_argument('--primer-max-tm', type=float, default=63.0,
                       help='最大退火温度(°C)|Maximum annealing temperature (°C)')

    # 产物大小范围|Product size range
    parser.add_argument('--product-min-size', type=int, default=100,
                       help='最小产物大小(bp)|Minimum product size (bp)')
    parser.add_argument('--product-max-size', type=int, default=300,
                       help='最大产物大小(bp)|Maximum product size (bp)')

    # 其他设计参数|Other design parameters
    parser.add_argument('--primer-num-return', type=int, default=5,
                       help='返回引物对数量|Number of primer pairs to return')
    parser.add_argument('--primer-max-ns', type=int, default=0,
                       help='允许的N碱基数量|Number of N bases accepted')
    parser.add_argument('--primer-gc-clamp', type=int, default=1,
                       help='GC clamp数量|GC clamp count')

    # 输出格式|Output format
    parser.add_argument('--output-format', choices=['csv', 'tsv', 'xlsx'],
                       default='csv',
                       help='输出文件格式|Output file format')
    parser.add_argument('--output-header-lang', choices=['zh', 'en'],
                       default='zh',
                       help='输出表头语言(zh:中文, en:英文)|Output header language (zh: Chinese, en: English)')
    # 引物设计策略|Primer design strategy
    parser.add_argument('--method', '-m', choices=['all', 'random'],
                       default='all',
                       help='引物设计策略: all=覆盖头尾(默认), random=随机设计|Primer design strategy: all=cover ends (default), random=random design')
    parser.add_argument('--primer-end-margin', type=int, default=200,
                       help='两端允许的引物位置范围bp,仅用于method=all|Allowed margin at ends in bp (only for method=all)')
    parser.add_argument('--auto-product-size', action='store_true', default=True,
                       help='自动根据序列长度设置产物大小范围(默认开启)|Auto set product size range based on sequence length (enabled by default)')
    parser.add_argument('--no-auto-product-size', dest='auto_product_size', action='store_false',
                       help='禁用自动产物大小范围|Disable automatic product size range')
    parser.add_argument('--product-size-min-ratio', type=float, default=0.5,
                       help='产物最小长度占序列长度的比例|Min product size ratio to sequence length (default: 0.5)')
    parser.add_argument('--product-size-max-ratio', type=float, default=1.0,
                       help='产物最大长度占序列长度的比例|Max product size ratio to sequence length (default: 1.0)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建评估器|Create evaluator
        evaluator = Primer3Evaluator(
            input_fasta=args.input_fasta,
            output_dir=args.output_dir,
            primer3_core_path=args.primer3_core_path,
            primer_min_size=args.primer_min_size,
            primer_opt_size=args.primer_opt_size,
            primer_max_size=args.primer_max_size,
            primer_min_tm=args.primer_min_tm,
            primer_opt_tm=args.primer_opt_tm,
            primer_max_tm=args.primer_max_tm,
            primer_product_size_range=(args.product_min_size, args.product_max_size),
            primer_num_return=args.primer_num_return,
            primer_max_ns_accepted=args.primer_max_ns,
            primer_gc_clamp=args.primer_gc_clamp,
            output_format=args.output_format,
            output_header_lang=args.output_header_lang,
            method=args.method,
            primer_end_margin=args.primer_end_margin,
            auto_product_size=args.auto_product_size,
            product_size_min_ratio=args.product_size_min_ratio,
            product_size_max_ratio=args.product_size_max_ratio
        )

        # 运行引物设计|Run primer design
        success = evaluator.run_design()

        if success:
            sys.exit(0)
        else:
            sys.exit(1)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
