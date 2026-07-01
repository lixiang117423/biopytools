"""基因密度计算主程序|Gene density calculation main program"""

import argparse
import sys
from pathlib import Path

from .config import GeneDensityConfig
from .utils import GeneDensityLogger
from .calculator import GeneDensityCalculator


def parse_arguments(args=None):
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='基因密度计算工具|Gene density calculation tool\n'
                    '按固定大小窗口统计每条染色体各区间的基因数量与基因密度(基因/Mb)|'
                    'Count genes and density (genes/Mb) per fixed-size window along each chromosome',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-i', '--gff',
                        required=True,
                        help='GFF3注释文件|GFF3 annotation file')
    parser.add_argument('-o', '--output-dir',
                        default='./gene_density_output',
                        help='输出目录|Output directory')
    parser.add_argument('-w', '--window-size',
                        type=int,
                        default=100000,
                        help='窗口大小(bp)|Window size (bp)')
    parser.add_argument('--feature-type',
                        default='gene',
                        help='统计的GFF feature类型(GFF第三列)|GFF feature type to count (column 3)')
    parser.add_argument('-g', '--genome',
                        default=None,
                        help='染色体长度来源(.fai或FASTA,可选,提升末尾窗口精度)|'
                        'Chromosome length source (.fai or FASTA, optional, improves last-window accuracy)')
    parser.add_argument('--prefix',
                        default=None,
                        help='输出文件前缀(默认GFF文件名stem)|Output file prefix (default: GFF stem)')
    parser.add_argument('--no-plot',
                        action='store_true',
                        help='不绘制密度图|Skip density plot')

    return parser.parse_args(args)


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建并校验配置|Create and validate config
        config = GeneDensityConfig(
            gff_file=args.gff,
            output_dir=args.output_dir,
            window_size=args.window_size,
            feature_type=args.feature_type,
            genome_file=args.genome,
            prefix=args.prefix,
            generate_plot=not args.no_plot,
        )
        config.validate()

        # 初始化日志|Initialize logger (99_logs/, §12)
        log_file = str(Path(config.output_dir) / "99_logs" / "gene_density.log")
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        logger_manager = GeneDensityLogger(log_file=log_file)
        logger = logger_manager.get_logger()

        # 运行计算|Run calculation
        calculator = GeneDensityCalculator(config, logger)
        calculator.run()

        sys.exit(0)

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}")
        sys.exit(1)

    except Exception as e:
        print(f"程序执行出错|Program execution error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
