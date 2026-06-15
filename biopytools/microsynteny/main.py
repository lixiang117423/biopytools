"""
微观共线性分析命令行入口|Microsynteny Analysis Command Line Interface
"""

import argparse
import sys
from pathlib import Path

from .config import MicrosyntenyConfig
from .analyzer import MicrosyntenyAnalyzer
from .utils import MicrosyntenyLogger


def parse_arguments():
    """解析命令行参数|Parse command line arguments

    Returns:
        argparse.Namespace: 解析后的参数|Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="微观共线性分析工具|Microsynteny Analysis Tool - Automated microsynteny analysis and visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i ./genome_data -g genes.txt -o output

步骤说明|Step information:
  Step 1: GFF/BED预处理|GFF/BED preprocessing (using JCVI)
  Step 2: 共线性分析|Synteny analysis (using JCVI)
  Step 3: 提取blocks|Extract synteny blocks
  Step 4: 可视化|Visualization (using pyCirclize)

依赖要求|Dependencies:
  - JCVI: 用于共线性分析|For synteny analysis (steps 1-3)
  - pyCirclize: 用于可视化|For visualization (step 4)
    安装|Install: pip install pycirclize
    或|or: conda install -c conda-forge pycirclize
        """
    )

    # 必需参数|Required parameters
    parser.add_argument(
        '-i', '--genome-folder',
        required=True,
        help='基因组文件夹路径|Genome folder path (containing A.fa, A.gff, B.fa, B.gff, ...)'
    )

    parser.add_argument(
        '-g', '--gene-list',
        required=True,
        help='目标基因列表文件|Target gene list file (two columns: species_id<tab>gene_id)'
    )

    # 可选参数|Optional parameters
    parser.add_argument(
        '-o', '--output-dir',
        default='./microsynteny_output',
        help='输出目录路径|Output directory path (default: ./microsynteny_output)'
    )

    parser.add_argument(
        '-j', '--jcvi-path',
        default='~/miniforge3/envs/jcvi_v.1.5.7',
        help='JCVI环境路径|JCVI environment path (default: ~/miniforge3/envs/jcvi_v.1.5.7)'
    )

    parser.add_argument(
        '--pycirclize-path',
        default='~/miniforge3/envs/pycirclize_v.1.10.1',
        help='pyCirclize环境路径|pyCirclize environment path (default: ~/miniforge3/envs/pycirclize_v.1.10.1)'
    )

    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=12,
        help='线程数|Number of threads (default: 12)'
    )

    parser.add_argument(
        '--extend-genes',
        type=int,
        default=30,
        help='延伸基因数|Number of genes to extend on each side (default: 30)'
    )

    parser.add_argument(
        '--cscore',
        type=float,
        default=0.99,
        help='共线性分数阈值(0-1)|Synteny score threshold 0-1 (default: 0.99)'
    )

    parser.add_argument(
        '--step',
        type=str,
        choices=['1', '2', '3', '4'],
        help='运行指定步骤|Run specific step only:\n'
             '1: 数据预处理|Data preprocessing\n'
             '2: 共线性分析|Synteny analysis\n'
             '3: 区块提取|Block extraction\n'
             '4: 绘图|Plotting'
    )

    parser.add_argument(
        '--log-level',
        type=str,
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help='日志级别|Log level (default: INFO)'
    )

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置对象|Create configuration object
        config = MicrosyntenyConfig(
            genome_folder=args.genome_folder,
            gene_list=args.gene_list,
            output_dir=args.output_dir,
            jcvi_path=args.jcvi_path,
            pycirclize_path=args.pycirclize_path,
            threads=args.threads,
            extend_genes=args.extend_genes,
            cscore=args.cscore,
            step=args.step
        )

        # 验证配置|Validate configuration
        config.validate()

        # 初始化日志器|Initialize logger
        logger_manager = MicrosyntenyLogger(
            log_file=config.log_file,
            log_level=args.log_level
        )
        logger = logger_manager.get_logger()

        # 打印配置信息|Print configuration info
        logger.info("=" * 60)
        logger.info("微观共线性分析工具|Microsynteny Analysis Tool")
        logger.info("=" * 60)
        logger.info(f"基因组文件夹|Genome folder: {config.genome_folder}")
        logger.info(f"基因列表文件|Gene list file: {config.gene_list}")
        logger.info(f"物种数量|Number of species: {len(config.get_species_list())}")
        logger.info(f"输出目录|Output directory: {config.output_dir}")
        logger.info(f"JCVI路径|JCVI path: {config.jcvi_path}")
        logger.info(f"pyCirclize路径|pyCirclize path: {config.pycirclize_path}")
        logger.info(f"线程数|Threads: {config.threads}")

        if args.step:
            logger.info(f"运行步骤|Run step: {args.step}")

        logger.info("=" * 60)

        # 创建分析器并运行|Create analyzer and run
        analyzer = MicrosyntenyAnalyzer(config, logger)
        analyzer.run_pipeline()

        logger.info("=" * 60)
        logger.info("所有步骤成功完成|All steps completed successfully!")
        logger.info(f"结果保存在|Results saved in: {config.output_dir}")
        logger.info("=" * 60)

        sys.exit(0)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)

    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
