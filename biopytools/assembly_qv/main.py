"""
Merqury QV计算主程序模块|Merqury QV Calculation Main Module
"""

import argparse
import sys
from pathlib import Path

from .config import MerquryQVConfig
from .utils import MerquryQVLogger
from .qv_calculator import MerquryQVCalculator


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Merqury QV值计算工具|Merqury QV Value Calculation Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i fastq_dir/ -g genome.fa
        '''
    )

    parser.add_argument('-i', '--input',
                        required=True,
                        help='FASTQ文件或目录|FASTQ file or directory')

    parser.add_argument('-g', '--genome',
                        required=True,
                        help='基因组FASTA文件|Genome FASTA file')

    parser.add_argument('-o', '--output-dir',
                        default='./merqury_qv_output',
                        help='输出目录|Output directory')

    parser.add_argument('-k', '--kmer-size',
                        type=int,
                        default=None,
                        help='K-mer大小|K-mer size')

    parser.add_argument('-t', '--threads',
                        type=int,
                        default=24,
                        help='线程数|Number of threads')

    parser.add_argument('--conda-env',
                        default="~/miniforge3/envs/merqury_v.1.3/bin/",
                        help='Conda环境路径|Conda environment path')

    parser.add_argument('--data-type',
                        choices=['auto', 'illumina', 'hifi'],
                        default='auto',
                        help='数据类型|Data type')

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s 1.0.0')

    return parser.parse_args()


class MerquryQVRunner:
    """Merqury QV计算运行器|Merqury QV Calculation Runner"""

    def __init__(self, config: MerquryQVConfig, logger):
        self.config = config
        self.logger = logger
        self.calculator = MerquryQVCalculator(config, logger)

    def run(self):
        """运行分析|Run analysis"""
        self.logger.info("开始Merqury QV分析流程|Starting Merqury QV analysis pipeline")
        self.logger.info(f"输入FASTQ目录|Input FASTQ directory: {self.config.fastq_dir}")
        self.logger.info(f"基因组文件|Genome file: {self.config.genome_file}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

        # 运行完整分析|Run full analysis
        qv_stats = self.calculator.run_full_analysis()

        if qv_stats:
            # 打印摘要|Print summary
            summary = self.calculator.get_summary()
            print(summary)

            self.logger.info("分析完成|Analysis completed successfully")
            return 0
        else:
            self.logger.error("分析失败|Analysis failed")
            return 1


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = MerquryQVConfig(
            fastq_dir=args.input,
            genome_file=args.genome,
            output_dir=args.output_dir,
            kmer_size=args.kmer_size,
            threads=args.threads,
            conda_env=args.conda_env,
            data_type=args.data_type
        )
        config.validate()

        # 创建日志|Create logger
        log_file = Path(config.output_dir) / "merqury_qv.log"
        logger_manager = MerquryQVLogger(str(log_file))
        logger = logger_manager.get_logger()

        # 运行分析|Run analysis
        runner = MerquryQVRunner(config, logger)
        exit_code = runner.run()

        sys.exit(exit_code)

    except KeyboardInterrupt:
        print("用户中断|User interrupted")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
