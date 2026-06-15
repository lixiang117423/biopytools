"""
基因组挂载率统计主程序模块|Genome Mount Rate Main Module
"""

import argparse
import sys
from .config import GenomeMountRateConfig
from .utils import GenomeMountRateLogger
from .calculator import GenomeMountRateCalculator


class GenomeMountRateRunner:
    """基因组挂载率统计运行器|Genome Mount Rate Runner"""

    def __init__(self, **kwargs):
        """初始化运行器|Initialize runner

        Args:
            **kwargs: 配置参数|Configuration parameters
                - fasta_file: FASTA文件路径|FASTA file path
                - number: 序列数量|Number of sequences
                - sort_by_length: 是否按长度排序|Whether to sort by length
        """
        # 初始化配置|Initialize configuration
        self.config = GenomeMountRateConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = GenomeMountRateLogger()
        self.logger = self.logger_manager.get_logger()

        # 初始化计算器|Initialize calculator
        self.calculator = GenomeMountRateCalculator(self.config, self.logger)

    def run(self):
        """运行分析|Run analysis"""
        self.logger.info("开始计算基因组挂载率|Starting genome mount rate calculation")

        try:
            # 计算挂载率|Calculate mount rate
            results = self.calculator.calculate()

            # 输出结果|Output results
            self._print_results(results)

            self.logger.info("计算完成|Calculation completed")
            return results

        except Exception as e:
            self.logger.error(f"计算失败|Calculation failed: {str(e)}")
            raise

    def _print_results(self, results: dict):
        """打印结果|Print results

        Args:
            results: 计算结果字典|Calculation results dictionary
        """
        separator = "-" * 40

        self.logger.info(separator)
        self.logger.info(f"总序列数|Total sequences: {results['total_seqs']}")
        self.logger.info(f"总基因组大小|Total genome size: {self.calculator.format_number(results['total_bp'])} bp")
        self.logger.info(
            f"统计目标|Target: {'最长|longest' if results['sorted'] else '前|top'} {results['target_n']} 条序列|sequences"
        )
        self.logger.info(f"目标序列总长|Target sequences total length: {self.calculator.format_number(results['target_bp'])} bp")
        self.logger.info(separator)
        self.logger.info(f"占比|Mount rate: {results['percentage']:.2f}%")
        self.logger.info(separator)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="计算FASTA文件中前N条（或最长N条）序列占总基因组长度的百分比|"
                   "Calculate the percentage of top N (or longest N) sequences in total genome length",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--input', required=True,
                       help='输入的FASTA文件路径|Input FASTA file path')
    parser.add_argument('-n', '--number', type=int, required=True,
                       help='要计算的前N条序列数量|Number of top N sequences to calculate')

    # 可选参数|Optional parameters
    parser.add_argument('--sort', action='store_true',
                       help='[推荐]开启此选项后，会先按长度从大到小排序，再计算前N条（即计算最长N条的占比）|'
                            '[Recommended] Sort sequences by length (descending) before calculating (i.e., calculate longest N)')

    args = parser.parse_args()

    # 创建运行器并运行|Create runner and run
    try:
        runner = GenomeMountRateRunner(
            fasta_file=args.input,
            number=args.number,
            sort_by_length=args.sort
        )
        runner.run()
        sys.exit(0)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
