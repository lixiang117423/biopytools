"""
Hi-C数据质量控制评估主程序模块|Hi-C QC Assessment Main Module
"""

import argparse
import sys
from .config import PairtoolsQCConfig
from .utils import PairtoolsQCLogger
from .calculator import PairtoolsQCCalculator


class PairtoolsQCRunner:
    """Hi-C数据质量控制评估运行器|Hi-C QC Assessment Runner"""

    def __init__(self, **kwargs):
        """初始化运行器|Initialize runner

        Args:
            **kwargs: 配置参数|Configuration parameters
                - pairs_file: Pairs文件路径|Pairs file path
                - pairtools_path: Pairtools可执行文件路径|Pairtools executable path
                - output_dir: 输出目录|Output directory
                - max_unmapped_rate: 未比对阈值|Unmapped threshold
                - max_single_sided_rate: 单端比对阈值|Single-sided threshold
                - min_mapped_rate: 双端比对阈值|Mapped threshold
                - max_dup_rate: PCR重复阈值|Duplication threshold
                - min_cis_trans_ratio: cis/trans阈值|cis/trans threshold
        """
        # 初始化配置|Initialize configuration
        self.config = PairtoolsQCConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        log_file = None
        if hasattr(self.config, 'output_dir'):
            log_file = f"{self.config.output_dir}/pairtools_qc.log"
        self.logger_manager = PairtoolsQCLogger(log_file=log_file)
        self.logger = self.logger_manager.get_logger()

        # 初始化计算器|Initialize calculator
        self.calculator = PairtoolsQCCalculator(self.config, self.logger)

    def run(self):
        """运行分析|Run analysis"""
        self.logger.info("开始Hi-C数据质量评估|Starting Hi-C data quality assessment")

        try:
            # 执行计算|Execute calculation
            results = self.calculator.calculate()

            # 输出结果|Output results
            self.print_report(results)

            # 判断是否通过|Check if passed
            if results['passed']:
                self.logger.info("质量评估通过|Quality assessment passed")
                return 0
            else:
                self.logger.warning("质量评估未通过|Quality assessment not passed")
                return 1

        except Exception as e:
            self.logger.error(f"评估失败|Assessment failed: {str(e)}")
            raise

    def print_report(self, results: dict):
        """打印评估报告|Print assessment report

        Args:
            results: 评估结果字典|Assessment results dictionary
        """
        separator = "-" * 60
        stats = results['stats']
        assessment = results['assessment']

        self.logger.info("")
        self.logger.info(separator)
        self.logger.info("Hi-C数据质量评估报告|Hi-C Data Quality Assessment Report")
        self.logger.info(separator)
        self.logger.info("")

        # 统计数据|Statistics
        self.logger.info("统计数据|Statistics:")
        self.logger.info(f"  总reads数|Total reads: {self.calculator.format_number(stats.get('total_pairs', 0))}")
        self.logger.info(f"  未比对reads|Unmapped reads: {self.calculator.format_number(stats.get('total_unmapped', 0))} ({stats.get('unmapped_rate', 0):.2f}%)")
        self.logger.info(f"  单端比对|Single-sided mapped: {self.calculator.format_number(stats.get('total_single_sided_mapped', 0))} ({stats.get('single_sided_rate', 0):.2f}%)")
        self.logger.info(f"  双端比对|Paired mapped: {self.calculator.format_number(stats.get('total_mapped', 0))} ({stats.get('mapped_rate', 0):.2f}%)")
        self.logger.info(f"  PCR重复|PCR duplicates: {self.calculator.format_number(stats.get('total_dups', 0))} ({stats.get('dup_rate', 0):.2f}%)")
        self.logger.info("")
        self.logger.info(f"  染色体内pairs|Cis pairs: {self.calculator.format_number(stats.get('cis', 0))}")
        self.logger.info(f"  染色体间pairs|Trans pairs: {self.calculator.format_number(stats.get('trans', 0))}")
        self.logger.info(f"  Cis/Trans比例|Cis/Trans ratio: {stats.get('cis_trans_ratio', 0):.2f}")
        self.logger.info("")

        # 质量评估|Quality assessment
        self.logger.info("质量评估|Quality Assessment:")

        metric_names = {
            'unmapped_rate': '未比对reads比例|Unmapped rate',
            'single_sided_rate': '单端比对比例|Single-sided rate',
            'mapped_rate': '双端比对率|Paired mapping rate',
            'dup_rate': 'PCR重复率|Duplication rate',
            'cis_trans_ratio': 'Cis/Trans比例|Cis/Trans ratio'
        }

        for metric_key, metric_name in metric_names.items():
            if metric_key in assessment:
                result = assessment[metric_key]
                status = "通过|PASS" if result['passed'] else "未通过|FAIL"
                self.logger.info(f"  [{status}] {metric_name}: {result['description']}")

        self.logger.info("")

        # 总体结果|Overall result
        if results['passed']:
            self.logger.info("总体评估|Overall: 通过|PASSED")
        else:
            self.logger.warning("总体评估|Overall: 未通过|NOT PASSED - 请检查上述未通过的指标|Please check the failed metrics above")

        self.logger.info("")
        self.logger.info(separator)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="使用pairtools评估Hi-C mapping数据质量|"
                   "Assess Hi-C mapping data quality using pairtools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--input', required=True,
                       help='输入的pairs或BAM文件路径|Input pairs or BAM file path')

    # 可选参数|Optional parameters
    parser.add_argument('-p', '--pairtools-path',
                       default='~/miniforge3/envs/pairtools_v.1.1.3/bin/pairtools',
                       help='Pairtools可执行文件路径|Pairtools executable path')
    parser.add_argument('-o', '--output-dir',
                       default='./pairtools_qc_output',
                       help='输出目录|Output directory')
    parser.add_argument('-c', '--chroms-path',
                       help='Chromosome sizes文件路径（BAM输入时必需）|'
                            'Chromosome sizes file path (required for BAM input)')

    # 质量阈值|Quality thresholds
    parser.add_argument('--max-unmapped-rate',
                       type=float, default=20.0,
                       help='未比对reads比例阈值(%%)|Threshold for unmapped reads rate (%%)')
    parser.add_argument('--max-single-sided-rate',
                       type=float, default=10.0,
                       help='单端比对比例阈值(%%)|Threshold for single-sided mapping rate (%%)')
    parser.add_argument('--min-mapped-rate',
                       type=float, default=80.0,
                       help='双端比对率阈值(%%)|Threshold for paired mapping rate (%%)')
    parser.add_argument('--max-dup-rate',
                       type=float, default=30.0,
                       help='PCR重复率阈值(%%)|Threshold for PCR duplication rate (%%)')
    parser.add_argument('--min-cis-trans-ratio',
                       type=float, default=4.0,
                       help='Cis/Trans比例阈值|Threshold for cis/trans ratio')

    args = parser.parse_args()

    # 创建运行器并运行|Create runner and run
    try:
        runner = PairtoolsQCRunner(
            pairs_file=args.input,
            pairtools_path=args.pairtools_path,
            output_dir=args.output_dir,
            chroms_path=args.chroms_path,
            max_unmapped_rate=args.max_unmapped_rate,
            max_single_sided_rate=args.max_single_sided_rate,
            min_mapped_rate=args.min_mapped_rate,
            max_dup_rate=args.max_dup_rate,
            min_cis_trans_ratio=args.min_cis_trans_ratio
        )
        exit_code = runner.run()
        sys.exit(exit_code)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
