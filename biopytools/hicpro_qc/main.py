"""
HiC-Pro质量控制评估主程序模块|HiC-Pro QC Assessment Main Module
"""

import argparse
import sys
from .config import HiCProQCConfig
from .utils import HiCProQCLogger
from .calculator import HiCProQCCalculator


class HiCProQCRunner:
    """HiC-Pro质量控制评估运行器|HiC-Pro QC Assessment Runner"""

    def __init__(self, **kwargs):
        """初始化运行器|Initialize runner

        Args:
            **kwargs: 配置参数|Configuration parameters
                - hicpro_dir: HiC-Pro输出目录|HiC-Pro output directory
                - output_dir: 输出目录|Output directory
                - sample_name: 样本名称|Sample name (optional)
                - min_mapping_rate: 最低比对率阈值|Minimum mapping rate threshold
                - min_unique_rate: 最低唯一比对率阈值|Minimum unique mapping rate threshold
                - min_valid_pairs_rate: 最低valid pairs比例阈值|Minimum valid pairs rate threshold
                - max_dangling_ends_rate: 最高dangling ends比例阈值|Maximum dangling ends rate threshold
                - max_self_ligation_rate: 最高self-ligation比例阈值|Maximum self-ligation rate threshold
                - max_religation_rate: 最高religation比例阈值|Maximum religation rate threshold
                - min_cis_trans_ratio: 最低cis/trans比例阈值|Minimum cis/trans ratio threshold
                - max_duplication_rate: 最高PCR重复率阈值|Maximum PCR duplication rate threshold
        """
        # 初始化配置|Initialize configuration
        self.config = HiCProQCConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        log_file = None
        if hasattr(self.config, 'output_dir'):
            log_file = f"{self.config.output_dir}/hicpro_qc.log"
        self.logger_manager = HiCProQCLogger(log_file=log_file)
        self.logger = self.logger_manager.get_logger()

        # 初始化计算器|Initialize calculator
        self.calculator = HiCProQCCalculator(self.config, self.logger)

    def run(self):
        """运行分析|Run analysis"""
        self.logger.info("开始HiC-Pro数据质量评估|Starting HiC-Pro data quality assessment")

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
        metrics = results['metrics']
        assessment = results['assessment']

        self.logger.info("")
        self.logger.info(separator)
        self.logger.info("HiC-Pro数据质量评估报告|HiC-Pro Data Quality Assessment Report")
        self.logger.info(separator)
        self.logger.info("")

        # 统计数据|Statistics
        self.logger.info("统计数据|Statistics:")

        # Mapping统计|Mapping statistics
        if 'total_reads' in metrics:
            self.logger.info(f"  总reads数|Total reads: {self.calculator.format_number(metrics.get('total_reads', 0))}")
            self.logger.info(f"  比对reads|Mapped reads: {self.calculator.format_number(metrics.get('total_mapped', 0))} ({metrics.get('mapping_rate', 0):.2f}%)")
            self.logger.info(f"  唯一比对|Uniquely mapped: {self.calculator.format_number(metrics.get('total_global_unique', 0))} ({metrics.get('unique_mapping_rate', 0):.2f}%)")
            self.logger.info("")

        # Valid pairs统计|Valid pairs statistics
        if 'total_pairs' in metrics:
            self.logger.info(f"  总pairs数|Total pairs: {self.calculator.format_number(metrics.get('total_pairs', 0))}")
            self.logger.info(f"  Valid pairs: {self.calculator.format_number(metrics.get('valid_pairs', 0))} ({metrics.get('valid_pairs_rate', 0):.2f}%)")
            self.logger.info(f"  Dangling ends: {self.calculator.format_number(metrics.get('dangling_ends', 0))} ({metrics.get('dangling_ends_rate', 0):.2f}%)")
            self.logger.info(f"  Self-ligation: {self.calculator.format_number(metrics.get('self_ligation', 0))} ({metrics.get('self_ligation_rate', 0):.2f}%)")
            self.logger.info(f"  Religation: {self.calculator.format_number(metrics.get('religation', 0))} ({metrics.get('religation_rate', 0):.2f}%)")
            self.logger.info("")

        # Cis/Trans|Cis/Trans ratio
        if 'cis_trans_ratio' in metrics:
            self.logger.info(f"  Cis/Trans比例|Cis/Trans ratio: {metrics.get('cis_trans_ratio', 0):.2f}")
            self.logger.info("")

        # 质量评估|Quality assessment
        self.logger.info("质量评估|Quality Assessment:")

        metric_names = {
            'mapping_rate': '比对率|Mapping rate',
            'unique_mapping_rate': '唯一比对率|Unique mapping rate',
            'valid_pairs_rate': 'Valid pairs比例|Valid pairs rate',
            'dangling_ends_rate': 'Dangling ends比例|Dangling ends rate',
            'self_ligation_rate': 'Self-ligation比例|Self-ligation rate',
            'religation_rate': 'Religation比例|Religation rate',
            'duplication_rate': 'PCR重复率|Duplication rate',
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
        description="使用HiC-Pro统计文件评估Hi-C数据质量|"
                   "Assess Hi-C data quality using HiC-Pro statistics files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--input', required=True,
                       help='HiC-Pro输出目录路径|HiC-Pro output directory path')

    # 可选参数|Optional parameters
    parser.add_argument('-o', '--output-dir',
                       default='./hicpro_qc_output',
                       help='输出目录|Output directory')
    parser.add_argument('-s', '--sample-name',
                       help='样本名称（可选，默认自动检测）|Sample name (optional, auto-detect by default)')

    # 质量阈值|Quality thresholds
    parser.add_argument('--min-mapping-rate',
                       type=float, default=70.0,
                       help='最低比对率阈值(%%)|Minimum mapping rate threshold (%%)')
    parser.add_argument('--min-unique-rate',
                       type=float, default=60.0,
                       help='最低唯一比对率阈值(%%)|Minimum unique mapping rate threshold (%%)')
    parser.add_argument('--min-valid-pairs-rate',
                       type=float, default=50.0,
                       help='最低valid pairs比例阈值(%%)|Minimum valid pairs rate threshold (%%)')
    parser.add_argument('--max-dangling-ends-rate',
                       type=float, default=15.0,
                       help='最高dangling ends比例阈值(%%)|Maximum dangling ends rate threshold (%%)')
    parser.add_argument('--max-self-ligation-rate',
                       type=float, default=5.0,
                       help='最高self-ligation比例阈值(%%)|Maximum self-ligation rate threshold (%%)')
    parser.add_argument('--max-religation-rate',
                       type=float, default=10.0,
                       help='最高religation比例阈值(%%)|Maximum religation rate threshold (%%)')
    parser.add_argument('--min-cis-trans-ratio',
                       type=float, default=5.0,
                       help='最低cis/trans比例阈值|Minimum cis/trans ratio threshold')
    parser.add_argument('--max-duplication-rate',
                       type=float, default=30.0,
                       help='最高PCR重复率阈值(%%)|Maximum PCR duplication rate threshold (%%)')

    args = parser.parse_args()

    # 创建运行器并运行|Create runner and run
    try:
        runner = HiCProQCRunner(
            hicpro_dir=args.input,
            output_dir=args.output_dir,
            sample_name=args.sample_name,
            min_mapping_rate=args.min_mapping_rate,
            min_unique_rate=args.min_unique_rate,
            min_valid_pairs_rate=args.min_valid_pairs_rate,
            max_dangling_ends_rate=args.max_dangling_ends_rate,
            max_self_ligation_rate=args.max_self_ligation_rate,
            max_religation_rate=args.max_religation_rate,
            min_cis_trans_ratio=args.min_cis_trans_ratio,
            max_duplication_rate=args.max_duplication_rate
        )
        exit_code = runner.run()
        sys.exit(exit_code)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
