"""
覆盖度过滤主程序|Coverage Filter Main Program
"""

import argparse
import sys
import os
from datetime import datetime
from pathlib import Path

# 确保可以导入模块|Ensure module can be imported
sys.path.insert(0, str(Path(__file__).parent.parent))

from .config import CoverageFilterConfig
from .utils import CoverageFilterLogger, check_dependencies
from .calculator import CoverageCalculator


class CoverageFilter:
    """覆盖度过滤器|Coverage Filter"""

    def __init__(self, **kwargs):
        """初始化过滤器|Initialize filter"""
        self.config = CoverageFilterConfig(**kwargs)
        self.config.validate()

        # 设置日志|Setup logging
        log_dir = Path(self.config.output_prefix).parent
        self.logger_manager = CoverageFilterLogger(log_dir, "coverage_filter.log")
        self.logger = self.logger_manager.get_logger()

        # 初始化计算器|Initialize calculator
        self.calculator = CoverageCalculator(self.config, self.logger)

    def run_filter(self):
        """运行过滤流程|Run filter pipeline"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("覆盖度过滤流程|Coverage Filter Pipeline")
            self.logger.info("=" * 60)

            # 检查依赖|Check dependencies
            if not check_dependencies(self.logger):
                return False

            # 步骤1: 计算覆盖度|Step 1: Calculate coverage
            self.logger.info("[1/4] 计算覆盖度|Calculating coverage")
            coverage_file = self.calculator.calculate_coverage()
            if not coverage_file:
                return False

            # 步骤2: 分类序列|Step 2: Classify sequences
            self.logger.info("[2/4] 分类序列质量|Classifying sequence quality")
            high_list, medium_list, low_list = self.calculator.classify_sequences(coverage_file)
            if not high_list:
                return False

            # 获取统计数量|Get count statistics
            high_count = self.calculator._count_lines(high_list)
            medium_count = self.calculator._count_lines(medium_list)
            low_count = self.calculator._count_lines(low_list)

            # 步骤3: 提取序列|Step 3: Extract sequences
            self.logger.info("[3/4] 提取过滤后的序列|Extracting filtered sequences")
            if not self.calculator.extract_sequences(high_list, medium_list, low_list):
                return False

            # 步骤4: 生成报告|Step 4: Generate report
            self.logger.info("[4/4] 生成统计报告|Generating statistics report")
            stats = self.calculator.generate_statistics()

            # 生成详细报告文件|Generate detailed report file
            report_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_report.txt")
            self._generate_report(high_count, medium_count, low_count, stats, report_file)

            self.logger.info("=" * 60)
            self.logger.info("过滤流程完成|Filter pipeline completed")
            self.logger.info("=" * 60)

            return True

        except KeyboardInterrupt:
            self.logger.info("用户中断操作|User interrupted operation")
            return False
        except Exception as e:
            self.logger.error(f"过滤过程中发生错误|Error occurred during filtering: {str(e)}")
            return False

    def _generate_report(self, high_count, medium_count, low_count, stats, report_file):
        """生成报告文件|Generate report file"""
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 60 + "\n")
            f.write("覆盖度过滤报告|Coverage Filter Report\n")
            f.write("=" * 60 + "\n")
            f.write(f"生成时间|Generated time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            f.write("-" * 60 + "\n")
            f.write("配置参数|Configuration Parameters\n")
            f.write("-" * 60 + "\n")
            f.write(f"BAM文件|BAM file: {self.config.bam_file}\n")
            f.write(f"FASTA文件|FASTA file: {self.config.fasta_file}\n")
            f.write(f"输出目录|Output directory: {self.config.output_dir}\n")
            f.write(f"高质量阈值|High quality threshold: "
                   f"覆盖度≥{self.config.high_coverage}%\n")
            f.write(f"中等质量阈值|Medium quality thresholds: "
                   f"覆盖度{self.config.medium_cov_min}-{self.config.high_coverage}%\n")
            f.write(f"低质量阈值|Low quality threshold: "
                   f"覆盖度<{self.config.medium_cov_min}%\n\n")

            f.write("-" * 60 + "\n")
            f.write("过滤结果|Filter Results\n")
            f.write("-" * 60 + "\n")
            f.write(f"高质量序列|High quality sequences: {high_count}\n")
            f.write(f"中等质量序列|Medium quality sequences: {medium_count}\n")
            f.write(f"低质量序列|Low quality sequences: {low_count}\n\n")

            f.write("-" * 60 + "\n")
            f.write("统计信息|Statistics\n")
            f.write("-" * 60 + "\n\n")

            f.write("原始序列|Original sequences:\n")
            f.write(stats.get('original', '统计失败|Statistics failed') + "\n\n")

            f.write("高质量序列|High quality sequences:\n")
            f.write(stats.get('high_quality', '统计失败|Statistics failed') + "\n\n")

            f.write("中等质量序列|Medium quality sequences:\n")
            f.write(stats.get('medium_quality', '统计失败|Statistics failed') + "\n\n")

            f.write("低质量序列|Low quality sequences:\n")
            f.write(stats.get('low_quality', '统计失败|Statistics failed') + "\n\n")

            f.write("-" * 60 + "\n")
            f.write("输出文件|Output Files\n")
            f.write("-" * 60 + "\n")
            f.write(f"- {os.path.join(self.config.output_dir, self.config.output_prefix + '_high_quality.fa')}    : 推荐用于后续分析|Recommended for downstream analysis\n")
            f.write(f"- {os.path.join(self.config.output_dir, self.config.output_prefix + '_medium_quality.fa')}  : 中等质量序列|Medium quality sequences\n")
            f.write(f"- {os.path.join(self.config.output_dir, self.config.output_prefix + '_low_quality.fa')}     : 可疑序列，建议进一步检查|Suspicious sequences, further inspection recommended\n")
            f.write(f"- {os.path.join(self.config.output_dir, self.config.output_prefix + '_coverage.txt')}       : 详细覆盖度数据|Detailed coverage data\n")
            f.write(f"- {os.path.join(self.config.output_dir, self.config.output_prefix + '_*_quality.list')}     : 各类序列ID列表|Sequence ID lists for each category\n")
            f.write(f"- {report_file}                                  : 本报告|This report\n")

        self.logger.info(f"报告已保存到|Report saved to: {report_file}")


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="覆盖度过滤工具|Coverage Filter Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i sample.bam -f genome.fa -o filtered
        '''
    )

    parser.add_argument('-i', '--bam-file',
                        required=True,
                        help='BAM文件路径|BAM file path')

    parser.add_argument('-f', '--fasta-file',
                        required=True,
                        help='基因组FASTA文件|Genome FASTA file')

    parser.add_argument('-o', '--output-prefix',
                        required=True,
                        help='输出文件前缀|Output file prefix')

    parser.add_argument('-d', '--output-dir',
                        default='.',
                        help='输出目录|Output directory (default: current directory)')

    parser.add_argument('-t', '--threads',
                        type=int,
                        default=12,
                        help='线程数|Number of threads (default: 12)')

    parser.add_argument('--high-cov',
                        type=float,
                        default=90.0,
                        help='高质量覆盖度阈值|High quality coverage threshold (default: 90.0)')

    parser.add_argument('--medium-cov-min',
                        type=float,
                        default=50.0,
                        help='中等质量最小覆盖度|Medium quality minimum coverage (default: 50.0)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    try:
        # 解析参数|Parse arguments
        args = parse_arguments()

        # 配置过滤器|Configure filter
        filter_config = {
            'bam_file': args.bam_file,
            'fasta_file': args.fasta_file,
            'output_prefix': args.output_prefix,
            'output_dir': args.output_dir,
            'threads': args.threads,
            'high_coverage': args.high_cov,
            'medium_cov_min': args.medium_cov_min,
        }

        # 创建并运行过滤器|Create and run filter
        filter_tool = CoverageFilter(**filter_config)
        success = filter_tool.run_filter()

        if success:
            report_path = os.path.join(args.output_dir, f"{args.output_prefix}_report.txt")
            print(f"\n过滤完成|Filter completed! 报告已保存至|Report saved to: {report_path}")
            sys.exit(0)
        else:
            print("\n过滤失败|Filter failed")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n用户中断操作|User interrupted operation")
        sys.exit(1)
    except Exception as e:
        print(f"程序错误|Program error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
