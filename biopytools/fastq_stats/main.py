"""
FASTQ文件统计主程序模块|FASTQ File Statistics Main Module
"""

import argparse
import sys
import os
from .config import FastqStatsConfig
from .utils import FastqStatsLogger, SeqkitChecker, FastqFileFinder
from .calculator import FastqStatsCalculator


class FastqStatsRunner:
    """FASTQ文件统计运行器|FASTQ File Statistics Runner"""

    def __init__(self, **kwargs):
        """
        初始化运行器|Initialize runner

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = FastqStatsConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = FastqStatsLogger(log_level="INFO")
        self.logger = self.logger_manager.get_logger()

        # 初始化计算器|Initialize calculator
        self.calculator = FastqStatsCalculator(self.config, self.logger)

        # 初始化文件查找器|Initialize file finder
        self.file_finder = FastqFileFinder()

    def check_dependencies(self):
        """检查依赖|Check dependencies"""
        if not SeqkitChecker.check_seqkit():
            self.logger.error("未找到seqkit|seqkit not found")
            self.logger.error("请安装seqkit|Please install seqkit:")
            self.logger.error("  conda install -c bioconda seqkit")
            self.logger.error("  或访问|or visit: https://bioinf.shenwei.me/seqkit/download/")
            sys.exit(1)

        self.logger.info("seqkit依赖检查通过|seqkit dependency check passed")

    def find_samples(self):
        """查找样品文件|Find sample files"""
        self.logger.info(f"正在查找样品文件|Searching for sample files in: {self.config.input_path}")

        samples = self.file_finder.find_paired_files(
            self.config.input_path,
            self.config.pattern
        )

        if not samples:
            self.logger.error("未找到匹配的FASTQ文件|No matching FASTQ files found")
            sys.exit(1)

        self.logger.info(f"找到{len(samples)}个样品|Found {len(samples)} samples")

        for sample, files in samples.items():
            r1_info = f"R1={os.path.basename(files['R1'])}"
            if 'R2' in files:
                r2_info = f", R2={os.path.basename(files['R2'])}"
                self.logger.info(f"  {sample}: {r1_info}{r2_info}")
            else:
                self.logger.info(f"  {sample}: {r1_info} (单端|single-end)")

        return samples

    def run_analysis(self):
        """运行分析|Run analysis"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("FASTQ文件统计工具|FASTQ File Statistics Tool")
            self.logger.info("=" * 60)

            # 检查依赖|Check dependencies
            self.check_dependencies()

            # 查找样品|Find samples
            samples = self.find_samples()

            # 处理样品|Process samples
            self.logger.info(f"使用{self.config.threads}个线程处理|Processing with {self.config.threads} threads")
            results = self.calculator.process_samples(samples)

            if not results:
                self.logger.error("处理失败|Processing failed")
                sys.exit(1)

            # 按样品名称排序|Sort by sample name
            results.sort(key=lambda x: x['sample_name'])

            # 写入结果|Write results
            self.calculator.write_results(results)

            # 打印摘要|Print summary
            self.calculator.print_summary(results)

            self.logger.info("=" * 60)
            self.logger.info(f"完成! 共处理{len(results)}个样品|Completed! Processed {len(results)} samples")
            self.logger.info("=" * 60)

            return True

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='基于seqkit的高性能FASTQ文件统计工具|High-performance FASTQ file statistics tool based on seqkit',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i /data/fastq/ -o results.csv -p "*_1.clean.fq.gz"
  %(prog)s -i sample_R1.fastq.gz -o stats.xlsx -t 8
  %(prog)s -i /data/ -o output.csv -p "*_R1.fastq.gz" -t 16

依赖|Dependencies:
  需要安装seqkit|Requires seqkit: conda install -c bioconda seqkit
  或访问|or visit: https://bioinf.shenwei.me/seqkit/download/
        """
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        help='输入FASTQ文件或包含FASTQ文件的目录|Input FASTQ file or directory containing FASTQ files'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='输出文件路径 (.csv 或 .xlsx)|Output file path (.csv or .xlsx)'
    )
    parser.add_argument(
        '-p', '--pattern',
        default=None,
        help='FASTQ文件匹配模式，如 "*_1.clean.fq.gz"，*代表样品名称|'
             'FASTQ file matching pattern, e.g., "*_1.clean.fq.gz", * represents sample name'
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=12,
        help='线程数|Number of threads (默认|default: 12)'
    )

    args = parser.parse_args()

    # 创建运行器并运行|Create runner and run
    runner = FastqStatsRunner(
        input_path=args.input,
        output_file=args.output,
        pattern=args.pattern,
        threads=args.threads
    )

    runner.run_analysis()


if __name__ == '__main__':
    main()
