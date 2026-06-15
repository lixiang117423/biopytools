"""
VCF抽样主程序模块|VCF Sampling Main Module
"""

import argparse
import sys
import time
from .config import VCFSamplerConfig
from .utils import VCFSamplerLogger
from .sampler import VCFSamplerCore


class VCFSampler:
    """VCF抽样主类|VCF Sampler Main Class"""

    def __init__(self, **kwargs):
        """
        初始化VCF抽样器|Initialize VCF Sampler

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = VCFSamplerConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = VCFSamplerLogger()
        self.logger = self.logger_manager.get_logger()

        # 初始化抽样核心|Initialize sampling core
        self.sampler_core = VCFSamplerCore(self.config, self.logger)

    def run_sampling(self):
        """
        运行VCF抽样流程|Run VCF sampling pipeline
        """
        start_time = time.time()

        self.logger.info("=" * 60)
        self.logger.info("VCF抽样流程|VCF Sampling Pipeline")
        self.logger.info("=" * 60)
        self.logger.info(f"输入文件|Input file: {self.config.input_vcf}")
        self.logger.info(f"输出文件|Output file: {self.config.output_vcf}")
        self.logger.info(f"抽样比例|Sampling rate: {self.config.sample_rate:.1%}")
        if self.config.random_seed is not None:
            self.logger.info(f"随机种子|Random seed: {self.config.random_seed}")
        self.logger.info("=" * 60)

        try:
            # 步骤1: 统计每条染色体的SNP数量|Step 1: Count SNPs by chromosome
            self.logger.info("步骤1/3|Step 1/3: 统计SNP数量|Count SNPs by chromosome")
            self.logger.info("=" * 60)
            chr_counts = self.sampler_core.count_snp_by_chromosome()
            self.logger.info("步骤1完成|Step 1 completed")

            # 步骤2: 选择要抽取的SNP索引|Step 2: Select SNP indices
            self.logger.info("")
            self.logger.info("步骤2/3|Step 2/3: 选择SNP索引|Select SNP indices")
            self.logger.info("=" * 60)
            selected_indices = self.sampler_core.select_snp_indices(chr_counts)
            self.logger.info("步骤2完成|Step 2 completed")

            # 步骤3: 写入抽样的VCF文件|Step 3: Write sampled VCF
            self.logger.info("")
            self.logger.info("步骤3/3|Step 3/3: 写入抽样结果|Write sampled results")
            self.logger.info("=" * 60)
            total_snp, written_snp = self.sampler_core.write_sampled_vcf(selected_indices)
            self.logger.info("步骤3完成|Step 3 completed")

            # 输出总结信息|Output summary
            elapsed_time = time.time() - start_time
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("抽样总结|Sampling Summary")
            self.logger.info("=" * 60)
            self.logger.info(f"总SNP数|Total SNPs: {total_snp:,}")
            self.logger.info(f"抽取的SNP数|Sampled SNPs: {written_snp:,}")
            self.logger.info(
                f"实际抽样比例|Actual sampling rate: {written_snp/total_snp:.2%} "
                f"(目标|target: {self.config.sample_rate:.1%})"
            )
            self.logger.info(f"运行时间|Runtime: {elapsed_time:.2f} 秒|seconds")
            self.logger.info(f"输出文件|Output file: {self.config.output_vcf}")
            self.logger.info("=" * 60)
            self.logger.info("VCF抽样流程成功完成|VCF sampling pipeline completed successfully")

            return True

        except Exception as e:
            self.logger.error(f"VCF抽样流程失败|VCF sampling pipeline failed: {str(e)}")
            raise


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='VCF抽样工具|VCF Sampling Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i input.vcf.gz -o output.vcf.gz
        """
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|required arguments')
    required.add_argument(
        '-i', '--input', required=True,
        help='输入VCF文件路径|Input VCF file path'
    )
    required.add_argument(
        '-o', '--output', required=True,
        help='输出VCF文件路径|Output VCF file path'
    )

    # 可选参数|Optional arguments
    optional = parser.add_argument_group('可选参数|optional arguments')
    optional.add_argument(
        '-r', '--sample-rate', type=float, default=0.25,
        help='抽样比例|Sampling rate'
    )
    optional.add_argument(
        '-s', '--random-seed', type=int, default=1288,
        help='随机种子|Random seed'
    )

    # 日志参数|Logging parameters
    log_group = parser.add_argument_group('日志参数|logging options')
    log_group.add_argument(
        '-v', '--verbose', action='count', default=0,
        help='详细输出模式|Verbose mode'
    )
    log_group.add_argument(
        '--log-file',
        help='日志文件路径|Log file path'
    )

    # 版本信息|Version information
    parser.add_argument(
        '-V', '--version', action='version',
        version='%(prog)s 1.0.0'
    )

    args = parser.parse_args()

    # 设置日志级别|Set log level
    if args.verbose >= 2:
        log_level = 10  # DEBUG
    elif args.verbose == 1:
        log_level = 20  # INFO
    else:
        log_level = 30  # WARNING

    # 创建抽样器并运行|Create sampler and run
    try:
        sampler = VCFSampler(
            input_vcf=args.input,
            output_vcf=args.output,
            sample_rate=args.sample_rate,
            random_seed=args.random_seed
        )

        # 更新日志级别|Update log level
        if args.verbose >= 2:
            sampler.logger.setLevel(10)
        elif args.verbose == 1:
            sampler.logger.setLevel(20)

        # 如果指定了日志文件，重新初始化日志|Reinitialize logger if log file specified
        if args.log_file:
            sampler.logger_manager = VCFSamplerLogger(log_file=args.log_file, log_level=log_level)
            sampler.logger = sampler.logger_manager.get_logger()

        sampler.run_sampling()

    except KeyboardInterrupt:
        print("[WARNING] 程序被用户中断|Program interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"[ERROR] 程序执行失败|Program execution failed: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
