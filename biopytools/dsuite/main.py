"""
Dsuite Main Module
Dsuite主模块
"""

import os
import sys
import time
import argparse
from .config import DsuiteConfig
from .utils import DsuiteLogger, VCFStatsCollector, DsuiteRunner


def setup_logger(output_dir: str):
    """设置日志 | Setup logging"""
    logger_manager = DsuiteLogger(output_dir)
    return logger_manager.get_logger()


def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='Dsuite D-statistic Analysis - Dsuite D统计量分析',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例 | Examples:

  # 基本用法
  %(prog)s -i variants.vcf.gz -s sets.txt -o output_dir

  # 指定输出前缀
  %(prog)s -i variants.vcf.gz -s sets.txt -o output_dir -p analysis

  # 自定义过滤参数
  %(prog)s -i variants.vcf.gz -s sets.txt -o output_dir \\
      --min-alleles 2 --max-alleles 2 --variant-type snps
        '''
    )

    # 必需参数 | Required parameters
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--vcf', required=True,
                         help='[FILE] 输入VCF文件路径 | Input VCF file path')
    required.add_argument('-s', '--sets', required=True,
                         help='[FILE] SETS分组文件路径 | SETS file path')
    required.add_argument('-o', '--output', required=True,
                         help='[DIR] 输出目录 | Output directory')

    # 可选参数 | Optional parameters
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-p', '--prefix', default='dsuite',
                         help='[STR] 输出文件前缀 (默认: dsuite) | Output file prefix (default: dsuite)')
    optional.add_argument('--dsuite-bin',
                         default='/share/org/YZWL/yzwl_lixg/software/Dsuite/Build/Dsuite',
                         help='[FILE] Dsuite可执行文件路径 | Dsuite binary path')
    optional.add_argument('--min-alleles', type=int, default=2,
                         help='[INT] 最小等位基因数 (默认: 2) | Min number of alleles (default: 2)')
    optional.add_argument('--max-alleles', type=int, default=2,
                         help='[INT] 最大等位基因数 (默认: 2) | Max number of alleles (default: 2)')
    optional.add_argument('--variant-type', default='snps',
                         choices=['snps', 'indels', 'both', 'none'],
                         help='[STR] 变异类型 (默认: snps) | Variant type (default: snps)')
    optional.add_argument('--bcftools', default='bcftools',
                         help='[CMD] bcftools命令路径 | bcftools command path')
    optional.add_argument('--collect-stats', action='store_true',
                         help='[FLAG] 是否收集VCF统计信息 | Whether to collect VCF statistics')

    args = parser.parse_args()

    # 创建输出目录 | Create output directory
    os.makedirs(args.output, exist_ok=True)

    # 设置日志 | Setup logging
    logger = setup_logger(args.output)

    start_time = time.time()

    try:
        logger.info("=" * 60)
        logger.info("Dsuite D-trios Analysis Pipeline")
        logger.info("=" * 60)

        # 初始化配置 | Initialize config
        config = DsuiteConfig(
            vcf_file=args.vcf,
            sets_file=args.sets,
            output_dir=args.output,
            output_prefix=args.prefix,
            dsuite_bin=args.dsuite_bin,
            min_alleles=args.min_alleles,
            max_alleles=args.max_alleles,
            variant_type=args.variant_type,
            bcftools=args.bcftools,
            collect_stats=args.collect_stats
        )

        config.validate()

        # 初始化工具类 | Initialize utilities
        stats_collector = VCFStatsCollector(logger)
        dsuite_runner = DsuiteRunner(logger)

        # 收集VCF统计信息（可选）| Collect VCF statistics (optional)
        stats = {}
        if config.collect_stats:
            stats = stats_collector.collect_statistics(config.vcf_file, config.bcftools)

            # 统计过滤后的变异数 | Count filtered variants
            numlines = stats_collector.count_filtered_variants(
                config.vcf_file,
                config.bcftools,
                config.min_alleles,
                config.max_alleles,
                config.variant_type,
                verbose=True  # 收集统计时显示详细信息
            )

            if numlines == 0:
                logger.error("过滤后没有变异位点！请检查过滤参数")
                return 1

            # 计算过滤比例 | Calculate filter ratio
            if 'total_variants' in stats:
                filter_ratio = (numlines / stats['total_variants']) * 100
                logger.info(f"过滤保留比例: {filter_ratio:.2f}%")
        else:
            # 不收集统计信息时，静默统计过滤后的变异数 | Silently count filtered variants when not collecting stats
            numlines = stats_collector.count_filtered_variants(
                config.vcf_file,
                config.bcftools,
                config.min_alleles,
                config.max_alleles,
                config.variant_type,
                verbose=False  # 不显示详细日志
            )

            if numlines == 0:
                logger.error("过滤后没有变异位点！请检查过滤参数")
                return 1

        # 构建完整输出前缀 | Build full output prefix
        output_prefix = os.path.join(config.output_dir, config.output_prefix)

        # 运行Dsuite Dtrios | Run Dsuite Dtrios
        if not dsuite_runner.run_dtrios(
            config.vcf_file,
            config.sets_file,
            output_prefix,
            numlines,
            config.dsuite_bin,
            config.bcftools,
            config.min_alleles,
            config.max_alleles,
            config.variant_type
        ):
            logger.error("Dsuite运行失败")
            return 1

        # 检查输出文件 | Check output files
        dsuite_runner.check_output_files(output_prefix)

        # 完成统计 | Completion statistics
        elapsed_time = time.time() - start_time

        logger.info("")
        logger.info("=" * 60)
        logger.info("分析总结")
        logger.info("=" * 60)
        logger.info(f"输入VCF: {config.vcf_file}")
        logger.info(f"SETS文件: {config.sets_file}")
        logger.info(f"输出目录: {config.output_dir}")

        # 只在收集统计信息时显示详细统计 | Show detailed stats only when collected
        if config.collect_stats and stats:
            logger.info(f"样本数量: {stats.get('sample_count', 'N/A')}")
            logger.info(f"总变异数: {stats.get('total_variants', 'N/A')}")

        logger.info(f"过滤后变异数: {numlines}")
        logger.info(f"运行时间: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
        logger.info("")
        logger.info("下一步建议:")
        logger.info(f"  1. 查看D统计结果: cat {output_prefix}_BBAA.txt")
        logger.info(f"  2. 查看树文件: cat {output_prefix}_tree.txt")
        logger.info("  3. 使用 Dsuite Dinvestigate 进行进一步分析")
        logger.info("=" * 60)
        logger.info("分析完成！")
        logger.info("=" * 60)

        return 0

    except KeyboardInterrupt:
        logger.warning("用户中断操作")
        return 130
    except ValueError as e:
        logger.error(f"配置错误: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"Pipeline失败: {str(e)}", exc_info=True)
        return 1


class DsuiteAnalyzer:
    """Dsuite分析器类 | Dsuite Analyzer Class"""

    def __init__(self, **kwargs):
        """
        初始化分析器 | Initialize analyzer

        Args:
            **kwargs: 配置参数 | Configuration parameters
        """
        self.config = DsuiteConfig(**kwargs)
        self.config.validate()

        # 初始化日志 | Initialize logging
        self.logger_manager = DsuiteLogger(self.config.output_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化工具类 | Initialize utilities
        self.stats_collector = VCFStatsCollector(self.logger)
        self.dsuite_runner = DsuiteRunner(self.logger)

    def run(self):
        """
        运行分析 | Run analysis

        Returns:
            是否成功 | Whether successful
        """
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("Dsuite D-trios Analysis Pipeline")
            self.logger.info("=" * 60)

            # 收集VCF统计信息（可选）| Collect VCF statistics (optional)
            stats = {}
            if self.config.collect_stats:
                stats = self.stats_collector.collect_statistics(
                    self.config.vcf_file,
                    self.config.bcftools
                )

            # 统计过滤后的变异数 | Count filtered variants
            numlines = self.stats_collector.count_filtered_variants(
                self.config.vcf_file,
                self.config.bcftools,
                self.config.min_alleles,
                self.config.max_alleles,
                self.config.variant_type,
                verbose=self.config.collect_stats  # 只在收集统计时显示详细信息
            )

            if numlines == 0:
                self.logger.error("过滤后没有变异位点！")
                return False

            # 构建完整输出前缀 | Build full output prefix
            output_prefix = os.path.join(
                self.config.output_dir,
                self.config.output_prefix
            )

            # 运行Dsuite Dtrios | Run Dsuite Dtrios
            if not self.dsuite_runner.run_dtrios(
                self.config.vcf_file,
                self.config.sets_file,
                output_prefix,
                numlines,
                self.config.dsuite_bin,
                self.config.bcftools,
                self.config.min_alleles,
                self.config.max_alleles,
                self.config.variant_type
            ):
                return False

            # 检查输出文件 | Check output files
            self.dsuite_runner.check_output_files(output_prefix)

            # 完成统计 | Completion statistics
            elapsed_time = time.time() - start_time

            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("分析总结")
            self.logger.info("=" * 60)
            self.logger.info(f"运行时间: {elapsed_time:.2f} seconds")
            self.logger.info("分析完成！")

            return True

        except Exception as e:
            self.logger.error(f"分析失败: {str(e)}", exc_info=True)
            return False


if __name__ == "__main__":
    sys.exit(main())
