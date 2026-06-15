"""
VCF过滤主程序模块|VCF Filtering Main Module
"""

import argparse
import sys
import time
from .config import FilterConfig
from .utils import FilterLogger, CommandRunner, check_dependencies
from .separator import VCFSeparator
from .filter import VCFFilter
from .statistics import VCFStatistics
from .vcf_repair import VCFRepairer

# 版本信息|Version information
VERSION = "1.0.0"


class VCFFilterAnalyzer:
    """VCF过滤分析主类|Main VCF Filter Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = FilterConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = FilterLogger(
            self.config.output_path,
            log_name="vcf_filtering.log",
            log_level=self.config.log_level,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        # 初始化各个处理器|Initialize processors
        self.separator = VCFSeparator(self.config, self.logger, self.cmd_runner)
        self.filter = VCFFilter(self.config, self.logger, self.cmd_runner)
        self.statistics = VCFStatistics(self.config, self.logger)
        self.repairer = VCFRepairer(self.logger)

    def check_dependencies(self):
        """检查依赖软件|Check dependencies"""
        return check_dependencies(self.config, self.logger)

    def run_filtering(self):
        """运行完整的过滤流程|Run complete filtering pipeline"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始VCF SNP/INDEL过滤流程|Starting VCF SNP/INDEL Filtering Pipeline")
            self.logger.info("=" * 60)
            self.logger.info(f"输入VCF文件|Input VCF file: {self.config.vcf_file}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"线程数|Threads: {self.config.threads}")
            self.logger.info(f"BCFtools路径|BCFtools path: {self.config.bcftools_path}")
            self.logger.info(f"变异类型|Variant type: {self.config.variant_type}")
            self.logger.info("=" * 60)

            # 步骤0: 检查并修复VCF文件（如果需要）|Step 0: Check and repair VCF file (if needed)
            if self.config.auto_repair_vcf:
                self.logger.info("检查VCF文件完整性|Checking VCF file integrity")
                repaired_file = f"{self.config.base_name}.repaired.vcf.gz"

                needs_repair, repair_result = self.repairer.check_and_repair(
                    self.config.vcf_file,
                    repaired_file,
                    auto_repair=self.config.auto_repair_vcf
                )

                if needs_repair and repair_result:
                    self.logger.info(f"使用修复后的文件继续分析|Continuing with repaired file: {repair_result}")
                    self.config.vcf_file = repair_result
                elif needs_repair and not repair_result:
                    self.logger.error("VCF文件需要修复但修复失败|VCF file needs repair but repair failed")
                    return False

            # 步骤1: 检查依赖|Step 1: Check dependencies
            self.logger.info("检查依赖软件|Checking dependencies")
            if not self.check_dependencies():
                raise RuntimeError("依赖检查失败|Dependency check failed")

            # 步骤2: 分离SNP和INDEL|Step 2: Separate SNPs and INDELs
            if not self.separator.separate_variants():
                raise RuntimeError("SNP/INDEL分离失败|SNP/INDEL separation failed")

            # 步骤3a: 过滤SNP|Step 3a: Filter SNPs
            if self.config.variant_type in ['both', 'snp_only']:
                if not self.filter.filter_snps(self.separator.raw_snp_file):
                    raise RuntimeError("SNP过滤失败|SNP filtering failed")
            else:
                self.logger.info("跳过SNP过滤|Skipping SNP filtering")

            # 步骤3b: 过滤INDEL|Step 3b: Filter INDELs
            if self.config.variant_type in ['both', 'indel_only']:
                if not self.filter.filter_indels(self.separator.raw_indel_file):
                    raise RuntimeError("INDEL过滤失败|INDEL filtering failed")
            else:
                self.logger.info("跳过INDEL过滤|Skipping INDEL filtering")

            # 步骤3c: 过滤双等位位点SNP|Step 3c: Filter biallelic SNPs
            if self.config.snp_biallelic and self.config.variant_type in ['both', 'snp_only']:
                if not self.filter.filter_biallelic_snps():
                    raise RuntimeError("双等位位点SNP过滤失败|Biallelic SNP filtering failed")
            else:
                self.logger.info("跳过双等位位点SNP过滤|Skipping biallelic SNP filtering")

            # 步骤4: 合并过滤后的变异（仅当variant_type为both时需要合并）|Step 4: Merge filtered variants (only needed when variant_type is both)
            if self.config.variant_type == 'both':
                if not self.filter.merge_filtered_variants_with_biallelic():
                    raise RuntimeError("变异合并失败|Variant merging failed")
            else:
                self.logger.info(f"单一变异类型模式({self.config.variant_type})，跳过合并步骤|Single variant type mode ({self.config.variant_type}), skipping merge step")

            # 步骤5: 生成统计报告|Step 5: Generate statistics report
            if not self.statistics.generate_statistics_report(self.separator, self.filter):
                raise RuntimeError("统计报告生成失败|Statistics report generation failed")

            # 完成|Complete
            self.logger.info("=" * 60)
            self.logger.info("VCF过滤流程成功完成|VCF Filtering Pipeline Completed Successfully!")
            self.logger.info("=" * 60)
            self.logger.info(f"结果保存在|Results saved in: {self.config.output_dir}")
            self.logger.info("=" * 60 + "\n")

        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止|Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)


def main():
    """主函数|Main function"""
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description='VCF SNP/INDEL过滤脚本(模块化版本)|VCF SNP/INDEL Filtering Script (Modular Version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i variants.vcf -o filtered_output
        '''
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True, dest='vcf_file',
                       help='输入VCF文件路径(支持压缩和未压缩)|Input VCF file path (supports compressed and uncompressed)')

    # 可选参数|Optional arguments
    parser.add_argument('-o', '--output-dir', default='./filtered_vcf', dest='output_dir',
                       help='输出目录|Output directory (default: ./filtered_vcf)')
    parser.add_argument('-t', '--threads', type=int, default=64,
                       help='线程数|Number of threads (default: 64)')

    # 变异类型参数|Variant type parameter
    parser.add_argument('--variant-type', dest='variant_type',
                       choices=['both', 'snp_only', 'indel_only'],
                       default='both',
                       help='输入VCF文件的变异类型|Variant type in input VCF (default: both)')

    # SNP过滤参数|SNP filtering parameters
    snp_group = parser.add_argument_group('SNP过滤参数|SNP Filtering Parameters')
    snp_group.add_argument('--snp-qual', type=float, default=30.0,
                          help='SNP最小质量值|SNP minimum QUAL (default: 30.0)')
    snp_group.add_argument('--snp-dp', type=int, default=10,
                          help='SNP最小测序深度|SNP minimum DP (default: 10)')
    snp_group.add_argument('--snp-mq', type=float, default=40.0,
                          help='SNP最小比对质量|SNP minimum MQ (default: 40.0)')
    snp_group.add_argument('--snp-qd', type=float, default=2.0,
                          help='SNP最小质量/深度比|SNP minimum QD (default: 2.0)')
    snp_group.add_argument('--snp-fs', type=float, default=60.0,
                          help='SNP最大FisherStrand值|SNP maximum FS (default: 60.0)')
    snp_group.add_argument('--snp-sor', type=float, default=3.0,
                          help='SNP最大StrandOddsRatio|SNP maximum SOR (default: 3.0)')
    snp_group.add_argument('--snp-mqrs', type=float, default=-12.5,
                          help='SNP最小MappingQualityRankSum|SNP minimum MQRankSum (default: -12.5)')
    snp_group.add_argument('--snp-rprs', type=float, default=-8.0,
                          help='SNP最小ReadPosRankSum|SNP minimum ReadPosRankSum (default: -8.0)')
    snp_group.add_argument('--snp-maf', type=float, default=0.05,
                          help='SNP最小次等位基因频率|SNP minimum MAF (default: 0.05)')
    snp_group.add_argument('--snp-biallelic', action='store_true', default=True,
                          help='只保留双等位位点SNP|Keep only biallelic SNPs (default: True)')
    snp_group.add_argument('--no-snp-biallelic', action='store_false', dest='snp_biallelic',
                          help='不禁用双等位点过滤|Do not filter for biallelic sites')

    # INDEL过滤参数|INDEL filtering parameters
    indel_group = parser.add_argument_group('INDEL过滤参数|INDEL Filtering Parameters')
    indel_group.add_argument('--indel-qual', type=float, default=30.0,
                            help='INDEL最小质量值|INDEL minimum QUAL (default: 30.0)')
    indel_group.add_argument('--indel-dp', type=int, default=10,
                            help='INDEL最小测序深度|INDEL minimum DP (default: 10)')
    indel_group.add_argument('--indel-mq', type=float, default=40.0,
                            help='INDEL最小比对质量|INDEL minimum MQ (default: 40.0)')
    indel_group.add_argument('--indel-qd', type=float, default=2.0,
                            help='INDEL最小质量/深度比|INDEL minimum QD (default: 2.0)')
    indel_group.add_argument('--indel-fs', type=float, default=200.0,
                            help='INDEL最大FisherStrand值|INDEL maximum FS (default: 200.0)')
    indel_group.add_argument('--indel-sor', type=float, default=10.0,
                            help='INDEL最大StrandOddsRatio|INDEL maximum SOR (default: 10.0)')
    indel_group.add_argument('--indel-rprs', type=float, default=-20.0,
                            help='INDEL最小ReadPosRankSum|INDEL minimum ReadPosRankSum (default: -20.0)')

    # 工具路径|Tool paths
    parser.add_argument('--bcftools-path', default='bcftools',
                       help='BCFtools软件路径|BCFtools software path (default: bcftools)')

    # 日志参数|Logging parameters
    log_group = parser.add_argument_group('日志选项|Logging options')
    log_group.add_argument("-v", "--verbose", action="count", default=0,
                          help="详细输出模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)")
    log_group.add_argument("--quiet", action="store_true",
                          help="静默模式(只输出ERROR)|Quiet mode (ERROR only)")
    log_group.add_argument("--log-level",
                          help="日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)|Log level (default: INFO)")
    log_group.add_argument("--log-file",
                          help="日志文件路径|Log file path")

    # 执行控制|Execution control
    exec_group = parser.add_argument_group('执行选项|Execution options')
    exec_group.add_argument("-f", "--force", action="store_true",
                           help="强制覆盖已存在文件|Force overwrite existing files")
    exec_group.add_argument("--dry-run", action="store_true",
                           help="模拟运行(不实际执行)|Dry run without execution")
    exec_group.add_argument("--repair-vcf", action="store_true", dest="auto_repair_vcf",
                           help="自动修复损坏的VCF文件（列数不匹配等问题）|Auto-repair corrupted VCF files (column mismatch, etc.)")

    # 版本信息|Version information
    parser.add_argument("-V", "--version", action="version",
                       version=f'%(prog)s {VERSION}')

    args = parser.parse_args()

    # 确定日志级别|Determine log level
    if args.log_level:
        log_level = args.log_level
    elif args.verbose >= 2:
        log_level = "DEBUG"
    elif args.verbose == 1:
        log_level = "INFO"
    elif args.quiet:
        log_level = "ERROR"
    else:
        log_level = "INFO"

    # 创建分析器并运行|Create analyzer and run
    analyzer = VCFFilterAnalyzer(
        vcf_file=args.vcf_file,
        output_dir=args.output_dir,
        threads=args.threads,
        bcftools_path=args.bcftools_path,
        variant_type=args.variant_type,
        # SNP parameters
        snp_qual=args.snp_qual,
        snp_dp=args.snp_dp,
        snp_mq=args.snp_mq,
        snp_qd=args.snp_qd,
        snp_fs=args.snp_fs,
        snp_sor=args.snp_sor,
        snp_mqrs=args.snp_mqrs,
        snp_rprs=args.snp_rprs,
        snp_maf=args.snp_maf,
        snp_biallelic=args.snp_biallelic,
        # INDEL parameters
        indel_qual=args.indel_qual,
        indel_dp=args.indel_dp,
        indel_mq=args.indel_mq,
        indel_qd=args.indel_qd,
        indel_fs=args.indel_fs,
        indel_sor=args.indel_sor,
        indel_rprs=args.indel_rprs,
        # Logging and execution control
        log_level=log_level,
        quiet=args.quiet,
        verbose=args.verbose,
        force=args.force,
        dry_run=args.dry_run,
        auto_repair_vcf=args.auto_repair_vcf
    )

    try:
        # 输出程序信息|Output program information
        analyzer.logger.info("=" * 60)
        analyzer.logger.info("Program: VCF SNP/INDEL Filtering")
        analyzer.logger.info(f"Version: {VERSION}")
        analyzer.logger.info("=" * 60)

        if args.dry_run:
            analyzer.logger.info("模拟运行模式-不会实际执行命令|DRY RUN mode - commands will not be executed")

        # 执行过滤|Run filtering
        analyzer.run_filtering()

        # 输出总结信息|Output summary
        elapsed_time = time.time() - start_time
        analyzer.logger.info("=" * 60)
        analyzer.logger.info("Pipeline Summary")
        analyzer.logger.info("=" * 60)
        analyzer.logger.info(f"Total runtime: {elapsed_time:.2f} seconds")
        analyzer.logger.info(f"Output directory: {args.output_dir}")
        analyzer.logger.info("Pipeline completed successfully")

    except KeyboardInterrupt:
        analyzer.logger.warning("用户中断程序执行|Process interrupted by user")
        sys.exit(130)
    except Exception as e:
        analyzer.logger.critical(f"Pipeline failed: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
