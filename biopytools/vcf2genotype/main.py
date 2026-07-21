"""
VCF基因型提取主程序模块|VCF Genotype Extraction Main Module
"""

import argparse
import sys
from collections import defaultdict
from .config import VCFConfig
from .utils import VCFLogger, check_dependencies
from .vcf_parser import VCFParserFactory
from .genotype_processor import GenotypeProcessor
from .output_formatter import OutputFormatter
from pathlib import Path

# 标准列名顺序|Standard column name order
STANDARD_FIELDS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL']


class VCFGenotypeExtractor:
    """VCF基因型提取主类|Main VCF Genotype Extractor Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = VCFConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = VCFLogger(
            Path(self.config.output_dir),
            verbose=self.config.verbose,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 检查依赖|Check dependencies
        self.dependencies = check_dependencies()
        self._log_dependencies()

        # 初始化各个处理器|Initialize processors
        self.parser = VCFParserFactory.create_parser(
            self.config.vcf_file,
            self.logger,
            use_fast=self.dependencies['cyvcf2']
        )
        self.processor = GenotypeProcessor(self.config, self.logger)
        self.formatter = OutputFormatter(self.config, self.logger)

    def _log_dependencies(self):
        """记录依赖信息|Log dependencies info"""
        self.logger.info("依赖检查|Dependency check:")
        for dep, available in self.dependencies.items():
            status = "可用|Available" if available else "不可用|Not available"
            self.logger.info(f"  {dep}: {status}")

    def _get_target_samples(self):
        """获取目标样本列表|Get target sample list"""
        all_samples = self.parser.get_samples()

        if self.config.samples == "all":
            return all_samples
        elif isinstance(self.config.samples, list):
            # 验证样本是否存在|Validate sample existence
            valid_samples = []
            for sample in self.config.samples:
                if sample in all_samples:
                    valid_samples.append(sample)
                else:
                    self.logger.warning(f"样本 '{sample}' 在VCF文件中未找到|Sample '{sample}' not found in VCF file")
            return valid_samples
        else:
            return all_samples

    def _build_fieldnames(self, target_samples, gt_columns):
        """构建输出列名列表|Build output fieldnames list"""
        fields = list(STANDARD_FIELDS)
        # QUAL后紧跟纯合/杂合比例|Homozygous/Heterozygous ratio right after QUAL
        fields.extend(['Homozygous_Ratio', 'Heterozygous_Ratio'])
        # 基因型统计列|Genotype count columns
        fields.extend(gt_columns)
        # 样本列|Sample columns
        fields.extend(target_samples)
        return fields

    def _collect_gt_types(self, target_samples):
        """全扫描收集基因型类型,并统计通过过滤的变异数/染色体分布|Full scan: collect GT types + count kept variants per chrom"""
        gt_types = set()
        total_variants = 0
        chrom_counts = defaultdict(int)
        for row in self.parser.parse_records(target_samples):
            if not self.processor.should_keep(row):
                continue
            total_variants += 1
            chrom_counts[row['CHROM']] += 1
            for key in row:
                if key.startswith('GT_'):
                    gt_types.add(key)
        return sorted(gt_types), total_variants, chrom_counts

    def run_extraction(self):
        """运行基因型提取（两遍扫描：第一遍收集GT类型，第二遍流式写出）|Run genotype extraction (two-pass: collect GT types, then stream)"""
        try:
            self.logger.info("开始VCF基因型提取|Starting VCF genotype extraction")
            self.logger.info(f"输入文件|Input file: {self.config.vcf_file}")
            self.logger.info(f"输出前缀|Output prefix: {self.config.output_prefix}")
            self.logger.info(f"输出格式|Output format: {self.config.output_type}")

            # 获取目标样本|Get target samples
            target_samples = self._get_target_samples()
            if not target_samples:
                self.logger.error("没有有效的目标样本|No valid target samples")
                return False

            self.logger.info(f"目标样本数|Number of target samples: {len(target_samples)}")

            # 第一遍：全扫描收集基因型类型 + 统计|First pass: collect GT types + counts
            self.logger.info("扫描基因型类型|Scanning genotype types")
            gt_columns, total_variants, chrom_counts = self._collect_gt_types(target_samples)
            self.logger.info(f"基因型类型|Genotype types: {gt_columns}")
            self.logger.info(f"通过过滤的变异数|Variants passing filter: {total_variants}")

            # 构建完整的列名列表|Build complete fieldnames list
            fieldnames = self._build_fieldnames(target_samples, gt_columns)

            # 判断是否按染色体拆分(config已归一为bool)|Check split (config normalizes to bool)
            should_split = bool(self.config.split_by_chromosome)

            # 试运行模式:只扫描统计,不写数据文件|Dry run: scan + stats only, no data files
            if self.config.dry_run:
                self.logger.info("试运行(DRY RUN):仅扫描统计,不写出数据文件|Dry run: scan + stats only, no data files written")
                stats = {
                    'total_variants': total_variants,
                    'chromosomes': sorted(chrom_counts.keys()),
                    'chromosome_counts': dict(chrom_counts),
                }
                self.formatter.write_summary(stats)
                self.logger.info("试运行完成(未生成数据文件)|Dry run completed (no data files generated)")
                return True

            # 汇总统计：仅用计数器|Summary stats: counters only
            total_variants = 0
            chrom_counts = defaultdict(int)

            # 打开输出文件|Open output files
            if should_split:
                chrom_streams = {}  # chrom -> [formatter_stream, row_count]
            else:
                self.formatter.open_stream(suffix="", fieldnames=fieldnames)

            # 第二遍：重新创建parser（cyvcf2迭代器不可重用）|Second pass: recreate parser (cyvcf2 iterator not reusable)
            self.parser = VCFParserFactory.create_parser(
                self.config.vcf_file,
                self.logger,
                use_fast=self.dependencies['cyvcf2']
            )

            # 第二遍：流式写出|Second pass: streaming output
            for row in self.parser.parse_records(target_samples):
                if not self.processor.should_keep(row):
                    continue

                chrom = row['CHROM']
                total_variants += 1
                chrom_counts[chrom] += 1

                # 补全未出现的GT_列为0|Fill missing GT_ columns with 0
                for gt_col in gt_columns:
                    row.setdefault(gt_col, 0)

                if should_split:
                    if chrom not in chrom_streams:
                        stream = OutputFormatter(self.config, self.logger)
                        stream.open_stream(suffix=chrom, fieldnames=fieldnames)
                        chrom_streams[chrom] = [stream, 0]
                    chrom_streams[chrom][0].write_row(row)
                    chrom_streams[chrom][1] += 1
                else:
                    self.formatter.write_row(row)

                if total_variants % 10000 == 0:
                    self.logger.info(f"已处理|Processed {total_variants} variants")

            # 关闭输出文件|Close output files
            if should_split:
                for chrom, (stream, count) in chrom_streams.items():
                    stream.close_stream(count)
            else:
                self.formatter.close_stream(total_variants)

            self.logger.info(f"处理完成，共处理|Processing completed, total processed: {total_variants} variants")

            # 生成汇总|Generate summary
            stats = {
                'total_variants': total_variants,
                'chromosomes': sorted(chrom_counts.keys()),
                'chromosome_counts': dict(chrom_counts),
            }
            self.formatter.write_summary(stats)

            self.logger.info("VCF基因型提取完成|VCF genotype extraction completed")
            return True

        except Exception as e:
            self.logger.error(f"提取过程中发生错误|Error during extraction: {e}")
            return False


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='VCF基因型提取工具|VCF Genotype Extraction Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i variants.vcf -o genotypes
        '''
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-i', '--input',
                         dest='vcf_file',
                         required=True,
                         help='VCF文件路径(支持.gz压缩格式)|VCF file path (supports .gz compressed format)')

    # 输出配置|Output configuration
    output = parser.add_argument_group('输出配置|Output configuration')
    output.add_argument('-o', '--output', default='vcf_genotype',
                       help='输出文件前缀|Output file prefix')
    output.add_argument('-t', '--output-type', choices=['txt', 'csv'], default='txt',
                       help='输出文件格式|Output file format')
    output.add_argument('--output-dir', default='./',
                       help='输出目录|Output directory')

    # 样本选择|Sample selection
    sample = parser.add_argument_group('样本选择|Sample selection')
    sample.add_argument('-s', '--samples', default='all',
                       help='样本选择：all（所有样本）或逗号分隔的样本名称|Sample selection: all (all samples) or comma-separated sample names')

    # 过滤选项|Filtering options
    filtering = parser.add_argument_group('过滤选项|Filtering options')
    filtering.add_argument('--biallelic-only', action='store_true',
                          help='只保留双等位位点|Keep only biallelic sites')
    filtering.add_argument('--min-length', type=int, default=None,
                          help='最小变异长度（包含），默认不限制|Minimum variant length (inclusive), default no limit')
    filtering.add_argument('--max-length', type=int, default=None,
                          help='最大变异长度（包含），默认不限制|Maximum variant length (inclusive), default no limit')

    # 高级选项|Advanced options
    advanced = parser.add_argument_group('高级选项|Advanced options')
    advanced.add_argument('-e', '--each', choices=['yes', 'y', 'no', 'n'], default='n',
                         help='按染色体拆分输出文件：yes/y（是）或no/n（否）|Split output files by chromosome: yes/y or no/n')
    advanced.add_argument('--dry-run', action='store_true',
                         help='试运行模式，不实际执行操作|Dry run mode, no actual operations performed')

    # 日志选项|Logging options
    logging_grp = parser.add_argument_group('日志选项|Logging options')
    logging_grp.add_argument('-v', '--verbose', action='count', default=0,
                        help='增加输出详细程度 (-v, -vv, -vvv)|Increase output verbosity (-v, -vv, -vvv)')
    logging_grp.add_argument('--quiet', action='store_true',
                        help='静默模式，仅输出错误信息|Quiet mode, only output error messages')
    logging_grp.add_argument('--log-file', type=str,
                        help='日志文件路径|Log file path')
    logging_grp.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO',
                        help='日志级别|Log level')

    args = parser.parse_args()

    # 创建提取器并运行|Create extractor and run
    extractor = VCFGenotypeExtractor(
        vcf_file=args.vcf_file,
        output_prefix=args.output,
        samples=args.samples,
        biallelic_only=args.biallelic_only,
        min_length=args.min_length,
        max_length=args.max_length,
        split_by_chromosome=args.each,
        output_type=args.output_type,
        output_dir=args.output_dir,
        verbose=(args.verbose > 0),
        quiet=args.quiet,
        log_file=args.log_file,
        log_level=args.log_level,
        dry_run=args.dry_run
    )

    success = extractor.run_extraction()
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
