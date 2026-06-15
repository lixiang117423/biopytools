"""
NCBI分类学注释主程序模块|NCBI Taxonomy Annotation Main Module
"""

import argparse
import sys
from .config import NCBITaxoConfig
from .utils import NCBITaxoLogger, CommandRunner, check_dependencies
from .processor import BlastTaxonomyProcessor
from .stats import TaxonomyStatsCalculator


class NCBITaxoAnnotator:
    """NCBI分类学注释主类|Main NCBI Taxonomy Annotator Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = NCBITaxoConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = NCBITaxoLogger(self.config.output_prefix)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化处理器|Initialize processors
        self.processor = BlastTaxonomyProcessor(self.config, self.logger, self.cmd_runner)
        self.stats_calculator = TaxonomyStatsCalculator(self.config, self.logger)

    def run_analysis(self):
        """运行完整的分析流程|Run complete analysis pipeline"""
        try:
            # 检查依赖|Check dependencies
            self.logger.info("检查必需工具|Checking required tools")
            check_dependencies()
            self.logger.info("所有依赖检查通过|All dependencies checked")

            # 运行处理流程|Run processing pipeline
            accession_lineage, blast_hits, accessions_list = self.processor.run_pipeline()

            # 检查是否有有效lineage数据|Check if there is valid lineage data
            if not accession_lineage:
                self.logger.error("未获取到有效的lineage数据，无法生成统计|No valid lineage data retrieved, cannot generate statistics")
                return 1

            # 计算统计信息|Calculate statistics
            self.logger.info("开始计算统计信息|Starting statistics calculation")
            stats = self.stats_calculator.calculate_statistics(accession_lineage, blast_hits)

            if not stats:
                self.logger.error("统计计算失败|Statistics calculation failed")
                return 1

            # 写入统计文件|Write statistics file
            stats_file = self.stats_calculator.write_statistics(
                stats, accession_lineage, blast_hits
            )

            # 生成汇总报告|Generate summary report
            self.stats_calculator.generate_summary_report(
                accession_lineage, blast_hits, stats_file
            )

            self.logger.info("所有分析完成|All analyses completed successfully")
            return 0

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            import traceback
            self.logger.debug(f"错误详情|Error details:\n{traceback.format_exc()}")
            return 1


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='NCBI分类学注释工具|NCBI Taxonomy Annotation Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='输入文件路径（BLAST结果或accession列表）|Input file path (BLAST results or accession list)')
    parser.add_argument('-o', '--output-prefix', required=True,
                       help='输出文件前缀|Output file prefix')

    # 输入类型配置|Input type configuration
    parser.add_argument('--input-type',
                       choices=['auto', 'blast', 'accession'],
                       default='auto',
                       help='输入文件类型|Input file type (auto/blast/accession)')
    parser.add_argument('--blast-column',
                       type=int,
                       default=2,
                       help='BLAST结果中accession所在的列（从1开始）|Column containing accession in BLAST results (1-based)')
    parser.add_argument('-l', '--min-length',
                       type=int,
                       default=1000,
                       help='最小比对长度过滤（bp）|Minimum alignment length filter (bp)')
    parser.add_argument('--fetch-titles',
                       action='store_true',
                       help='获取accession的序列描述|Fetch accession titles (requires edirect)')


    # 数据库配置|Database configuration
    parser.add_argument('--taxid-db',
                       default='~/database/ncbi_taxonomy/nucl_gb.accession2taxid.gz',
                       help='TaxID数据库路径|TaxID database path')

    # 分类学配置|Taxonomy configuration
    parser.add_argument('--lineage-format',
                       default='{k};{p};{c};{o};{f};{g};{s}',
                       help='分类层级格式|Lineage format string')
    parser.add_argument('--no-full-lineage',
                       action='store_true',
                       help='不保留完整lineage|Do not keep full lineage')

    # 统计配置|Statistics configuration
    parser.add_argument('--stats-by',
                       nargs='+',
                       choices=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
                       default=['genus', 'species'],
                       help='统计层级|Statistics levels')
    parser.add_argument('--stats-target',
                       choices=['blast_hits', 'unique_accessions', 'both'],
                       default='both',
                       help='统计对象|Statistics target')
    parser.add_argument('--stats-output',
                       choices=['txt', 'csv'],
                       default='txt',
                       help='统计输出格式|Statistics output format')


    # 性能配置|Performance configuration
    parser.add_argument('-t', '--threads',
                       type=int,
                       default=4,
                       help='线程数|Number of threads')

    args = parser.parse_args()

    # 创建注释器并运行|Create annotator and run
    annotator = NCBITaxoAnnotator(
        input_file=args.input,
        output_prefix=args.output_prefix,
        input_type=args.input_type,
        blast_column=args.blast_column,
        min_alignment_length=args.min_length,
        fetch_titles=args.fetch_titles,
        taxid_db=args.taxid_db,
        lineage_format=args.lineage_format,
        keep_full_lineage=not args.no_full_lineage,
        stats_by=args.stats_by,
        stats_target=args.stats_target,
        stats_output=args.stats_output,
        threads=args.threads
    )

    exit_code = annotator.run_analysis()
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
