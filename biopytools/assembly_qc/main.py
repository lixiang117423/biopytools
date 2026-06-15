"""
基因组组装质量评估主程序模块|Genome Assembly Quality Control Main Module
"""

import sys
import argparse
from pathlib import Path
from .config import AssemblyQCConfig
from .utils import AssemblyQCLogger, CommandRunner, calculate_genome_stats
from .busco_evaluator import BUSCOEvaluator
from .lai_evaluator import LAIEvaluator
from .qv_evaluator import QVEvaluator
from .mapping_evaluator import MappingEvaluator
from .long_read_mapping_evaluator import LongReadMappingEvaluator
from .report_generator import ReportGenerator


class AssemblyQC:
    """基因组组装质量评估主类|Genome Assembly Quality Control Main Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = AssemblyQCConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = AssemblyQCLogger(Path(self.config.output_dir))
        self.logger = self.logger_manager.get_logger()

        # 初始化评估器|Initialize evaluators (不再使用cmd_runner|No longer using cmd_runner)
        self.busco_evaluator = BUSCOEvaluator(self.config, self.logger)
        self.lai_evaluator = LAIEvaluator(self.config, self.logger)
        self.qv_evaluator = QVEvaluator(self.config, self.logger)
        self.mapping_evaluator = MappingEvaluator(self.config, self.logger)
        self.long_read_mapping_evaluator = LongReadMappingEvaluator(self.config, self.logger)

        # 初始化报告生成器|Initialize report generator
        self.report_generator = ReportGenerator(self.config, self.logger)

    def run_analysis(self):
        """运行完整的质量评估流程|Run complete quality control pipeline"""
        try:
            self._print_header()
            self._print_config()

            # 存储所有评估结果|Store all evaluation results
            results = {}

            # 跟踪步骤状态|Track step status
            step_status = {}

            # 1. 基因组统计|Genome statistics
            self.logger.info("=" * 80)
            self.logger.info("步骤1: 基因组统计|Step 1: Genome Statistics")
            self.logger.info("=" * 80)
            results['genome_stats'] = calculate_genome_stats(self.config.genome, self.logger)

            # 2. BUSCO评估|BUSCO evaluation
            self.logger.info("=" * 80)
            self.logger.info("步骤2: BUSCO完整性评估|Step 2: BUSCO Completeness Evaluation")
            self.logger.info("=" * 80)
            if self.config.skip_busco:
                step_status['busco'] = 'skipped'
                self.logger.info("跳过BUSCO评估|Skipping BUSCO evaluation (user specified)")
            else:
                busco_results = self.busco_evaluator.evaluate()
                if busco_results:
                    results['busco'] = busco_results
                    step_status['busco'] = 'success'
                else:
                    step_status['busco'] = 'failed'

            # 3. LAI评估|LAI evaluation
            self.logger.info("=" * 80)
            self.logger.info("步骤3: LAI指数评估|Step 3: LAI Index Evaluation")
            self.logger.info("=" * 80)
            if self.config.skip_lai:
                step_status['lai'] = 'skipped'
                self.logger.info("跳过LAI评估|Skipping LAI evaluation (user specified)")
            else:
                lai_results = self.lai_evaluator.evaluate()
                if lai_results:
                    results['lai'] = lai_results
                    step_status['lai'] = 'success'
                else:
                    step_status['lai'] = 'failed'

            # 4. QV评估（可选）|QV evaluation (optional)
            if self.config.enable_qv:
                self.logger.info("=" * 80)
                self.logger.info("步骤4: QV质量值评估|Step 4: QV Quality Value Evaluation")
                self.logger.info("=" * 80)
                qv_results = self.qv_evaluator.evaluate()
                if qv_results:
                    results['qv'] = qv_results
                    step_status['qv'] = 'success'
                else:
                    step_status['qv'] = 'failed'
            else:
                step_status['qv'] = 'disabled'

            # 5. Mapping评估（可选）|Mapping evaluation (optional)
            if self.config.enable_mapping:
                self.logger.info("=" * 80)
                self.logger.info("步骤5: Mapping评估|Step 5: Mapping Evaluation")
                self.logger.info("=" * 80)
                mapping_results = self.mapping_evaluator.evaluate()
                if mapping_results:
                    results['mapping'] = mapping_results
                    step_status['mapping'] = 'success'
                else:
                    step_status['mapping'] = 'failed'
            else:
                step_status['mapping'] = 'disabled'

            # 6. 三代数据Mapping评估（可选）|Long-read mapping evaluation (optional)
            if self.config.enable_long_read_mapping:
                self.logger.info("=" * 80)
                self.logger.info("步骤6: 三代数据Mapping评估|Step 6: Long-read Mapping Evaluation")
                self.logger.info("=" * 80)
                long_read_mapping_results = self.long_read_mapping_evaluator.evaluate()
                if long_read_mapping_results:
                    results['long_read_mapping'] = long_read_mapping_results
                    step_status['long_read_mapping'] = 'success'
                else:
                    step_status['long_read_mapping'] = 'failed'
            else:
                step_status['long_read_mapping'] = 'disabled'

            # 7. 生成报告|Generate reports
            self.logger.info("=" * 80)
            self.logger.info("步骤7: 生成报告|Step 7: Generate Reports")
            self.logger.info("=" * 80)
            self.report_generator.generate_reports(results)

            # 7. 打印摘要|Print summary
            self._print_summary(results, step_status)

            self.logger.info("=" * 80)
            self.logger.info("基因组组装质量评估完成|Genome assembly quality control completed")
            self.logger.info("=" * 80)

            return True

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False

    def _print_header(self):
        """打印标题|Print header"""
        self.logger.info("")
        self.logger.info("=" * 80)
        self.logger.info("基因组组装质量评估流程|Genome Assembly Quality Control Workflow")
        self.logger.info("=" * 80)
        self.logger.info("")

    def _print_config(self):
        """打印配置信息|Print configuration"""
        self.logger.info("工作流配置|Workflow Configuration:")
        self.logger.info("-" * 80)
        summary = self.config.get_summary()

        self.logger.info(f"  样品名称|Sample name: {summary['sample_name']}")
        self.logger.info(f"  基因组|Genome: {summary['genome']}")
        self.logger.info(f"  BUSCO谱系|BUSCO lineage: {summary['lineage']}")
        self.logger.info("")
        self.logger.info("流程控制|Workflow control:")
        self.logger.info(f"  BUSCO评估|BUSCO evaluation: {'启用|Enabled' if summary['busco_enabled'] else '禁用|Disabled'}")
        self.logger.info(f"  LAI评估|LAI evaluation: {'启用|Enabled' if summary['lai_enabled'] else '禁用|Disabled'}")
        self.logger.info(f"  QV评估|QV evaluation: {'启用|Enabled' if summary['qv_enabled'] else '禁用|Disabled'}")
        self.logger.info(f"  Mapping评估|Mapping evaluation: {'启用|Enabled' if summary['mapping_enabled'] else '禁用|Disabled'}")
        self.logger.info(f"  三代数据Mapping评估|Long-read mapping evaluation: {'启用|Enabled' if self.config.enable_long_read_mapping else '禁用|Disabled'}")
        self.logger.info("")
        self.logger.info(f"  输出目录|Output directory: {summary['output_dir']}")
        self.logger.info("")

    def _print_summary(self, results: dict, step_status: dict):
        """打印评估摘要|Print evaluation summary"""
        self.logger.info("")
        self.logger.info("=" * 80)
        self.logger.info("评估摘要|Evaluation Summary")
        self.logger.info("=" * 80)

        # 基因组统计|Genome statistics
        genome_stats = results.get('genome_stats', {})
        self.logger.info(f"基因组大小|Genome size: {genome_stats.get('total_size_mb', 'N/A'):.2f} Mb")
        self.logger.info(f"Contig N50: {genome_stats.get('contig_n50_mb', 'N/A'):.2f} Mb")
        self.logger.info(f"GC含量|GC content: {genome_stats.get('gc_content', 'N/A'):.2f}%")

        # 步骤状态|Step status
        self.logger.info("")
        self.logger.info("评估步骤状态|Evaluation Step Status:")
        self.logger.info("-" * 80)

        status_text = {
            'success': '成功|Success',
            'failed': '失败|Failed',
            'skipped': '跳过|Skipped',
            'disabled': '未启用|Disabled'
        }

        for step, status in step_status.items():
            status_msg = status_text.get(status, status)
            step_name = {
                'busco': 'BUSCO完整性评估|BUSCO Completeness',
                'lai': 'LAI指数评估|LAI Index',
                'qv': 'QV质量值评估|QV Quality Value',
                'mapping': 'Mapping评估|Mapping Evaluation',
                'long_read_mapping': '三代数据Mapping评估|Long-read Mapping Evaluation'
            }.get(step, step)
            self.logger.info(f"  {step_name}: {status_msg}")

        # 显示成功步骤的结果|Display results of successful steps
        if step_status.get('busco') == 'success':
            busco_results = results.get('busco', {})
            self.logger.info(f"BUSCO完整度|BUSCO completeness: {busco_results.get('complete', 'N/A')}%")

        if step_status.get('lai') == 'success':
            lai_results = results.get('lai', {})
            if lai_results.get('error'):
                self.logger.warning(f"LAI指数|LAI index: 不适用|Not applicable - {lai_results.get('error', 'Unknown reason')}")
            else:
                self.logger.info(f"LAI指数|LAI index: {lai_results.get('lai_score', 'N/A')}")

        if step_status.get('qv') == 'success':
            qv_results = results.get('qv', {})
            # 支持多种数据类型的QV值和错误率|Support QV values and error rates for multiple data types
            if 'ngs_qv_value' in qv_results:
                error_rate = qv_results.get('ngs_error_rate', 'N/A')
                error_str = f"{error_rate:.2e}" if error_rate != 'N/A' else 'N/A'
                self.logger.info(f"NGS数据QV值|NGS QV value: {qv_results.get('ngs_qv_value', 'N/A')} (错误率|Error rate: {error_str})")
            if 'long_read_qv_value' in qv_results:
                error_rate = qv_results.get('long_read_error_rate', 'N/A')
                error_str = f"{error_rate:.2e}" if error_rate != 'N/A' else 'N/A'
                self.logger.info(f"三代数据QV值|Long-read QV value: {qv_results.get('long_read_qv_value', 'N/A')} (错误率|Error rate: {error_str})")
            # 如果只有一个qv_value（向后兼容）|If only qv_value exists (backward compatibility)
            if 'qv_value' in qv_results and 'ngs_qv_value' not in qv_results and 'long_read_qv_value' not in qv_results:
                self.logger.info(f"QV值|QV value: {qv_results.get('qv_value', 'N/A')}")

        if step_status.get('mapping') == 'success':
            mapping_results = results.get('mapping', {})
            mapping_rate = mapping_results.get('overall_mapping_rate', mapping_results.get('mapping_rate', 'N/A'))
            mean_coverage = mapping_results.get('mean_coverage', 'N/A')
            self.logger.info(f"比对率|Mapping rate: {mapping_rate}")
            self.logger.info(f"平均覆盖度|Mean coverage: {mean_coverage}")

        if step_status.get('long_read_mapping') == 'success':
            long_read_mapping_results = results.get('long_read_mapping', {})
            overall_mapping_rate = long_read_mapping_results.get('overall_mapping_rate', 'N/A')
            total_reads = long_read_mapping_results.get('total_reads', 'N/A')
            mapped_reads = long_read_mapping_results.get('mapped_reads', 'N/A')
            self.logger.info(f"三代数据总reads数|Long-read total reads: {total_reads}")
            self.logger.info(f"三代数据比对reads数|Long-read mapped reads: {mapped_reads}")
            self.logger.info(f"三代数据比对率|Long-read mapping rate: {overall_mapping_rate}")

        self.logger.info("")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info("")


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="基因组组装质量评估工具|Genome Assembly Quality Control Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Example:
  %(prog)s --genome genome.fa --lineage embryophyta_odb10 --ngs-reads ./illumina_reads --long-reads ./hifi_reads --long-read-type hifi -o qc_results
        '''
    )

    # 必需参数|Required parameters
    parser.add_argument('--genome', '-g',
                        required=True,
                        help='基因组FASTA文件|Genome FASTA file')

    parser.add_argument('--lineage', '-l',
                        required=True,
                        help='BUSCO数据集名称（如embryophyta_odb10）或完整路径|BUSCO dataset name (e.g., embryophyta_odb10) or full path')

    parser.add_argument('--output-dir', '-o',
                        default='./assembly_qc_output',
                        help='输出目录|Output directory')

    # 样品信息|Sample information
    parser.add_argument('--sample-name', '-s',
                        default='genome_sample',
                        help='样品名称|Sample name')

    # 核心评估参数|Core evaluation parameters
    parser.add_argument('--skip-busco',
                        action='store_true',
                        help='跳过BUSCO评估|Skip BUSCO evaluation')

    parser.add_argument('--skip-lai',
                        action='store_true',
                        help='跳过LAI评估|Skip LAI evaluation')

    parser.add_argument('--lai-full-mode',
                        action='store_true',
                        help='LAI完整模式（不使用-qq，运行blastn计算，用于种间比较）|LAI full mode (no -qq, run blastn for interspecies comparison)')

    # QV评估参数|QV evaluation parameters
    parser.add_argument('--enable-qv',
                        action='store_true',
                        help='启用QV评估（默认启用）|Enable QV evaluation (default: enabled)')

    parser.add_argument('--qv-kmer-size',
                        type=int,
                        help='k-mer大小（None表示自动选择）|K-mer size (None for auto)')

    # Mapping评估参数|Mapping evaluation parameters
    parser.add_argument('--enable-mapping',
                        action='store_true',
                        help='启用NGS Mapping评估（默认启用）|Enable NGS mapping evaluation (default: enabled)')

    parser.add_argument('--enable-long-read-mapping',
                        action='store_true',
                        help='启用三代数据Mapping评估（默认启用）|Enable long-read mapping evaluation (default: enabled)')

    parser.add_argument('--mapping-pattern',
                        default='_1.clean.fq.gz',
                        help='FASTQ文件匹配模式|FASTQ file pattern')

    # Reads数据|Reads data
    parser.add_argument('--ngs-reads',
                        help='NGS reads目录（用于QV和mapping）|NGS reads directory (for QV and mapping)')

    parser.add_argument('--long-reads',
                        help='Long-reads目录（用于QV和mapping）|Long-reads directory (for QV and mapping)')

    parser.add_argument('--long-read-type',
                        choices=['ont', 'pacbio', 'hifi'],
                        default='hifi',
                        help='Long-read数据类型|Long-read data type')

    # 报告参数|Report parameters
    parser.add_argument('--no-html',
                        action='store_true',
                        help='不生成HTML报告|Do not generate HTML report')

    parser.add_argument('--no-table',
                        action='store_true',
                        help='不生成表格|Do not generate table')

    parser.add_argument('--table-format',
                        choices=['tsv', 'xlsx', 'both'],
                        default='both',
                        help='表格格式|Table format')

    # 线程参数|Thread parameter
    parser.add_argument('--threads', '-t',
                        type=int,
                        default=12,
                        help='线程数（自动分配给各子模块）|Threads (automatically distributed to sub-modules)')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    try:
        # 解析参数|Parse arguments
        args = parse_arguments()

        # 构建配置字典|Build configuration dictionary
        config = {
            'genome': args.genome,
            'lineage': args.lineage,
            'output_dir': args.output_dir,
            'sample_name': args.sample_name,
            'threads': args.threads,

            # 核心评估|Core evaluation
            'skip_busco': args.skip_busco,
            'skip_lai': args.skip_lai,
            'lai_quick_mode': not args.lai_full_mode,  # 默认True使用快速模式，--lai-full-mode设为False

            # QV评估|QV evaluation
            'enable_qv': args.enable_qv,
            'qv_kmer_size': args.qv_kmer_size,

            # Mapping评估|Mapping evaluation
            'enable_mapping': args.enable_mapping,
            'mapping_pattern': args.mapping_pattern,

            # 三代数据Mapping评估|Long-read mapping evaluation
            'enable_long_read_mapping': args.enable_long_read_mapping,

            # Reads数据|Reads data
            'ngs_reads': args.ngs_reads,
            'long_reads': args.long_reads,
            'long_read_type': args.long_read_type,

            # 报告参数|Report parameters
            'generate_html': not args.no_html,
            'generate_table': not args.no_table,
            'table_format': args.table_format,
        }

        # 创建并运行评估器|Create and run evaluator
        qc = AssemblyQC(**config)
        success = qc.run_analysis()

        sys.exit(0 if success else 1)

    except KeyboardInterrupt:
        print("\n用户中断操作|User interrupted operation")
        sys.exit(1)
    except Exception as e:
        print(f"程序错误|Program error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
