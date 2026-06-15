"""
BWA-GATK流程主控制模块|BWA-GATK Pipeline Main Control Module
"""

import sys
import argparse
from pathlib import Path

from .config import PipelineConfig
from .logger import PipelineLogger
from .utils import CommandRunner, check_dependencies
from .reference import ReferenceIndexer
from .sample_detector import SampleDetector
from .quality_control import QualityController
from .alignment import BWAAligner
from .preprocessing import BAMPreprocessor
from .variant_calling import VariantCaller
from .filtering import VariantFilter
from .statistics import StatisticsGenerator

class BWAGATKPipeline:
    """BWA-GATK流程主类|BWA-GATK Pipeline Main Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = PipelineConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志|Initialize logging
        self.logger_manager = PipelineLogger(
            self.config.logs_dir, 
            verbose=self.config.verbose
        )
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, dry_run=self.config.dry_run)
        
        # 初始化各个模块|Initialize modules
        self.ref_indexer = ReferenceIndexer(self.config, self.logger, self.cmd_runner)
        self.quality_controller = QualityController(self.config, self.logger, self.cmd_runner)
        self.sample_detector = SampleDetector(self.config, self.logger)
        self.aligner = BWAAligner(self.config, self.logger, self.cmd_runner)
        self.preprocessor = BAMPreprocessor(self.config, self.logger, self.cmd_runner)
        self.variant_caller = VariantCaller(self.config, self.logger, self.cmd_runner)
        self.filter = VariantFilter(self.config, self.logger, self.cmd_runner)
        self.stats_generator = StatisticsGenerator(self.config, self.logger, self.cmd_runner)
    
    def run(self):
        """运行完整流程|Run complete pipeline"""
        try:
            self.logger.info("开始BWA-GATK变异检测流程|Starting BWA-GATK variant calling pipeline")

            # 1. 检查依赖|Check dependencies
            check_dependencies(self.config, self.logger)

            # 2. 质量控制（如果启用）|Quality control (if enabled)
            if not self.config.skip_qc:
                qc_success = self.quality_controller.run_quality_control()
                if not qc_success:
                    self.logger.error("质控失败|QC failed")
                    sys.exit(1)

            # 3. 准备参考基因组|Prepare reference genome
            self.ref_indexer.check_and_build_indices()

            # 4. 检测样本|Detect samples
            # 使用清洁FASTQ目录（如果质控已执行）或原始输入目录
            detection_path = self.config.clean_fastq_dir if not self.config.skip_qc else self.config.input_path
            samples = self.sample_detector.detect_samples()

            # 5. 比对和预处理|Alignment and preprocessing
            bam_files = []
            for sample in samples:
                # 生成Read Group|Generate Read Group
                read_group = self.sample_detector.get_read_group(sample['name'])

                # BWA比对|BWA alignment
                sam_file = self.aligner.align_sample(sample, read_group)

                # BAM预处理|BAM preprocessing
                bam_file = self.preprocessor.process_sample(sample['name'], sam_file)
                bam_files.append(bam_file)

            # 6. 变异检测（GVCF模式）| Variant calling (GVCF mode)
            raw_vcf = self.variant_caller.call_variants_gvcf_mode(samples, bam_files)

            # 7. 变异过滤|Variant filtering
            filtered_vcfs = self.filter.filter_variants(raw_vcf)

            # 8. 生成统计|Generate statistics
            self.stats_generator.generate_statistics(samples, bam_files, filtered_vcfs)

            # 完成|Complete
            self.logger.info("=" * 80)
            self.logger.info("流程完成|Pipeline completed successfully!")
            self.logger.info("=" * 80)
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            if not self.config.skip_qc:
                self.logger.info(f"清洁数据目录|Clean data directory: {self.config.clean_fastq_dir}")
            self.logger.info(f"统计报告|Statistics report: {self.config.stats_dir / 'pipeline_statistics.txt'}")
            self.logger.info(f"日志文件|Log file: {self.logger_manager.log_file}")

        except Exception as e:
            self.logger.error(f"流程执行失败|Pipeline execution failed: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            sys.exit(1)

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='BWA-GATK变异检测流程(模块化版本)|BWA-GATK Variant Calling Pipeline (Modular Version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples: %(prog)s -i fastq_dir/ -g reference.fa -o output/
        '''
    )
    
    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-g', '--genome', required=True,
                         help='参考基因组FASTA文件|Reference genome FASTA file')
    required.add_argument('-i', '--input', required=True, dest='input_path',
                         help='输入FASTQ文件目录（包含原始或清洁的FASTQ文件）|Input FASTQ directory (containing raw or clean FASTQ files)')

    # 输出参数|Output arguments
    output = parser.add_argument_group('输出参数|Output arguments')
    output.add_argument('-o', '--output-dir', default='./bwa_gatk_output', dest='output_dir',
                       help='输出目录|Output directory')
    
    # 分析参数|Analysis arguments
    analysis = parser.add_argument_group('分析参数|Analysis arguments')
    analysis.add_argument('-t', '--threads', type=int, default=88,
                         help='线程数|Number of threads')
    analysis.add_argument('-p', '--ploidy', type=int, default=2,
                         help='倍性|Ploidy')
    analysis.add_argument('-m', '--mem-per-thread', type=int, default=10,
                         help='每线程内存(GB)|Memory per thread in GB')
    analysis.add_argument('--intervals',
                         help='区间限定BED文件|Intervals BED file for targeted analysis')

    # 过滤参数|Filtering arguments (optional, GATK Best Practices used by default)
    filtering = parser.add_argument_group('过滤参数|Filtering arguments (optional, GATK Best Practices used by default)')
    filtering.add_argument('--snp-qd', type=float, default=2.0,
                          help='SNP QD阈值|SNP QD threshold')
    filtering.add_argument('--snp-fs', type=float, default=60.0,
                          help='SNP FS阈值|SNP FS threshold')
    filtering.add_argument('--indel-qd', type=float, default=2.0,
                          help='InDel QD阈值|InDel QD threshold')
    
    # 运行选项|Run options
    run_opts = parser.add_argument_group('运行选项|Run options')
    run_opts.add_argument('--force-restart', action='store_true',
                         help='强制从头开始，忽略已有文件|Force restart, ignore existing files')
    run_opts.add_argument('--dry-run', action='store_true',
                         help='干运行模式，仅显示命令不执行|Dry-run mode, show commands without execution')
    run_opts.add_argument('-v', '--verbose', action='store_true',
                         help='详细日志模式|Verbose logging mode')

    # 质控选项|Quality control options
    qc_opts = parser.add_argument_group('质控选项|Quality control options')
    qc_opts.add_argument('--skip-qc', action='store_true',
                        help='跳过质控步骤，直接使用输入数据|Skip QC step, use input data directly')
    qc_opts.add_argument('--fastp-path', default='fastp',
                        help='fastp可执行文件路径|fastp executable path')
    qc_opts.add_argument('--qc-threads', type=int, default=12,
                        help='质控线程数|QC threads')
    qc_opts.add_argument('--qc-quality-threshold', type=int, default=20,
                        help='质控质量阈值|QC quality threshold')
    qc_opts.add_argument('--qc-min-length', type=int, default=50,
                        help='质控最小长度|QC minimum length')
    qc_opts.add_argument('--qc-unqualified-percent', type=int, default=40,
                        help='不合格碱基百分比阈值|Unqualified base percentage threshold')
    qc_opts.add_argument('--qc-n-base-limit', type=int, default=10,
                        help='N碱基数量限制|N base count limit')
    qc_opts.add_argument('--qc-read1-suffix',
                        help='R1文件后缀|R1 file suffix')
    qc_opts.add_argument('--qc-read2-suffix',
                        help='R2文件后缀|R2 file suffix')
    qc_opts.add_argument('--qc-single-end', action='store_true',
                        help='单末端测序模式|Single-end sequencing mode')

    # 工具路径|Tool paths (optional)
    tools = parser.add_argument_group('工具路径|Tool paths (optional)')
    tools.add_argument('--bwa-path', default='bwa',
                      help='BWA路径|BWA path')
    tools.add_argument('--samtools-path', default='samtools',
                      help='SAMtools路径|SAMtools path')
    tools.add_argument('--gatk-path', default='gatk',
                      help='GATK路径|GATK path')
    
    args = parser.parse_args()

    # 创建并运行流程|Create and run pipeline
    pipeline = BWAGATKPipeline(
        reference=args.genome,
        input_path=args.input_path,
        output_dir=args.output_dir,
        threads=args.threads,
        ploidy=args.ploidy,
        mem_per_thread=args.mem_per_thread,
        intervals=args.intervals,
        force_restart=args.force_restart,
        dry_run=args.dry_run,
        verbose=args.verbose,
        bwa_path=args.bwa_path,
        samtools_path=args.samtools_path,
        gatk_path=args.gatk_path,
        snp_qd=args.snp_qd,
        snp_fs=args.snp_fs,
        indel_qd=args.indel_qd,
        # 质控参数|QC parameters
        skip_qc=args.skip_qc,
        fastp_path=args.fastp_path,
        qc_threads=args.qc_threads,
        qc_quality_threshold=args.qc_quality_threshold,
        qc_min_length=args.qc_min_length,
        qc_unqualified_percent=args.qc_unqualified_percent,
        qc_n_base_limit=args.qc_n_base_limit,
        qc_read1_suffix=args.qc_read1_suffix,
        qc_read2_suffix=args.qc_read2_suffix,
        qc_single_end=args.qc_single_end
    )
    
    pipeline.run()

if __name__ == "__main__":
    main()
