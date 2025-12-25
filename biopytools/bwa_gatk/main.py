"""
🧬 BWA-GATK流程主控制模块 | BWA-GATK Pipeline Main Control Module
"""

import sys
import argparse
from pathlib import Path

from .config import PipelineConfig
from .logger import PipelineLogger
from .utils import CommandRunner, check_dependencies
from .reference import ReferenceIndexer
from .sample_detector import SampleDetector
from .alignment import BWAAligner
from .preprocessing import BAMPreprocessor
from .variant_calling import VariantCaller
from .filtering import VariantFilter
from .statistics import StatisticsGenerator

class BWAGATKPipeline:
    """BWA-GATK流程主类 | BWA-GATK Pipeline Main Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = PipelineConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = PipelineLogger(
            self.config.logs_dir, 
            verbose=self.config.verbose
        )
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, dry_run=self.config.dry_run)
        
        # 初始化各个模块 | Initialize modules
        self.ref_indexer = ReferenceIndexer(self.config, self.logger, self.cmd_runner)
        self.sample_detector = SampleDetector(self.config, self.logger)
        self.aligner = BWAAligner(self.config, self.logger, self.cmd_runner)
        self.preprocessor = BAMPreprocessor(self.config, self.logger, self.cmd_runner)
        self.variant_caller = VariantCaller(self.config, self.logger, self.cmd_runner)
        self.filter = VariantFilter(self.config, self.logger, self.cmd_runner)
        self.stats_generator = StatisticsGenerator(self.config, self.logger, self.cmd_runner)
    
    def run(self):
        """运行完整流程 | Run complete pipeline"""
        try:
            self.logger.info("🚀 开始BWA-GATK变异检测流程 | Starting BWA-GATK variant calling pipeline")
            
            # 1. 检查依赖 | Check dependencies
            check_dependencies(self.config, self.logger)
            
            # 2. 准备参考基因组 | Prepare reference genome
            self.ref_indexer.check_and_build_indices()
            
            # 3. 检测样本 | Detect samples
            samples = self.sample_detector.detect_samples()
            
            # 4. 比对和预处理 | Alignment and preprocessing
            bam_files = []
            for sample in samples:
                # 生成Read Group | Generate Read Group
                read_group = self.sample_detector.get_read_group(sample['name'])
                
                # BWA比对 | BWA alignment
                sam_file = self.aligner.align_sample(sample, read_group)
                
                # BAM预处理 | BAM preprocessing
                bam_file = self.preprocessor.process_sample(sample['name'], sam_file)
                bam_files.append(bam_file)
            
            # 5. 变异检测（GVCF模式）| Variant calling (GVCF mode)
            raw_vcf = self.variant_caller.call_variants_gvcf_mode(samples, bam_files)
            
            # 6. 变异过滤 | Variant filtering
            filtered_vcfs = self.filter.filter_variants(raw_vcf)
            
            # 7. 生成统计 | Generate statistics
            self.stats_generator.generate_statistics(samples, bam_files, filtered_vcfs)
            
            # 完成 | Complete
            self.logger.info("=" * 80)
            self.logger.info("🎉 流程完成！| Pipeline completed successfully!")
            self.logger.info("=" * 80)
            self.logger.info(f"📁 输出目录 | Output directory: {self.config.output_dir}")
            self.logger.info(f"📊 统计报告 | Statistics report: {self.config.stats_dir / 'pipeline_statistics.txt'}")
            self.logger.info(f"📝 日志文件 | Log file: {self.logger_manager.log_file}")
            
        except Exception as e:
            self.logger.error(f"❌ 流程执行失败 | Pipeline execution failed: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 BWA-GATK变异检测流程 (模块化版本) | BWA-GATK Variant Calling Pipeline (Modular Version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  # 基本用法 | Basic usage
  %(prog)s -r reference.fa -i fastq_dir/ -o output/
  
  # 指定倍性和线程 | Specify ploidy and threads
  %(prog)s -r ref.fa -i fastq/ -o out/ -p 4 -t 64
  
  # 使用区间限定和详细日志 | Use intervals and verbose logging
  %(prog)s -r ref.fa -i fastq/ -o out/ --intervals target.bed -v
  
  # 干运行模式（仅显示命令）| Dry-run mode (show commands only)
  %(prog)s -r ref.fa -i fastq/ -o out/ --dry-run
  
  # 强制重新开始 | Force restart
  %(prog)s -r ref.fa -i fastq/ -o out/ --force-restart
        """
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('必需参数 | Required arguments')
    required.add_argument('-r', '--reference', required=True,
                         help='参考基因组FASTA文件 | Reference genome FASTA file')
    required.add_argument('-i', '--input', required=True, dest='input_path',
                         help='输入FASTQ文件或目录 | Input FASTQ file or directory')
    
    # 输出参数 | Output arguments
    output = parser.add_argument_group('输出参数 | Output arguments')
    output.add_argument('-o', '--output', default='./bwa_gatk_output', dest='output_dir',
                       help='输出目录 (默认: ./bwa_gatk_output) | Output directory (default: ./bwa_gatk_output)')
    
    # 分析参数 | Analysis arguments
    analysis = parser.add_argument_group('分析参数 | Analysis arguments')
    analysis.add_argument('-t', '--threads', type=int, default=88,
                         help='线程数 (默认: 88) | Number of threads (default: 88)')
    analysis.add_argument('-p', '--ploidy', type=int, default=2,
                         help='倍性 (默认: 2) | Ploidy (default: 2)')
    analysis.add_argument('-m', '--mem-per-thread', type=int, default=10,
                         help='每线程内存(GB) (默认: 10) | Memory per thread in GB (default: 10)')
    analysis.add_argument('--intervals',
                         help='区间限定BED文件 | Intervals BED file for targeted analysis')
    
    # 过滤参数 | Filtering arguments (可选，使用默认GATK Best Practices)
    filtering = parser.add_argument_group('过滤参数 | Filtering arguments (optional, GATK Best Practices used by default)')
    filtering.add_argument('--snp-qd', type=float, default=2.0,
                          help='SNP QD阈值 (默认: 2.0) | SNP QD threshold (default: 2.0)')
    filtering.add_argument('--snp-fs', type=float, default=60.0,
                          help='SNP FS阈值 (默认: 60.0) | SNP FS threshold (default: 60.0)')
    filtering.add_argument('--indel-qd', type=float, default=2.0,
                          help='InDel QD阈值 (默认: 2.0) | InDel QD threshold (default: 2.0)')
    
    # 运行选项 | Run options
    run_opts = parser.add_argument_group('运行选项 | Run options')
    run_opts.add_argument('--force-restart', action='store_true',
                         help='强制从头开始，忽略已有文件 | Force restart, ignore existing files')
    run_opts.add_argument('--dry-run', action='store_true',
                         help='干运行模式，仅显示命令不执行 | Dry-run mode, show commands without execution')
    run_opts.add_argument('-v', '--verbose', action='store_true',
                         help='详细日志模式 | Verbose logging mode')
    
    # 工具路径 | Tool paths
    tools = parser.add_argument_group('工具路径 | Tool paths (optional)')
    tools.add_argument('--bwa-path', default='bwa',
                      help='BWA路径 (默认: bwa) | BWA path (default: bwa)')
    tools.add_argument('--samtools-path', default='samtools',
                      help='SAMtools路径 (默认: samtools) | SAMtools path (default: samtools)')
    tools.add_argument('--gatk-path', default='gatk',
                      help='GATK路径 (默认: gatk) | GATK path (default: gatk)')
    
    args = parser.parse_args()
    
    # 创建并运行流程 | Create and run pipeline
    pipeline = BWAGATKPipeline(
        reference=args.reference,
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
        indel_qd=args.indel_qd
    )
    
    pipeline.run()

if __name__ == "__main__":
    main()
