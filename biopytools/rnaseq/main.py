# """
# RNA-seq分析主程序模块 | RNA-seq Analysis Main Module
# """

# import argparse
# import sys
# import os
# from .config import RNASeqConfig
# from .utils import RNASeqLogger, CommandRunner
# from .data_processing import SampleParser
# from .alignment import HISAT2Indexer, HISAT2Aligner
# from .quantification import StringTieQuantifier, GTFValueExtractor
# from .results import ExpressionMatrixMerger, SummaryGenerator

# class RNASeqAnalyzer:
#     """RNA-seq分析主类 | Main RNA-seq Analyzer Class"""
    
#     def __init__(self, **kwargs):
#         # 初始化配置 | Initialize configuration
#         self.config = RNASeqConfig(**kwargs)
#         self.config.validate()
        
#         # 初始化日志 | Initialize logging
#         self.logger_manager = RNASeqLogger(self.config.output_path)
#         self.logger = self.logger_manager.get_logger()
        
#         # 初始化命令执行器 | Initialize command runner
#         self.cmd_runner = CommandRunner(self.logger)
        
#         # 初始化各个处理器 | Initialize processors
#         self.sample_parser = SampleParser(self.logger)
#         self.hisat2_indexer = HISAT2Indexer(self.config, self.logger, self.cmd_runner)
#         self.hisat2_aligner = HISAT2Aligner(self.config, self.logger, self.cmd_runner)
#         self.stringtie_quantifier = StringTieQuantifier(self.config, self.logger, self.cmd_runner)
#         self.gtf_extractor = GTFValueExtractor(self.config, self.logger)
#         self.matrix_merger = ExpressionMatrixMerger(self.config, self.logger)
#         self.summary_generator = SummaryGenerator(self.config, self.logger)
    
#     def process_single_sample(self, sample_info: dict, index_prefix: str) -> str:
#         """处理单个样本的完整工作流程 | Process complete workflow for a single sample"""
#         sample_name = sample_info["name"]
#         fastq1 = sample_info["fastq1"]
#         fastq2 = sample_info["fastq2"]

#         self.logger.info(f"\n{'=' * 60}")
#         self.logger.info(f"处理样本 | Processing sample: {sample_name}")
#         self.logger.info(f"{'=' * 60}")

#         # 设置文件路径 | Set file paths
#         bam_file = os.path.join(self.config.output_dir, f"{sample_name}.sorted.bam")
#         stringtie_output = os.path.join(self.config.output_dir, "stringtie_output", f"{sample_name}.gtf")
#         fpkm_output = os.path.join(self.config.output_dir, "fpkm_output", f"{sample_name}.fpkm.txt")

#         # 1. HISAT2比对 | HISAT2 alignment
#         if not self.hisat2_aligner.run_hisat2_mapping(index_prefix, fastq1, fastq2, bam_file):
#             return None

#         # 2. StringTie定量 | StringTie quantification
#         if not self.stringtie_quantifier.run_stringtie(bam_file, stringtie_output):
#             return None

#         # 3. 提取FPKM值 | Extract FPKM values
#         if not self.gtf_extractor.extract_gtf_values(stringtie_output, sample_name, fpkm_output):
#             return None

#         # 4. 如果请求，删除BAM文件 | Remove BAM file if requested
#         if self.config.remove_bam.lower() in ["yes", "y"]:
#             if os.path.exists(bam_file):
#                 os.remove(bam_file)
#                 self.logger.info(f"✓ 已删除BAM文件 | Removed BAM file: {bam_file}")
#         else:
#             self.logger.info(f"✓ 保留BAM文件 | BAM file retained: {bam_file}")

#         return fpkm_output
    
#     def run_analysis(self):
#         """运行完整的RNA-seq分析流程 | Run complete RNA-seq analysis pipeline"""
#         try:
#             self.logger.info("RNA-seq分析流程开始 | RNA-seq analysis pipeline started")
#             self.logger.info(f"基因组文件 | Genome file: {self.config.genome_file}")
#             self.logger.info(f"输入路径 | Input path: {self.config.input_path}")
#             self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
#             self.logger.info(f"GTF文件 | GTF file: {self.config.gtf_file}")
#             self.logger.info(f"线程数 | Threads: {self.config.threads}")
#             self.logger.info(f"删除BAM文件 | Remove BAM files: {self.config.remove_bam}")
#             if self.config.fastq_pattern:
#                 self.logger.info(f"文件模式 | File pattern: {self.config.fastq_pattern}")

#             # 步骤1: 构建HISAT2索引 | Step 1: Build HISAT2 index
#             self.logger.info(f"\n{'=' * 60}")
#             self.logger.info("步骤1: 构建HISAT2索引 | Step 1: Building HISAT2 index")
#             self.logger.info(f"{'=' * 60}")
#             index_prefix = self.hisat2_indexer.build_hisat2_index()

#             # 步骤2: 解析输入样本 | Step 2: Parse input samples
#             self.logger.info(f"\n{'=' * 60}")
#             self.logger.info("步骤2: 解析输入样本 | Step 2: Parsing input samples")
#             self.logger.info(f"{'=' * 60}")
#             samples = self.sample_parser.parse_input_samples(self.config.input_path, self.config.fastq_pattern)
            
#             if not samples:
#                 self.logger.error("错误：未找到有效的样本文件 | Error: No valid sample files found")
#                 sys.exit(1)

#             self.logger.info(f"找到 {len(samples)} 个样本 | Found {len(samples)} samples:")
#             for sample in samples:
#                 self.logger.info(f"  - {sample['name']}: {sample['fastq1']}, {sample['fastq2']}")
            
#             # 保存样本信息到配置 | Save sample info to config
#             self.config.samples = samples

#             # 步骤3: 处理所有样本 | Step 3: Process all samples
#             self.logger.info(f"\n{'=' * 60}")
#             self.logger.info("步骤3: 处理所有样本 | Step 3: Processing all samples")
#             self.logger.info(f"{'=' * 60}")
#             fpkm_files = []

#             for sample_info in samples:
#                 fpkm_file = self.process_single_sample(sample_info, index_prefix)
#                 if fpkm_file:
#                     fpkm_files.append(fpkm_file)
#                 else:
#                     self.logger.error(f"样本处理失败 | Sample processing failed: {sample_info['name']}")

#             # 步骤4: 合并表达矩阵 | Step 4: Merge expression matrix
#             self.logger.info(f"\n{'=' * 60}")
#             self.logger.info("步骤4: 合并表达矩阵 | Step 4: Merging expression matrix")
#             self.logger.info(f"{'=' * 60}")
            
#             if not self.matrix_merger.merge_expression_matrix(fpkm_files):
#                 self.logger.error("表达矩阵合并失败 | Expression matrix merging failed")
#                 sys.exit(1)

#             # 生成总结报告 | Generate summary report
#             self.summary_generator.generate_summary_report()

#             self.logger.info(f"\n{'=' * 60}")
#             self.logger.info("分析完成！| Analysis completed!")
#             self.logger.info(f"{'=' * 60}")
#             self.logger.info(f"输出文件位于 | Output files in: {self.config.output_dir}")
#             self.logger.info("  - all.fpkm.tpm.txt: 所有样本的FPKM和TPM矩阵 | FPKM and TPM matrix for all samples.")

#         except Exception as e:
#             self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
#             sys.exit(1)

# def main():
#     """主函数 | Main function"""
#     parser = argparse.ArgumentParser(
#         description="RNA-seq分析流程：HISAT2 + StringTie (模块化版本) | RNA-seq analysis pipeline: HISAT2 + StringTie (Modular Version)",
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter
#     )

#     # 必需参数 | Required parameters
#     parser.add_argument("-g", "--genome", required=True, 
#                        help="基因组fasta文件路径 | Genome fasta file path")
#     parser.add_argument("-f", "--gtf", required=True, 
#                        help="基因注释GTF文件路径 | Gene annotation GTF file path")
#     parser.add_argument("-i", "--input", required=True,
#                        help="输入fastq文件目录或样本信息文件 | Input fastq file directory or sample information file")
#     parser.add_argument("-o", "--output", required=True, 
#                        help="输出目录 | Output directory")

#     # 可选参数 | Optional parameters
#     parser.add_argument("-p", "--pattern", default=None,
#                        help='Fastq文件命名模式，例如 "*.R1.fastq.gz" 或 "*_1.fq.gz"，*代表样本名 | '
#                             'Fastq file naming pattern, e.g., "*.R1.fastq.gz" or "*_1.fq.gz", * represents sample name')
#     parser.add_argument("-r", "--remove", default="no", choices=["yes", "y", "no", "n"],
#                        help="处理后删除BAM文件 | Remove BAM files after processing")
#     parser.add_argument("-t", "--threads", type=int, default=8, 
#                        help="线程数 | Number of threads")

#     args = parser.parse_args()

#     # 创建分析器并运行 | Create analyzer and run
#     analyzer = RNASeqAnalyzer(
#         genome_file=args.genome,
#         gtf_file=args.gtf,
#         input_path=args.input,
#         output_dir=args.output,
#         threads=args.threads,
#         fastq_pattern=args.pattern,
#         remove_bam=args.remove
#     )
    
#     analyzer.run_analysis()

# if __name__ == "__main__":
#     main()

"""
RNA-seq分析主程序模块 | RNA-seq Analysis Main Module
"""

import argparse
import sys
import os
from .config import RNASeqConfig
from .utils import RNASeqLogger, CommandRunner
from .data_processing import SampleParser
from .alignment import HISAT2Indexer, HISAT2Aligner
from .quantification import StringTieQuantifier, GTFValueExtractor
from .results import ExpressionMatrixMerger, SummaryGenerator

class RNASeqAnalyzer:
    """RNA-seq分析主类 | Main RNA-seq Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = RNASeqConfig(**kwargs)
        self.config.validate()

        # 初始化日志 | Initialize logging
        self.logger_manager = RNASeqLogger(
            self.config.output_path,
            verbose=self.config.verbose,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化各个处理器 | Initialize processors
        self.sample_parser = SampleParser(self.logger)
        self.hisat2_indexer = HISAT2Indexer(self.config, self.logger, self.cmd_runner)
        self.hisat2_aligner = HISAT2Aligner(self.config, self.logger, self.cmd_runner)
        self.stringtie_quantifier = StringTieQuantifier(self.config, self.logger, self.cmd_runner)
        self.gtf_extractor = GTFValueExtractor(self.config, self.logger)
        self.matrix_merger = ExpressionMatrixMerger(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def process_single_sample(self, sample_info: dict, index_prefix: str) -> str:
        """处理单个样本的完整工作流程|Process complete workflow for a single sample"""
        sample_name = sample_info["name"]
        fastq1 = sample_info["fastq1"]
        fastq2 = sample_info["fastq2"]

        self.logger.info("=" * 60)
        self.logger.info(f"处理样本|Processing sample: {sample_name}")
        self.logger.info("=" * 60)

        # 设置文件路径|Set file paths
        bam_file = os.path.join(self.config.output_dir, f"{sample_name}.sorted.bam")
        stringtie_output = os.path.join(self.config.output_dir, "stringtie_output", f"{sample_name}.gtf")
        fpkm_output = os.path.join(self.config.output_dir, "fpkm_output", f"{sample_name}.fpkm.txt")

        # 1. HISAT2比对|HISAT2 alignment
        if not self.hisat2_aligner.run_hisat2_mapping(index_prefix, fastq1, fastq2, bam_file):
            self.logger.warning(f"样本比对失败或超时，跳过该样本|Sample alignment failed or timed out, skipping: {sample_name}")
            return None

        # 2. StringTie定量|StringTie quantification
        if not self.stringtie_quantifier.run_stringtie(bam_file, stringtie_output):
            self.logger.warning(f"样本定量失败，跳过该样本|Sample quantification failed, skipping: {sample_name}")
            return None

        # 3. 提取FPKM值|Extract FPKM values
        if not self.gtf_extractor.extract_gtf_values(stringtie_output, sample_name, fpkm_output):
            self.logger.warning(f"FPKM提取失败，跳过该样本|FPKM extraction failed, skipping: {sample_name}")
            return None

        # 4. 如果请求，删除BAM文件|Remove BAM file if requested
        if self.config.remove_bam.lower() in ["yes", "y"]:
            if os.path.exists(bam_file):
                os.remove(bam_file)
                self.logger.info(f"已删除BAM文件|Removed BAM file: {bam_file}")
        else:
            self.logger.info(f"保留BAM文件|BAM file retained: {bam_file}")

        return fpkm_output
    
    def run_analysis(self):
        """运行完整的RNA-seq分析流程|Run complete RNA-seq analysis pipeline"""
        import time
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("RNA-seq分析流程开始|RNA-seq analysis pipeline started")
            self.logger.info("=" * 60)
            self.logger.info(f"基因组文件|Genome file: {self.config.genome_file}")
            self.logger.info(f"输入路径|Input path: {self.config.input_path}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"GTF文件|GTF file: {self.config.gtf_file}")
            self.logger.info(f"线程数|Threads: {self.config.threads}")
            self.logger.info(f"删除BAM文件|Remove BAM files: {self.config.remove_bam}")
            if self.config.fastq_pattern:
                self.logger.info(f"文件模式|File pattern: {self.config.fastq_pattern}")

            # 步骤1: 构建HISAT2索引|Step 1: Build HISAT2 index
            self.logger.info("=" * 60)
            self.logger.info("步骤1: 构建HISAT2索引|Step 1: Building HISAT2 index")
            self.logger.info("=" * 60)
            index_prefix = self.hisat2_indexer.build_hisat2_index()

            # 步骤2: 解析输入样本|Step 2: Parse input samples
            self.logger.info("=" * 60)
            self.logger.info("步骤2: 解析输入样本|Step 2: Parsing input samples")
            self.logger.info("=" * 60)
            samples = self.sample_parser.parse_input_samples(self.config.input_path, self.config.fastq_pattern)

            if not samples:
                self.logger.error("错误：未找到有效的样本文件|Error: No valid sample files found")
                sys.exit(1)

            self.logger.info(f"找到 {len(samples)} 个样本|Found {len(samples)} samples:")
            for sample in samples:
                self.logger.info(f"  - {sample['name']}: {sample['fastq1']}, {sample['fastq2']}")

            # 保存样本信息到配置|Save sample info to config
            self.config.samples = samples

            # 步骤3: 处理所有样本|Step 3: Process all samples
            self.logger.info("=" * 60)
            self.logger.info("步骤3: 处理所有样本|Step 3: Processing all samples")
            self.logger.info("=" * 60)
            fpkm_files = []

            for sample_info in samples:
                fpkm_file = self.process_single_sample(sample_info, index_prefix)
                if fpkm_file:
                    fpkm_files.append(fpkm_file)
                else:
                    self.logger.error(f"样本处理失败|Sample processing failed: {sample_info['name']}")

            # 步骤4: 合并表达矩阵|Step 4: Merge expression matrix
            self.logger.info("=" * 60)
            self.logger.info("步骤4: 合并表达矩阵|Step 4: Merging expression matrix")
            self.logger.info("=" * 60)

            if not self.matrix_merger.merge_expression_matrix(fpkm_files):
                self.logger.error("表达矩阵合并失败|Expression matrix merging failed")
                sys.exit(1)

            # 生成总结报告|Generate summary report
            self.summary_generator.generate_summary_report()

            # 输出总结信息|Output summary
            elapsed_time = time.time() - start_time
            self.logger.info("=" * 60)
            self.logger.info("分析流程总结|Pipeline Summary")
            self.logger.info("=" * 60)
            self.logger.info(f"总运行时间|Total runtime: {elapsed_time:.2f} 秒|seconds")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info("分析流程成功完成|Pipeline completed successfully")
            self.logger.info(f"主要输出文件|Main output file: all.fpkm.tpm.txt")
            self.logger.info("=" * 60)

        except KeyboardInterrupt:
            self.logger.warning("操作被用户中断|Operation interrupted by user")
            sys.exit(130)
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止|Analysis pipeline terminated unexpectedly: {e}", exc_info=True)
            sys.exit(1)

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="RNA-seq分析流程：HISAT2 + StringTie (模块化版本)|RNA-seq analysis pipeline: HISAT2 + StringTie (Modular Version)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument("-g", "--genome", required=True,
                       help="基因组fasta文件路径|Genome fasta file path")
    required.add_argument("-f", "--gtf", required=True,
                       help="基因注释GTF文件路径|Gene annotation GTF file path")
    required.add_argument("-i", "--input", required=True,
                       help="输入fastq文件目录或样本信息文件|Input fastq file directory or sample information file")
    required.add_argument("-o", "--output", required=True,
                       help="输出目录|Output directory")

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('可选参数|Optional parameters')
    optional.add_argument("-p", "--pattern", default=None,
                       help='Fastq文件命名模式，例如 "*.R1.fastq.gz" 或 "*_1.fq.gz"，*代表样本名|'
                            'Fastq file naming pattern, e.g., "*.R1.fastq.gz" or "*_1.fq.gz", * represents sample name')
    optional.add_argument("-r", "--remove", default="no", choices=["yes", "y", "no", "n"],
                       help="处理后删除BAM文件|Remove BAM files after processing")
    optional.add_argument("-t", "--threads", type=int, default=8,
                       help="线程数|Number of threads")
    optional.add_argument("--sample-timeout", type=int, default=21600,
                       help="单个样本处理超时时间（秒），默认6小时|Sample processing timeout in seconds (default: 21600 = 6 hours)")

    # 日志选项|Logging options
    logging = parser.add_argument_group('日志选项|Logging options')
    logging.add_argument('-v', '--verbose',
                        action='count',
                        default=0,
                        help='增加输出详细程度 (-v, -vv, -vvv)|Increase output verbosity')
    logging.add_argument('--quiet',
                        action='store_true',
                        help='静默模式，仅输出错误信息|Quiet mode, only output errors')
    logging.add_argument('--log-file',
                        type=str,
                        help='日志文件路径|Log file path')
    logging.add_argument('--log-level',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO',
                        help='日志级别|Log level')

    # 高级选项|Advanced options
    advanced = parser.add_argument_group('高级选项|Advanced options')
    advanced.add_argument('--dry-run',
                         action='store_true',
                         help='试运行模式，不实际执行|Dry run mode, no actual execution')

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    try:
        analyzer = RNASeqAnalyzer(
            genome_file=args.genome,
            gtf_file=args.gtf,
            input_path=args.input,
            output_dir=args.output,
            threads=args.threads,
            fastq_pattern=args.pattern,
            remove_bam=args.remove,
            sample_timeout=args.sample_timeout,
            verbose=(args.verbose > 0),
            quiet=args.quiet,
            log_file=args.log_file,
            log_level=args.log_level,
            dry_run=args.dry_run
        )

        analyzer.run_analysis()

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n操作被用户中断|Operation interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()