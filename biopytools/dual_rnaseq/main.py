"""
互作转录组分析主程序模块|Dual RNA-seq Analysis Main Module
"""

import argparse
import sys
import os
from .config import DualRNASeqConfig
from .utils import DualRNASeqLogger, CommandRunner
from .data_processing import SampleParser
from .indexing import DualIndexBuilder
from .classification import ReadsClassifier
from .quantification import DualQuantifier
from .results import DualMatrixMerger, SummaryGenerator
from .bam_to_fastq import BamToFastqExtractor
from .mapping_stats import MappingStatistics


class DualRNASeqAnalyzer:
    """互作转录组分析主类|Main Dual RNA-seq Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = DualRNASeqConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = DualRNASeqLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化各个处理器|Initialize processors
        self.sample_parser = SampleParser(self.logger)
        self.index_builder = DualIndexBuilder(self.config, self.logger, self.cmd_runner)
        self.reads_classifier = ReadsClassifier(self.config, self.logger, self.cmd_runner)
        self.dual_quantifier = DualQuantifier(self.config, self.logger, self.cmd_runner)
        self.matrix_merger = DualMatrixMerger(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
        self.bam_to_fastq_extractor = BamToFastqExtractor(self.config, self.logger, self.cmd_runner)
        self.mapping_stats = MappingStatistics(self.config, self.logger)

    def run_analysis(self):
        """运行完整的互作转录组分析流程|Run complete dual RNA-seq analysis pipeline"""
        try:
            self.logger.info("互作转录组分析流程开始|Dual RNA-seq analysis pipeline started")
            self.logger.info(f"物种1|Species 1 ({self.config.species1_name}):")
            self.logger.info(f"  基因组|Genome: {self.config.species1_genome}")
            self.logger.info(f"  GTF: {self.config.species1_gtf}")
            self.logger.info(f"物种2|Species 2 ({self.config.species2_name}):")
            self.logger.info(f"  基因组|Genome: {self.config.species2_genome}")
            self.logger.info(f"  GTF: {self.config.species2_gtf}")
            self.logger.info(f"输入路径|Input path: {self.config.input_path}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"线程数|Threads: {self.config.threads}")
            self.logger.info(f"最小MAPQ|Minimum MAPQ: {self.config.min_mapq}")

            # 步骤1: 构建双物种HISAT2索引|Step 1: Build dual-species HISAT2 indexes
            self.logger.info("=" * 60)
            self.logger.info("步骤1: 构建双物种HISAT2索引|Step 1: Building dual-species HISAT2 indexes")
            self.logger.info("=" * 60)

            species1_index, species2_index = self.index_builder.build_dual_indexes()

            if not species1_index or not species2_index:
                self.logger.error("索引构建失败|Index building failed")
                sys.exit(1)

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

            # 步骤3: Reads物种分类|Step 3: Reads species classification
            self.logger.info("=" * 60)
            self.logger.info("步骤3: Reads物种分类|Step 3: Reads species classification")
            self.logger.info("=" * 60)

            for sample in samples:
                sample_name = sample["name"]
                fastq1 = sample["fastq1"]
                fastq2 = sample["fastq2"]

                if not self.reads_classifier.process_sample_classification(
                    fastq1, fastq2, species1_index, species2_index, sample_name
                ):
                    self.logger.error(f"样本分类失败|Sample classification failed: {sample_name}")
                    sys.exit(1)

            # 步骤3: 比对统计|Step 3: Alignment statistics
            self.logger.info("=" * 60)
            self.logger.info("步骤3: 计算比对统计|Step 3: Calculating alignment statistics")
            self.logger.info("=" * 60)

            # 收集所有样本的BAM文件路径|Collect BAM file paths for all samples
            bam_files = []
            temp_dir = os.path.join(self.config.output_dir, "02.classification", ".temp")

            for sample in samples:
                sample_name = sample["name"]
                species1_bam = os.path.join(temp_dir, sample_name, f"{sample_name}.{self.config.species1_name}.raw.bam")
                species2_bam = os.path.join(temp_dir, sample_name, f"{sample_name}.{self.config.species2_name}.raw.bam")
                bam_files.append((sample_name, species1_bam, species2_bam))

            # 生成比对统计报告|Generate alignment statistics report
            stats_dir = os.path.join(self.config.output_dir, "03.alignment_statistics")
            os.makedirs(stats_dir, exist_ok=True)

            if not self.mapping_stats.process_all_samples(bam_files, stats_dir):
                self.logger.warning("比对统计生成失败|Alignment statistics generation failed")
            else:
                self.logger.info("比对统计生成成功|Alignment statistics generated successfully")

            # 步骤4: 双物种定量分析|Step 4: Dual-species quantification
            self.logger.info("=" * 60)
            self.logger.info("步骤4: 双物种定量分析|Step 4: Dual-species quantification")
            self.logger.info("=" * 60)

            if not self.dual_quantifier.quantify_all_samples(samples):
                self.logger.error("定量分析失败|Quantification analysis failed")
                sys.exit(1)

            # 步骤5: 生成双物种表达矩阵|Step 5: Generate dual-species expression matrix
            self.logger.info("=" * 60)
            self.logger.info("步骤5: 生成双物种表达矩阵|Step 5: Generating dual-species expression matrix")
            self.logger.info("=" * 60)

            # 获取定量结果文件|Get quantification result files
            quant_dir = os.path.join(self.config.output_dir, "04.quantification")
            species1_fpkm_files = []
            species2_fpkm_files = []

            for sample in samples:
                sample_name = sample["name"]
                species1_fpkm = os.path.join(quant_dir, self.config.species1_name, f"{sample_name}.fpkm.txt")
                species2_fpkm = os.path.join(quant_dir, self.config.species2_name, f"{sample_name}.fpkm.txt")

                if os.path.exists(species1_fpkm):
                    species1_fpkm_files.append(species1_fpkm)
                if os.path.exists(species2_fpkm):
                    species2_fpkm_files.append(species2_fpkm)

            # 合并表达矩阵|Merge expression matrices
            if species1_fpkm_files:
                self.matrix_merger.merge_expression_matrix(species1_fpkm_files, self.config.species1_name)
            if species2_fpkm_files:
                self.matrix_merger.merge_expression_matrix(species2_fpkm_files, self.config.species2_name)

            # 生成总结报告|Generate summary report
            self.summary_generator.generate_summary_report()

            # 步骤6: 从BAM文件提取FASTQ（可选）|Step 6: Extract FASTQ from BAM (optional)
            if self.config.extract_fastq:
                self.logger.info("=" * 60)
                self.logger.info("步骤6: 从BAM文件提取FASTQ|Step 6: Extracting FASTQ from BAM files")
                self.logger.info("=" * 60)

                sample_names = [sample["name"] for sample in samples]
                if not self.bam_to_fastq_extractor.extract_all_samples(sample_names):
                    self.logger.warning("FASTQ提取部分失败|FASTQ extraction partially failed")
                else:
                    self.logger.info("FASTQ提取完成|FASTQ extraction completed")

            self.logger.info("=" * 60)
            self.logger.info("分析完成！|Analysis completed!")
            self.logger.info("=" * 60)
            self.logger.info(f"输出文件位于|Output files in: {self.config.output_dir}")
            self.logger.info(f"  - {self.config.species1_name}_matrix.txt: {self.config.species1_name}表达矩阵|expression matrix")
            self.logger.info(f"  - {self.config.species2_name}_matrix.txt: {self.config.species2_name}表达矩阵|expression matrix")
            if self.config.extract_fastq:
                self.logger.info(f"  - 06.extracted_fastq/: 提取的FASTQ文件|Extracted FASTQ files")
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
        description="互作转录组分析流程：双物种转录组同时分析|Dual RNA-seq analysis pipeline: Simultaneous analysis of dual-species transcriptome",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s --species1-name host --species1-genome host.fa --species1-gtf host.gtf \\
           --species2-name pathogen --species2-genome pathogen.fa --species2-gtf pathogen.gtf \\
           -i ./fastq_data -o results
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')

    required.add_argument("--species1-name", required=True,
                         help="物种1名称|Species 1 name (e.g., host)")
    required.add_argument("--species1-genome", required=True,
                         help="物种1基因组FASTA文件|Species 1 genome FASTA file")
    required.add_argument("--species1-gtf", required=True,
                         help="物种1GTF注释文件|Species 1 GTF annotation file")

    required.add_argument("--species2-name", required=True,
                         help="物种2名称|Species 2 name (e.g., pathogen)")
    required.add_argument("--species2-genome", required=True,
                         help="物种2基因组FASTA文件|Species 2 genome FASTA file")
    required.add_argument("--species2-gtf", required=True,
                         help="物种2GTF注释文件|Species 2 GTF annotation file")

    required.add_argument("-i", "--input", required=True,
                         help="输入FASTQ文件目录或样本信息文件|Input FASTQ file directory or sample information file")
    required.add_argument("-o", "--output-dir", required=True,
                         help="输出目录|Output directory")

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('可选参数|Optional parameters')

    optional.add_argument("-p", "--pattern", default="*_1.clean.fq.gz",
                         help='FASTQ文件命名模式|FASTQ file naming pattern (e.g., "*.R1.fastq.gz")')
    optional.add_argument("-t", "--threads", type=int, default=12,
                         help="线程数|Number of threads")

    # 分类参数|Classification parameters
    classification = parser.add_argument_group('分类参数|Classification parameters')

    classification.add_argument("--min-mapq", type=int, default=20,
                               help="最小mapping quality值|Minimum mapping quality value")
    classification.add_argument("--no-unique-only", dest="unique_only", action="store_false",
                               help="不禁用非唯一比对|Do not disable non-unique mappings")
    classification.add_argument("--no-extract-fastq", dest="extract_fastq", action="store_false",
                               help="不提取FASTQ文件|Do not extract FASTQ files from BAM (default: extract)")

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    try:
        analyzer = DualRNASeqAnalyzer(
            species1_name=args.species1_name,
            species1_genome=args.species1_genome,
            species1_gtf=args.species1_gtf,
            species2_name=args.species2_name,
            species2_genome=args.species2_genome,
            species2_gtf=args.species2_gtf,
            input_path=args.input,
            output_dir=args.output_dir,
            threads=args.threads,
            fastq_pattern=args.pattern,
            min_mapq=args.min_mapq,
            unique_only=args.unique_only,
            extract_fastq=args.extract_fastq
        )

        analyzer.run_analysis()

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("操作被用户中断|Operation interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
