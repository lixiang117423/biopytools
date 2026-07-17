"""
最长转录本提取主程序模块|Longest mRNA Extraction Main Module
"""

import argparse
import sys
from .config import LongestMRNAConfig
from .utils import LongestMRNALogger, CommandRunner, TempFileManager
from .data_processing import CDSCalculator, GFFGenomeAligner
from .extraction import SequenceExtractor
from .results import GeneInfoGenerator, StatisticsCalculator, SummaryGenerator

class LongestMRNAExtractor:
    """最长转录本提取主类|Main Longest mRNA Extractor Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = LongestMRNAConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = LongestMRNALogger()
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化各个处理器|Initialize processors
        self.cds_calculator = CDSCalculator(self.logger)
        self.sequence_extractor = SequenceExtractor(self.config, self.logger, self.cmd_runner)
        self.gene_info_generator = GeneInfoGenerator(self.config, self.logger)
        self.stats_calculator = StatisticsCalculator(self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)

    def run_extraction(self):
        """运行提取流程|Run extraction pipeline"""
        # 管理对齐阶段产生的过滤GFF临时文件|Manages filtered GFF temp file produced by the alignment step
        gff_temp_manager = TempFileManager(self.logger)
        try:
            self.logger.info("开始最长转录本提取流程|Starting longest mRNA extraction pipeline")
            self.logger.info(f"基因组文件|Genome file: {self.config.genome_file}")
            self.logger.info(f"GFF3文件|GFF3 file: {self.config.gff3_file}")
            self.logger.info(f"输出文件|Output file: {self.config.output_file}")

            # 步骤0: 对齐GFF与基因组序列名，跳过基因组中不存在的序列|Step 0: Align GFF seqids with genome, skip seqids absent from genome
            self.logger.info("步骤0: 对齐GFF与基因组序列名|Step 0: Aligning GFF seqids with genome FASTA")
            aligner = GFFGenomeAligner(self.config.genome_file, self.config.gff3_file, self.logger)
            effective_gff = aligner.align(gff_temp_manager)
            # 重定向后续所有步骤使用(可能的)过滤后GFF，保证最长转录本/基因信息/序列输出三者一致|
            # Redirect all downstream steps to the (possibly filtered) GFF so that longest-transcript,
            # gene-info and sequence outputs stay consistent with each other
            self.config.gff3_file = effective_gff

            # 步骤1: 计算最长转录本|Step 1: Calculate longest transcripts
            self.logger.info("步骤1: 分析GFF3文件，计算最长转录本|Step 1: Analyzing GFF3 file, calculating longest transcripts")
            longest_transcripts = self.cds_calculator.calculate_from_gff(self.config.gff3_file)

            if not longest_transcripts:
                self.logger.error("未找到任何转录本信息|No transcript information found")
                sys.exit(1)

            # 步骤2: 生成基因信息文件|Step 2: Generate gene info file
            self.logger.info("步骤2: 生成基因信息文件|Step 2: Generating gene info file")
            self.gene_info_generator.generate_gene_info(longest_transcripts, self.cds_calculator.gene_metadata)

            # 步骤3: 提取蛋白质序列|Step 3: Extract protein sequences
            self.logger.info("步骤3: 提取蛋白质序列|Step 3: Extracting protein sequences")
            if not self.sequence_extractor.extract_protein_sequences(longest_transcripts):
                self.logger.error("序列提取失败|Sequence extraction failed")
                sys.exit(1)

            # 步骤4: 提取CDS序列(可选)|Step 4: Extract CDS sequences (optional)
            if self.config.cds_output_file:
                self.logger.info("步骤4: 提取CDS核苷酸序列|Step 4: Extracting CDS nucleotide sequences")
                if not self.sequence_extractor.extract_cds_sequences(longest_transcripts):
                    self.logger.error("CDS序列提取失败|CDS sequence extraction failed")
                    sys.exit(1)

            # 步骤5: 汇总统计|Step 5: Summary statistics
            self.logger.info("步骤5: 生成汇总统计|Step 5: Generating summary")
            stats = self.stats_calculator.calculate_statistics(
                self.cds_calculator.gene_transcripts, longest_transcripts
            )
            stats['noncoding_skipped'] = getattr(self.cds_calculator, 'skipped_noncoding', 0)
            self.summary_generator.print_summary(stats)

            self.logger.info("最长转录本提取流程完成|Longest mRNA extraction pipeline completed successfully")

        except Exception as e:
            self.logger.error(f"提取流程在执行过程中意外终止|Extraction pipeline terminated unexpectedly: {e}")
            sys.exit(1)
        finally:
            gff_temp_manager.cleanup()

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="最长转录本提取工具(模块化版本)|Longest mRNA Extraction Tool (Modular Version)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-g', '--genome', required=True,
                       help='输入基因组FASTA文件|Input genome FASTA file')
    parser.add_argument('-f', '--gff3', required=True,
                       help='输入GFF3注释文件|Input GFF3 annotation file')
    parser.add_argument('-o', '--output', required=True,
                       help='输出FASTA文件|Output FASTA file')

    # 可选参数|Optional arguments
    parser.add_argument('--gene-info',
                       help='基因信息输出文件|Gene info output file (auto-generated by default)')
    parser.add_argument('--cds-output',
                       help='CDS核苷酸序列输出文件|CDS nucleotide sequence output file')

    args = parser.parse_args()

    # 创建提取器并运行|Create extractor and run
    extractor = LongestMRNAExtractor(
        genome_file=args.genome,
        gff3_file=args.gff3,
        output_file=args.output,
        gene_info_file=args.gene_info,
        cds_output_file=args.cds_output
    )

    extractor.run_extraction()

if __name__ == "__main__":
    main()
