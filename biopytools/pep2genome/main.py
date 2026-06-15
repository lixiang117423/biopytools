"""
蛋白质到基因组比对主程序模块|Protein to Genome Alignment Main Module
"""

import argparse
import sys
import os
from .config import Pep2GenomeConfig
from .utils import Pep2GenomeLogger, CommandRunner
from .calculator import MiniprotAligner
from .parser import PAFParser
from .statistics import StatisticsGenerator
from .gff3_exporter import GFF3Exporter, BEDExporter
from .sequence_extractor import SequenceExtractor


class Pep2GenomeAnalyzer:
    """蛋白质到基因组比对分析主类|Main Protein to Genome Alignment Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = Pep2GenomeConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = Pep2GenomeLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化各个处理器|Initialize processors
        self.aligner = MiniprotAligner(self.config, self.logger, self.cmd_runner)
        self.parser = PAFParser(self.logger)
        self.stats_generator = StatisticsGenerator(self.config, self.logger)
        self.gff3_exporter = GFF3Exporter(self.config, self.logger)
        self.bed_exporter = BEDExporter(self.config, self.logger)
        self.sequence_extractor = SequenceExtractor(self.config, self.logger, self.cmd_runner)

    def run_analysis(self):
        """运行完整的蛋白质到基因组比对分析流程|Run complete protein to genome alignment pipeline"""
        try:
            self.logger.info("蛋白质到基因组比对分析流程开始|Protein to genome alignment analysis started")
            self.logger.info(f"基因组文件|Genome file: {self.config.genome_fa}")
            self.logger.info(f"蛋白质文件|Protein file: {self.config.protein_fa}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"线程数|Threads: {self.config.threads}")

            # 步骤1: 运行Miniprot对齐|Step 1: Run Miniprot alignment
            self.logger.info("=" * 60)
            self.logger.info("步骤1: 运行Miniprot对齐|Step 1: Running Miniprot alignment")
            self.logger.info("=" * 60)

            paf_file = self.aligner.run_alignment()

            # 步骤2: 解析PAF文件|Step 2: Parse PAF file
            self.logger.info("=" * 60)
            self.logger.info("步骤2: 解析PAF文件|Step 2: Parsing PAF file")
            self.logger.info("=" * 60)

            records = self.parser.parse_file(paf_file)

            if not records:
                self.logger.warning("警告：未解析到任何比对记录|Warning: No alignment records parsed")
                sys.exit(0)

            # 步骤3: 生成统计报告|Step 3: Generate statistics report
            if self.config.export_statistics:
                self.logger.info("=" * 60)
                self.logger.info("步骤3: 生成统计报告|Step 3: Generating statistics report")
                self.logger.info("=" * 60)

                self.stats_generator.generate_statistics(records)

            # 步骤4: 导出GFF3格式|Step 4: Export GFF3 format
            if self.config.export_gff3:
                self.logger.info("=" * 60)
                self.logger.info("步骤4: 导出GFF3格式|Step 4: Exporting GFF3 format")
                self.logger.info("=" * 60)

                self.gff3_exporter.export_gff3(records)

            # 步骤5: 导出BED格式|Step 5: Export BED format
            bed_file = None
            if self.config.export_bed:
                self.logger.info("=" * 60)
                self.logger.info("步骤5: 导出BED格式|Step 5: Exporting BED format")
                self.logger.info("=" * 60)

                bed_file = self.bed_exporter.export_bed(records)

            # 步骤6: 提取基因组序列|Step 6: Extract genome sequences
            if self.config.extract_sequences and bed_file:
                self.logger.info("=" * 60)
                self.logger.info("步骤6: 提取基因组序列|Step 6: Extracting genome sequences")
                self.logger.info("=" * 60)

                self.sequence_extractor.extract_sequences(bed_file)

            # 完成|Completion
            self.logger.info("=" * 60)
            self.logger.info("分析完成！|Analysis completed!")
            self.logger.info("=" * 60)
            self.logger.info(f"输出文件位于|Output files in: {self.config.output_dir}")
            self.logger.info(f"  - alignment.paf: Miniprot比对结果|Miniprot alignment results")
            if self.config.export_statistics:
                self.logger.info(f"  - alignment_statistics.txt: 统计报告|Statistics report")
            if self.config.export_gff3:
                self.logger.info(f"  - alignment.gff3: GFF3格式注释|GFF3 format annotations")
            if self.config.export_bed:
                self.logger.info(f"  - alignment.bed: BED格式注释|BED format annotations")
            if self.config.extract_sequences:
                self.logger.info(f"  - alignment.fa: 提取的基因组序列|Extracted genome sequences")
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
        description="蛋白质到基因组比对分析流程：使用Miniprot进行蛋白质与基因组的比对|Protein to genome alignment pipeline: Align proteins to genome using Miniprot",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s --genome genome.fa --protein protein.fa -o results
  %(prog)s --genome genome.fa --protein protein.fa -o results -t 24 --no-gff3
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')

    required.add_argument("--genome", required=True,
                         help="基因组FASTA文件|Genome FASTA file")
    required.add_argument("--protein", required=True,
                         help="蛋白质FASTA文件|Protein FASTA file")
    required.add_argument("-o", "--output-dir", required=True,
                         help="输出目录|Output directory")

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('可选参数|Optional parameters')

    optional.add_argument("-t", "--threads", type=int, default=12,
                         help="线程数|Number of threads")
    optional.add_argument("--miniprot-path", default="miniprot",
                         help="Miniprot工具路径|Miniprot tool path (default: miniprot)")

    # 输出选项|Output options
    output = parser.add_argument_group('输出选项|Output options')

    output.add_argument("--no-gff3", dest="export_gff3", action="store_false",
                       help="不导出GFF3格式|Do not export GFF3 format")
    output.add_argument("--no-bed", dest="export_bed", action="store_false",
                       help="不导出BED格式|Do not export BED format")
    output.add_argument("--no-statistics", dest="export_statistics", action="store_false",
                       help="不生成统计报告|Do not generate statistics report")
    output.add_argument("--no-extract", dest="extract_sequences", action="store_false",
                       help="不提取基因组序列|Do not extract genome sequences")

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    try:
        analyzer = Pep2GenomeAnalyzer(
            genome_fa=args.genome,
            protein_fa=args.protein,
            output_dir=args.output_dir,
            threads=args.threads,
            miniprot_path=args.miniprot_path,
            export_gff3=args.export_gff3,
            export_bed=args.export_bed,
            export_statistics=args.export_statistics,
            extract_sequences=args.extract_sequences
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
