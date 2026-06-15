"""
启动子提取器命令行接口模块|Promoter Extractor Command Line Interface Module
"""

import argparse
import sys
import os
import time
from pathlib import Path
from .config import PromoterExtractorConfig
from .utils import PromoterLogger
from .calculator import PromoterExtractor

# 版本信息|Version information
VERSION = "1.0.0"


class PromoterRunner:
    """启动子提取器运行器|Promoter Extractor Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = PromoterExtractorConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = PromoterLogger(
            self.config.output_path,
            log_name="promoter_extractor.log",
            log_level=self.config.log_level,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化提取器|Initialize extractor
        self.extractor = PromoterExtractor(self.config, self.logger)

    def run(self):
        """运行启动子提取流程|Run promoter extraction pipeline"""
        try:
            self.logger.info("=" * 80)
            self.logger.info("启动启动子提取器|Starting Promoter Extractor")
            self.logger.info("=" * 80)
            self.logger.info(f"输入GFF文件|Input GFF file: {self.config.gff_file}")
            self.logger.info(f"输入基因组文件|Input genome file: {self.config.genome_file}")
            self.logger.info(f"输出前缀|Output prefix: {self.config.output_prefix}")
            self.logger.info(f"启动子长度|Promoter length: {self.config.promoter_length}bp")
            self.logger.info("=" * 80)

            # 提取启动子|Extract promoters
            fasta_records, bed_records = self.extractor.extract_promoters()

            # 写入输出文件|Write output files
            self._write_outputs(fasta_records, bed_records)

            # 完成信息|Completion information
            self.logger.info("=" * 80)
            self.logger.info("启动子提取完成|Promoter extraction completed successfully!")
            self.logger.info("=" * 80)
            self.logger.info(f"提取的启动子数量|Promoters extracted: {self.extractor.calculator.stats['extracted']}")
            self.logger.info(f"总基因数|Total genes: {self.extractor.calculator.stats['total_genes']}")
            self.logger.info(f"边界截断|Boundary truncated: {self.extractor.calculator.stats['skipped_boundary']}")
            self.logger.info(f"长度不足|Too short: {self.extractor.calculator.stats['skipped_length']}")
            self.logger.info("=" * 80)
            self.logger.info("输出文件|Output files:")

            fasta_file = f"{self.config.output_prefix}.fa"
            self.logger.info(f"  - {fasta_file}: 启动子序列FASTA|Promoter sequences FASTA")

            if self.config.output_bed:
                bed_file = f"{self.config.output_prefix}.bed"
                self.logger.info(f"  - {bed_file}: 启动子位置BED|Promoter positions BED")

            if self.config.output_stats:
                stats_file = f"{self.config.output_prefix}_stats.txt"
                self.logger.info(f"  - {stats_file}: 统计信息|Statistics")

        except Exception as e:
            self.logger.error(f"提取过程中发生严重错误|A critical error occurred during extraction: {e}")
            sys.exit(1)

    def _write_outputs(self, fasta_records, bed_records):
        """
        写入输出文件|Write output files

        Args:
            fasta_records: FASTA记录列表|FASTA records list
            bed_records: BED记录列表|BED records list
        """
        # 写入FASTA文件|Write FASTA file
        fasta_file = f"{self.config.output_prefix}.fa"
        self._write_fasta(fasta_file, fasta_records)

        # 写入BED文件|Write BED file
        if self.config.output_bed:
            bed_file = f"{self.config.output_prefix}.bed"
            self._write_bed(bed_file, bed_records)

        # 写入统计文件|Write statistics file
        if self.config.output_stats:
            stats_file = f"{self.config.output_prefix}_stats.txt"
            self._write_stats(stats_file)

    def _write_fasta(self, filename: str, records: list):
        """
        写入FASTA文件|Write FASTA file

        Args:
            filename: 文件名|Filename
            records: 记录列表|Records list [(header, sequence), ...]
        """
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                for header, sequence in records:
                    f.write(f">{header}\n")
                    # 每行80个字符|80 characters per line
                    for i in range(0, len(sequence), 80):
                        f.write(sequence[i:i+80] + "\n")

            self.logger.info(f"FASTA文件已写入|FASTA file written: {filename}")

        except Exception as e:
            self.logger.error(f"写入FASTA文件失败|Failed to write FASTA file: {e}")
            raise

    def _write_bed(self, filename: str, records: list):
        """
        写入BED文件|Write BED file

        Args:
            filename: 文件名|Filename
            records: 记录列表|Records list [(seqid, start, end, name, score, strand), ...]
        """
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                for seqid, start, end, name, score, strand in records:
                    # BED格式：seqid start end name score strand
                    f.write(f"{seqid}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

            self.logger.info(f"BED文件已写入|BED file written: {filename}")

        except Exception as e:
            self.logger.error(f"写入BED文件失败|Failed to write BED file: {e}")
            raise

    def _write_stats(self, filename: str):
        """
        写入统计文件|Write statistics file

        Args:
            filename: 文件名|Filename
        """
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                stats = self.extractor.calculator.stats
                f.write("启动子提取统计|Promoter Extraction Statistics\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"总基因数|Total genes: {stats['total_genes']}\n")
                f.write(f"成功提取|Successfully extracted: {stats['extracted']}\n")
                f.write(f"边界截断|Boundary truncated: {stats['skipped_boundary']}\n")
                f.write(f"长度不足|Too short (< {self.config.min_length}bp): {stats['skipped_length']}\n")
                f.write(f"成功率|Success rate: ")
                if stats['total_genes'] > 0:
                    rate = stats['extracted'] / stats['total_genes'] * 100
                    f.write(f"{rate:.2f}%\n")
                else:
                    f.write("N/A\n")
                f.write(f"\n启动子长度|Promoter length: {self.config.promoter_length}bp\n")
                f.write(f"最小长度要求|Min length requirement: {self.config.min_length}bp\n")

            self.logger.info(f"统计文件已写入|Statistics file written: {filename}")

        except Exception as e:
            self.logger.error(f"写入统计文件失败|Failed to write statistics file: {e}")
            raise


def main():
    """主函数入口|Main function entry point"""
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="启动子提取器工具|Promoter Extractor Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i genes.gff3 -g genome.fa -o promoters
  %(prog)s -i genes.gff3 -g genome.fa -o promoters -p 1500 --gene-list gene_ids.txt
  %(prog)s -i genes.gff3 -g genome.fa -o promoters --no-bed
        '''
    )

    # 必需参数|Required arguments
    parser.add_argument("-i", "--gff", required=True,
                       help="输入GFF3文件路径|Input GFF3 file path")
    parser.add_argument("-g", "--genome", required=True,
                       help="输入基因组FASTA文件路径|Input genome FASTA file path")
    parser.add_argument("-o", "--output", default="promoters",
                       help="输出前缀|Output prefix (default: promoters)")

    # 启动子参数|Promoter parameters
    parser.add_argument("-p", "--promoter-length", type=int, default=2000,
                       help="启动子长度（bp）|Promoter length in bp (default: 2000)")
    parser.add_argument("--min-length", type=int, default=0,
                       help="最小接受长度（bp）|Minimum acceptable length in bp (default: 0)")

    # 基因选择|Gene selection
    parser.add_argument("--gene-list",
                       help="基因ID列表文件|Gene ID list file (one gene ID per line)")

    # 输出选项|Output options
    parser.add_argument("--no-bed", action="store_true",
                       help="不输出BED格式文件|Do not output BED format file")
    parser.add_argument("--no-stats", action="store_true",
                       help="不输出统计文件|Do not output statistics file")

    # 处理选项|Processing options
    parser.add_argument("-t", "--threads", type=int, default=1,
                       help="线程数|Number of threads (default: 1)")
    parser.add_argument("--keep-intermediate", action="store_true",
                       help="保留中间文件|Keep intermediate files")

    # 日志参数|Logging parameters
    log_group = parser.add_argument_group('日志选项|Logging options')
    log_group.add_argument("--verbose", action="count", default=0,
                          help="详细输出模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)")
    log_group.add_argument("--quiet", action="store_true",
                          help="静默模式(只输出ERROR)|Quiet mode (ERROR only)")
    log_group.add_argument("--log-level",
                          help="日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)|Log level (default: INFO)")

    # 执行控制|Execution control
    exec_group = parser.add_argument_group('执行选项|Execution options')
    exec_group.add_argument("-f", "--force", action="store_true",
                           help="强制覆盖已存在文件|Force overwrite existing files")
    exec_group.add_argument("--dry-run", action="store_true",
                           help="模拟运行(不实际执行)|Dry run without execution")

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

    # 创建运行器并运行|Create runner and run
    runner = PromoterRunner(
        gff_file=args.gff,
        genome_file=args.genome,
        output_prefix=args.output,
        promoter_length=args.promoter_length,
        min_length=args.min_length,
        gene_list=args.gene_list,
        output_bed=not args.no_bed,
        output_stats=not args.no_stats,
        threads=args.threads,
        keep_intermediate=args.keep_intermediate,
        log_level=log_level,
        quiet=args.quiet,
        verbose=args.verbose,
        force=args.force,
        dry_run=args.dry_run
    )

    try:
        if args.dry_run:
            runner.logger.info("模拟运行模式-不会实际执行命令|DRY RUN mode - commands will not be executed")

        # 执行提取|Run extraction
        runner.run()

        # 输出总结信息|Output summary
        elapsed_time = time.time() - start_time
        runner.logger.info("=" * 80)
        runner.logger.info("流程总结|Pipeline Summary")
        runner.logger.info("=" * 80)
        runner.logger.info(f"总运行时间|Total runtime: {elapsed_time:.2f} seconds")
        runner.logger.info(f"输出前缀|Output prefix: {args.output}")
        runner.logger.info("流程完成|Pipeline completed successfully")

    except KeyboardInterrupt:
        runner.logger.warning("用户中断程序执行|Process interrupted by user")
        sys.exit(130)
    except Exception as e:
        runner.logger.critical(f"Pipeline failed: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
