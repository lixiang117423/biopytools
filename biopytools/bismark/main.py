"""
Bismark流程主程序模块|Bismark Pipeline Main Module
"""

import argparse
import sys
from .config import BismarkConfig
from .utils import BismarkLogger, check_dependencies
from .core_analysis import CorePipeline
from . import __version__
from pathlib import Path
import os

class BismarkAnalyzer:
    def __init__(self, **kwargs):
        self.config = BismarkConfig(**kwargs)
        self.config.validate()
        
        os.makedirs(self.config.output_dir, exist_ok=True)
        os.makedirs(self.config.mapping_dir, exist_ok=True)
        os.makedirs(self.config.result_dir, exist_ok=True)
        os.makedirs(self.config.tmp_dir, exist_ok=True)

        self.logger_manager = BismarkLogger(Path(self.config.output_dir))
        self.logger = self.logger_manager.get_logger()
        self.core_pipeline = CorePipeline(self.config, self.logger)
    
    def run_analysis(self):
        try:
            self.logger.info(f"开始Bismark甲基化分析流程|Starting Bismark methylation analysis pipeline (v{__version__})")
            self.logger.info(f"主输出目录|Main output directory: {self.config.output_dir}")
            check_dependencies(self.config, self.logger)

            if not self.core_pipeline.build_bismark_index():
                raise RuntimeError("Bismark索引构建失败|Bismark index build failed")

            if not self.core_pipeline.run_mapping_and_extraction():
                raise RuntimeError("甲基化比对和提取失败|Methylation mapping and extraction failed")

            self.core_pipeline.results_manager.generate_summary_report()
            self.logger.info("=== Bismark分析流程全部完成|Bismark analysis pipeline completed successfully ===")

        except Exception as e:
            self.logger.error(f"分析流程意外终止|Analysis pipeline terminated unexpectedly: {e}")
            raise

def create_parser():
    parser = argparse.ArgumentParser(
        description=f"Bismark甲基化分析流程|Bismark Methylation Analysis Pipeline v{__version__}",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -g /path/to/genome.fa -r /path/to/raw_data -o /path/to/output_dir
        '''
    )

    req = parser.add_argument_group("必需参数|Required arguments")
    req.add_argument('-g', '--genome', required=True, help='基因组FASTA文件路径|Genome FASTA file')
    req.add_argument('-i', '--input', required=True, help='原始FASTQ数据目录|Raw FASTQ data directory')
    req.add_argument('-o', '--output-dir', required=True, help='主输出目录|Main output directory')

    opt = parser.add_argument_group("可选参数|Optional arguments")
    opt.add_argument('-p', '--pattern', default='_1_clean.fq.gz', help='R1文件的后缀模式|Suffix pattern for R1 files')
    opt.add_argument('-t', '--threads', type=int, default=88, help='使用的线程数|Number of threads')
    opt.add_argument('--sort-buffer', type=str, default='400G', help='提取步骤的排序缓存大小|Sort buffer size')
    opt.add_argument('--no-no-overlap', dest='no_overlap', action='store_false', help='不忽略重叠的reads|Do NOT ignore overlapping reads')

    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()

    try:
        analyzer = BismarkAnalyzer(
            raw_dir=args.input,
            genome_fa=args.genome,
            output_dir=args.output_dir,
            threads=args.threads,
            sort_buffer=args.sort_buffer,
            no_overlap=args.no_overlap,
            pattern=args.pattern
        )
        analyzer.run_analysis()

    except Exception as e:
        print(f"\n分析失败|Analysis failed: {e}", file=sys.stderr)
        sys.exit(1)
