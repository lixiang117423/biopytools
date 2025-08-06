"""
Bismark流程主程序模块 | Bismark Pipeline Main Module
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
            self.logger.info(f"🚀 开始Bismark甲基化分析流程 | Starting Bismark methylation analysis pipeline (v{__version__})")
            self.logger.info(f"主输出目录 | Main output directory: {self.config.output_dir}")
            check_dependencies(self.config, self.logger)
            
            if not self.core_pipeline.build_bismark_index():
                raise RuntimeError("❌ Bismark索引构建失败 | Bismark index build failed")
            
            if not self.core_pipeline.run_mapping_and_extraction():
                raise RuntimeError("❌ 甲基化比对和提取失败 | Methylation mapping and extraction failed")
            
            self.core_pipeline.results_manager.generate_summary_report()
            self.logger.info("🎉 === Bismark分析流程全部完成 | Bismark analysis pipeline completed successfully ===")
            
        except Exception as e:
            self.logger.error(f"❌ 分析流程意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            raise

def create_parser():
    parser = argparse.ArgumentParser(
        description=f"Bismark甲基化分析流程 | Bismark Methylation Analysis Pipeline v{__version__}",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
示例 | Examples:
  # 基础用法 | Basic usage
  %(prog)s -g /path/to/genome.fa -r /path/to/raw_data -o /path/to/output_dir

  # 使用默认值运行 (线程88, 内存400G, 模式_1_clean.fq.gz) | Run with defaults
  %(prog)s -g genome.fa -r ./cleandata -o ./results

目录结构说明 | Directory Structure Explanation:
  /path/to/genome_dir/       (基因组所在目录 | Directory containing the genome)
    ├── Bisulfite_Genome/    (自动创建的索引子目录 | Auto-created index subdirectory)
    └── genome.fa            (你的基因组文件 | Your genome file)

  /path/to/output_dir/       (你指定的输出目录 | Your specified output directory)
    ├── mapping/             (存放最终的BAM文件 | Stores final BAM files)
    ├── result/              (存放最终的甲基化结果, 包括拆分文件 | Stores final methylation results, including split files)
    ├── tmp/                 (存放所有中间文件 | Stores all intermediate files)
    └── logs/                (存放运行日志 | Stores run logs)
"""
    )
    
    req = parser.add_argument_group("必需参数 | Required Arguments")
    req.add_argument('-g', '--genome-fa', required=True, help='🧬 基因组FASTA文件路径 | Path to genome FASTA file')
    req.add_argument('-r', '--raw-dir', required=True, help='📁 原始FASTQ数据目录 | Raw FASTQ data directory')
    req.add_argument('-o', '--output-dir', required=True, help='📂 主输出目录 | Main output directory')
    
    opt = parser.add_argument_group("可选参数 | Optional Arguments")
    opt.add_argument('-p', '--pattern', default='_1_clean.fq.gz', help='R1文件的后缀模式 | Suffix pattern for R1 files (default: "_1_clean.fq.gz")')
    opt.add_argument('-j', '--threads', type=int, default=88, help='🔧 使用的线程数 | Number of threads to use (default: 88)')
    opt.add_argument('--sort-buffer', type=str, default='400G', help='💾 提取步骤的排序缓存大小 | Sort buffer size for extraction step (default: 400G)')
    opt.add_argument('--no-no-overlap', dest='no_overlap', action='store_false', help='🧪 在提取时【不】忽略重叠的reads (默认忽略) | Do NOT ignore overlapping reads during extraction (default is to ignore)')
    
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    
    try:
        analyzer = BismarkAnalyzer(
            raw_dir=args.raw_dir,
            genome_fa=args.genome_fa,
            output_dir=args.output_dir,
            threads=args.threads,
            sort_buffer=args.sort_buffer,
            no_overlap=args.no_overlap,
            pattern=args.pattern
        )
        analyzer.run_analysis()
        
    except Exception as e:
        print(f"\n❌ 分析失败 | Analysis failed: {e}", file=sys.stderr)
        sys.exit(1)
