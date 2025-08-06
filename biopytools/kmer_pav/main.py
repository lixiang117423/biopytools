"""
🚀 K-mer PAV分析主模块 | K-mer PAV Analysis Main Module
"""

import argparse
import sys
from pathlib import Path

from .config import PAVConfig
from .utils import PAVLogger, CommandRunner, check_dependencies
from .data_processing import GenomeKmerExtractor, SampleKmerExtractor, PresenceAnalyzer
from .matrix_generator import PresenceMatrixGenerator, WindowAnalyzer
from .reporter import AnalysisReporter

class KmerPAVAnalyzer:
    """🧬 K-mer PAV分析器 | K-mer PAV Analyzer"""
    
    def __init__(self, genome_file, fastq_dir, output_dir='./kmer_pav_output', 
                 kmer_size=51, threads=8, canonical=True, sort_output=True,
                 fastq_pattern="*_1.fq.gz", window_size=500000, step_size=None,
                 unikmer_path='unikmer'):
        """
        🔧 初始化分析器 | Initialize analyzer
        
        Args:
            genome_file: 🧬 基因组FASTA文件路径 | Genome FASTA file path
            fastq_dir: 📂 FASTQ文件目录 | FASTQ files directory
            output_dir: 📁 输出目录 | Output directory
            kmer_size: 🧮 K-mer长度 | K-mer size
            threads: 🔧 线程数 | Number of threads
            canonical: ✅ 是否使用canonical k-mer | Whether to use canonical k-mers
            sort_output: 📊 是否排序输出 | Whether to sort output
            fastq_pattern: 🎯 FASTQ文件匹配模式 | FASTQ file pattern
            window_size: 📏 窗口大小(bp) | Window size (bp)
            step_size: 👣 步长(bp)，None表示非重叠 | Step size (bp), None means non-overlapping
            unikmer_path: 🛠️ unikmer工具路径 | unikmer tool path
        """
        # 创建配置 | Create configuration
        self.config = PAVConfig(
            genome_file=genome_file,
            fastq_dir=fastq_dir,
            output_dir=output_dir,
            kmer_size=kmer_size,
            threads=threads,
            canonical=canonical,
            sort_output=sort_output,
            fastq_pattern=fastq_pattern,
            window_size=window_size,
            step_size=step_size,
            unikmer_path=unikmer_path
        )
        
        # 验证配置 | Validate configuration
        self.config.validate()
        
        # 设置日志 | Setup logging
        self.logger_manager = PAVLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 设置命令执行器 | Setup command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化处理模块 | Initialize processing modules
        self.genome_extractor = GenomeKmerExtractor(self.config, self.logger, self.cmd_runner)
        self.sample_extractor = SampleKmerExtractor(self.config, self.logger, self.cmd_runner)
        self.presence_analyzer = PresenceAnalyzer(self.config, self.logger, self.cmd_runner)
        self.matrix_generator = PresenceMatrixGenerator(self.config, self.logger)
        self.window_analyzer = WindowAnalyzer(self.config, self.logger)
        self.reporter = AnalysisReporter(self.config, self.logger)
    
    def run_analysis(self):
        """🚀 运行完整的K-mer PAV分析流程 | Run complete K-mer PAV analysis pipeline"""
        try:
            self.logger.info("🚀 开始K-mer PAV分析流程 | Starting K-mer PAV analysis pipeline")
            self.logger.info(f"🧬 基因组文件 | Genome file: {self.config.genome_file}")
            self.logger.info(f"📂 FASTQ目录 | FASTQ directory: {self.config.fastq_dir}")
            self.logger.info(f"📁 输出目录 | Output directory: {self.config.output_dir}")
            self.logger.info(f"🧮 K-mer长度 | K-mer size: {self.config.kmer_size}")
            
            # 检查依赖 | Check dependencies
            self.check_dependencies()
            
            # 步骤1: 从基因组提取k-mer | Step 1: Extract k-mers from genome
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("🧬 步骤1: 从基因组提取k-mer | Step 1: Extract k-mers from genome")
            self.logger.info(f"{'=' * 60}")
            
            if not self.genome_extractor.extract_genome_kmers():
                raise RuntimeError("❌ 基因组k-mer提取失败 | Genome k-mer extraction failed")
            
            # 步骤2: 获取k-mer位置信息 | Step 2: Get k-mer positions
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("📍 步骤2: 获取k-mer位置信息 | Step 2: Get k-mer positions")
            self.logger.info(f"{'=' * 60}")
            
            if not self.genome_extractor.locate_kmers():
                raise RuntimeError("❌ k-mer位置定位失败 | K-mer position location failed")
            
            # 步骤3: 从FASTQ样本提取k-mer | Step 3: Extract k-mers from FASTQ samples
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("🧪 步骤3: 从FASTQ样本提取k-mer | Step 3: Extract k-mers from FASTQ samples")
            self.logger.info(f"{'=' * 60}")
            
            samples, failed_extraction = self.sample_extractor.discover_and_process_samples()
            if not samples:
                raise RuntimeError("❌ 未发现任何样本 | No samples discovered")
            
            # 步骤4: 检查存在性 | Step 4: Check presence
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("🎯 步骤4: 检查基因组k-mer在样本中的存在情况 | Step 4: Check presence of genome k-mers in samples")
            self.logger.info(f"{'=' * 60}")
            
            failed_presence = self.presence_analyzer.analyze_presence(samples)
            
            # 步骤5: 生成存在性矩阵 | Step 5: Generate presence matrix
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("📊 步骤5: 生成存在性矩阵 | Step 5: Generate presence matrix")
            self.logger.info(f"{'=' * 60}")
            
            if not self.matrix_generator.generate_matrix():
                raise RuntimeError("❌ 存在性矩阵生成失败 | Presence matrix generation failed")
            
            # 步骤6: 窗口分析 | Step 6: Window analysis
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("🪟 步骤6: 窗口分析 | Step 6: Window analysis")
            self.logger.info(f"{'=' * 60}")
            
            if not self.window_analyzer.analyze_windows():
                self.logger.warning("⚠️ 窗口分析失败，继续进行报告生成 | Window analysis failed, continuing with report generation")
            
            # 步骤7: 生成分析报告 | Step 7: Generate analysis report
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("📝 步骤7: 生成分析报告 | Step 7: Generate analysis report")
            self.logger.info(f"{'=' * 60}")
            
            self.reporter.generate_report()
            
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("🎉 K-mer PAV分析完成！| K-mer PAV analysis completed!")
            self.logger.info(f"{'=' * 60}")
            
            # 显示失败的样本（如果有）| Show failed samples (if any)
            if failed_extraction or failed_presence:
                self.logger.warning("\n⚠️ 注意: 部分样本处理失败 | Note: Some samples failed processing")
                if failed_extraction:
                    self.logger.warning(f"❌ k-mer提取失败的样本: {', '.join(failed_extraction)} | Failed k-mer extraction: {', '.join(failed_extraction)}")
                if failed_presence:
                    self.logger.warning(f"❌ 存在性检查失败的样本: {', '.join(failed_presence)} | Failed presence check: {', '.join(failed_presence)}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 分析过程中出现错误: {e} | Error occurred during analysis: {e}")
            return False
    
    def check_dependencies(self):
        """🔍 检查依赖软件 | Check dependencies"""
        check_dependencies(self.config, self.logger)

def create_argument_parser():
    """⚙️ 创建命令行参数解析器 | Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='🧬 K-mer PAV (Presence/Absence Variation) 分析工具 | K-mer PAV Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
📖 示例 | Examples:
  # 🎯 基本分析 | Basic analysis
  python -m kmer_pav -g genome.fa -f /path/to/fastq -o results
  
  # ⚙️ 自定义参数 | Custom parameters  
  python -m kmer_pav -g genome_chr12.fa -f fastq_dir \\
                    -o kmer_analysis -k 31 -t 16 -p "*_R1.fastq.gz"
  
  # 🪟 窗口分析 | Window analysis
  python -m kmer_pav -g genome.fa -f fastq_dir \\
                    -o results -w 1000000 -s 200000
  
  # 🛠️ 指定unikmer路径 | Specify unikmer path
  python -m kmer_pav -g genome.fa -f fastq_data \\
                    -o output --unikmer /usr/local/bin/unikmer
"""
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('📋 必需参数 | Required arguments')
    required.add_argument('-g', '--genome', required=True,
                         help='🧬 基因组FASTA文件路径 | Genome FASTA file path')
    required.add_argument('-f', '--fastq-dir', required=True,
                         help='📂 FASTQ文件目录 | FASTQ files directory')
    
    # 输入参数 | Input arguments
    input_group = parser.add_argument_group('📥 输入参数 | Input arguments')
    input_group.add_argument('-p', '--pattern', default='*_1.fq.gz',
                           help='🎯 FASTQ文件匹配模式，*表示样本名 (默认: *_1.fq.gz) | FASTQ file pattern, * represents sample name (default: *_1.fq.gz)')
    
    # 输出参数 | Output arguments
    output = parser.add_argument_group('📤 输出参数 | Output arguments')
    output.add_argument('-o', '--output', default='./kmer_pav_output',
                       help='📁 输出目录 (默认: ./kmer_pav_output) | Output directory (default: ./kmer_pav_output)')
    
    # K-mer参数 | K-mer arguments
    kmer = parser.add_argument_group('🧮 K-mer参数 | K-mer arguments')
    kmer.add_argument('-k', '--kmer-size', type=int, default=51,
                     help='🧮 K-mer长度 (默认: 51) | K-mer size (default: 51)')
    kmer.add_argument('--no-canonical', action='store_true',
                     help='❌ 不使用canonical k-mer | Do not use canonical k-mers')
    kmer.add_argument('--no-sort', action='store_true',
                     help='❌ 不排序输出 | Do not sort output')
    
    # 窗口分析参数 | Window analysis arguments
    window = parser.add_argument_group('🪟 窗口分析参数 | Window analysis arguments')
    window.add_argument('-w', '--window', type=int, default=500000,
                       help='📏 窗口大小(bp) (默认: 500000) | Window size in bp (default: 500000)')
    window.add_argument('-s', '--step', type=int,
                       help='👣 步长(bp)，指定则使用重叠窗口 | Step size in bp, specify for overlapping windows')
    
    # 性能参数 | Performance arguments
    performance = parser.add_argument_group('⚡ 性能参数 | Performance arguments')
    performance.add_argument('-t', '--threads', type=int, default=8,
                           help='🔧 线程数 (默认: 8) | Number of threads (default: 8)')
    
    # 工具路径 | Tool paths
    tools = parser.add_argument_group('🛠️ 工具路径 | Tool paths')
    tools.add_argument('--unikmer', default='unikmer',
                      help='🛠️ unikmer工具路径 (默认: unikmer) | unikmer tool path (default: unikmer)')
    
    return parser

def main():
    """🚀 主函数 | Main function"""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # 创建分析器 | Create analyzer
    analyzer = KmerPAVAnalyzer(
        genome_file=args.genome,
        fastq_dir=args.fastq_dir,
        output_dir=args.output,
        kmer_size=args.kmer_size,
        threads=args.threads,
        canonical=not args.no_canonical,
        sort_output=not args.no_sort,
        fastq_pattern=args.pattern,
        window_size=args.window,
        step_size=args.step,
        unikmer_path=args.unikmer
    )
    
    # 运行分析 | Run analysis
    success = analyzer.run_analysis()
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
