"""
🚀 MetaGraph K-mer分析主程序模块 | MetaGraph K-mer Analysis Main Module
"""

import sys
import argparse
from .config import KmerConfig
from .logger import KmerLogger
from .utils import CommandRunner, FileFormatDetector, check_dependencies
from .metagraph_builder import MetaGraphBuilder
from .kmer_query import KmerQueryAnalyzer
from .merger import ResultsMerger

class KmerAnalyzer:
    """MetaGraph K-mer分析主类 | Main MetaGraph K-mer Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = KmerConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = KmerLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.metagraph_builder = MetaGraphBuilder(self.config, self.logger, self.cmd_runner)
        self.query_analyzer = KmerQueryAnalyzer(self.config, self.logger, self.cmd_runner)
        self.results_merger = ResultsMerger(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的MetaGraph k-mer分析流程 | Run complete MetaGraph k-mer analysis pipeline"""
        try:
            self.logger.info("=" * 80)
            self.logger.info("🧬 MetaGraph K-mer库构建与查询分析流程启动")
            self.logger.info("MetaGraph K-mer Library Construction and Query Analysis Pipeline Started")
            self.logger.info("=" * 80)
            
            # 检查依赖 | Check dependencies
            self.logger.info("🔍 步骤0: 检查依赖软件 | Step 0: Checking Dependencies")
            self.check_dependencies()
            
            # 获取参考文件 | Get reference files
            ref_files = FileFormatDetector.get_all_files(self.config.reference)
            if not ref_files:
                raise ValueError(f"❌ 未找到参考文件 | No reference files found: {self.config.reference}")
            
            self.logger.info(f"📁 找到 {len(ref_files)} 个参考文件 | Found {len(ref_files)} reference file(s)")
            
            # 检测文件格式 | Detect file format
            file_format = FileFormatDetector.detect_format(ref_files[0])
            self.logger.info(f"📄 文件格式 | File format: {file_format.upper()}")
            
            # 步骤1: 构建MetaGraph库 | Step 1: Build MetaGraph library
            self.logger.info("🏗️  步骤1: 构建MetaGraph参考库 | Step 1: Building MetaGraph Reference Library")
            self.logger.info(f"   📁 参考输入 | Reference input: {self.config.reference}")
            self.logger.info(f"   🔢 k-mer长度 | k-mer length: {self.config.kmer_length}")
            self.logger.info(f"   🧵 线程数 | Threads: {self.config.threads}")
            self.logger.info(f"   🔄 模式 | Mode: {'canonical' if self.config.canonical else 'basic'}")
            
            # 构建图 | Build graph
            self.metagraph_builder.build_graph(ref_files, file_format)
            
            # 添加坐标注释（仅FASTA） | Add coordinate annotations (FASTA only)
            coords_file = None
            if file_format == 'fasta':
                has_coords = self.metagraph_builder.annotate_with_coordinates(ref_files, file_format)
                
                if has_coords:
                    # 方案1: 使用MetaGraph query提取（可能有问题，使用方案2）
                    # coords_file = self.metagraph_builder.query_coordinates(ref_files)
                    
                    # 方案2: 直接从FASTA提取（更可靠）
                    self.logger.info("📝 使用直接提取方法获取k-mer坐标 | Using direct extraction method for k-mer coordinates")
                    coords_file = self.metagraph_builder.extract_kmers_from_fasta(ref_files)
            else:
                self.logger.warning("⚠️  FASTQ格式不支持坐标信息，将只输出k-mer和丰度")
                self.logger.warning("FASTQ format doesn't support coordinates, will output k-mer and abundance only")
                
                # FASTQ情况下，直接提取唯一k-mer列表（不带坐标）
                # 这里可以用MetaGraph dump或者KMC来做
                self.logger.info("使用MetaGraph构建完成，但FASTQ不保存坐标")
                # 暂时跳过，直接进行查询统计
            
            # 步骤2: 统计查询文件 | Step 2: Count query file
            self.logger.info("📊 步骤2: 统计查询文件k-mer丰度 | Step 2: Counting K-mer Abundance in Query")
            self.logger.info(f"   📁 查询输入 | Query input: {self.config.query}")
            
            counts_file = self.query_analyzer.count_kmers()
            
            # 步骤3: 合并结果 | Step 3: Merge results
            if coords_file:
                self.logger.info("🔗 步骤3: 合并结果 | Step 3: Merging Results")
                output_file = self.results_merger.merge_results(coords_file, counts_file)
            else:
                self.logger.warning("⚠️  无坐标信息，仅保存丰度统计")
                output_file = counts_file
            
            # 完成 | Complete
            self.logger.info("=" * 80)
            self.logger.info("🎉 分析流程完成 | Analysis Pipeline Completed Successfully!")
            self.logger.info("=" * 80)
            self.logger.info(f"📊 最终输出文件 | Final output file: {output_file}")
            self.logger.info(f"📁 输出目录 | Output directory: {self.config.output_dir}")
            self.logger.info(f"🗂️  MetaGraph图文件 | MetaGraph graph file: {self.config.dbg_file}")
            if file_format == 'fasta':
                self.logger.info(f"📍 坐标注释文件 | Coordinate annotation file: {self.config.anno_file}")
            self.logger.info("=" * 80)
            
        except Exception as e:
            self.logger.error(f"❌ 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly")
            self.logger.error(f"错误信息 | Error message: {str(e)}")
            import traceback
            self.logger.error(f"详细错误 | Detailed error:\n{traceback.format_exc()}")
            sys.exit(1)


def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 MetaGraph K-mer库构建与查询分析工具 | MetaGraph K-mer Library Construction and Query Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s -r reference.fasta -q query.fastq -o kmer_results
  %(prog)s -r ref_folder/ -q reads.fastq.gz -k 21 -t 64
  %(prog)s --reference genome.fa --query sample.fq --kmer 31 --threads 88 --output results/
  
更多信息 | More information:
  - 使用MetaGraph构建高效的k-mer图索引 | Uses MetaGraph to build efficient k-mer graph index
  - FASTA输入会保存坐标信息 | FASTA input saves coordinate information
  - 使用KMC统计k-mer丰度 | Uses KMC to count k-mer abundance
  - 自动处理canonical k-mer | Automatically handles canonical k-mers
  
文档参考 | Documentation:
  - MetaGraph: https://github.com/ratschlab/metagraph
  - KMC: https://github.com/refresh-bio/KMC
        """
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('必需参数 | Required Arguments')
    required.add_argument('-r', '--reference', required=True,
                         help='🧬 参考序列文件或文件夹 | Reference sequence file or folder (FASTA/FASTQ)')
    required.add_argument('-q', '--query', required=True,
                         help='🔍 查询序列文件 | Query sequence file (FASTA/FASTQ)')
    
    # 输出参数 | Output arguments
    output = parser.add_argument_group('输出参数 | Output Arguments')
    output.add_argument('-o', '--output', default='./kmer_output',
                       help='📁 输出目录 | Output directory (default: ./kmer_output)')
    
    # K-mer参数 | K-mer parameters
    kmer_params = parser.add_argument_group('K-mer参数 | K-mer Parameters')
    kmer_params.add_argument('-k', '--kmer', type=int, default=31,
                            help='🔢 k-mer长度 | K-mer length (default: 31)')
    kmer_params.add_argument('--no-canonical', action='store_true',
                            help='❌ 禁用canonical k-mer处理 | Disable canonical k-mer processing')
    kmer_params.add_argument('--min-count', type=int, default=1,
                            help='📊 最小k-mer计数阈值 | Minimum k-mer count threshold (default: 1)')
    
    # 性能参数 | Performance parameters
    performance = parser.add_argument_group('性能参数 | Performance Parameters')
    performance.add_argument('-t', '--threads', type=int, default=88,
                           help='🧵 线程数 | Number of threads (default: 88)')
    performance.add_argument('-m', '--memory', type=int, default=64,
                           help='💾 KMC内存限制(GB) | KMC memory limit in GB (default: 64)')
    
    # 工具路径 | Tool paths
    tools = parser.add_argument_group('工具路径 | Tool Paths')
    tools.add_argument('--metagraph-path', default='metagraph',
                      help='🔧 MetaGraph工具路径 | MetaGraph tool path (default: metagraph)')
    tools.add_argument('--kmc-path', default='kmc',
                      help='🔧 KMC工具路径 | KMC tool path (default: kmc)')
    tools.add_argument('--kmc-tools-path', default='kmc_tools',
                      help='🔧 KMC-tools路径 | KMC-tools path (default: kmc_tools)')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = KmerAnalyzer(
        reference=args.reference,
        query=args.query,
        output_dir=args.output,
        kmer_length=args.kmer,
        threads=args.threads,
        memory_gb=args.memory,
        canonical=not args.no_canonical,
        min_count=args.min_count,
        metagraph_path=args.metagraph_path,
        kmc_path=args.kmc_path,
        kmc_tools_path=args.kmc_tools_path
    )
    
    analyzer.run_analysis()


if __name__ == "__main__":
    main()
