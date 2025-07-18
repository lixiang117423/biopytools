"""
K-mer PAV分析主程序模块 (双阶段设计) | K-mer PAV Analysis Main Module (Two-Stage Design)
"""

import argparse
import sys
from .config import KmerConfig
from .utils import KmerLogger, cleanup_files
from .database_builder import DatabaseBuilder
from .query_processor import QueryProcessor
from .matrix_builder import KmerMatrixBuilder, StatisticsCalculator
from .results import ResultsWriter, SummaryGenerator

class KmerPAVAnalyzer:
    """K-mer PAV分析主类 (双阶段设计) | Main K-mer PAV Analyzer Class (Two-Stage Design)"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = KmerConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = KmerLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化各个处理器 | Initialize processors
        self.database_builder = DatabaseBuilder(self.config, self.logger)
        self.query_processor = QueryProcessor(self.config, self.logger)
        self.matrix_builder = KmerMatrixBuilder(self.config, self.logger)
        self.stats_calc = StatisticsCalculator(self.logger)
        self.results_writer = ResultsWriter(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的K-mer PAV分析流程 | Run complete K-mer PAV analysis pipeline"""
        try:
            self.logger.info("=" * 90)
            self.logger.info("开始K-mer PAV分析 (双阶段设计) | Starting K-mer PAV analysis (Two-Stage Design)")
            self.logger.info("=" * 90)
            
            # 阶段1: 构建k-mer数据库 | Phase 1: Build k-mer database
            database_name = self.database_builder.build_unified_database()
            
            # 阶段2: 处理查询样本 | Phase 2: Process query samples
            query_files = self.query_processor.get_query_files()
            samples = self.query_processor.prepare_samples(query_files)
            count_files, sample_names = self.query_processor.process_all_samples(samples, database_name)
            
            # 阶段3: 构建矩阵和统计 | Phase 3: Build matrices and statistics
            count_matrix, pa_matrix = self.matrix_builder.build_matrices(count_files, sample_names)
            stats = self.stats_calc.calculate_statistics(count_matrix, pa_matrix)
            
            # 保存结果 | Save results
            output_files = self.results_writer.save_matrices(count_matrix, pa_matrix)
            stats_file = self.results_writer.save_statistics(stats)
            output_files['statistics'] = stats_file
            
            # 生成总结报告 | Generate summary report
            self.summary_generator.generate_summary_report(output_files, stats)
            
            # 最终清理 | Final cleanup
            self._final_cleanup(database_name)
            
            self.logger.info("=" * 90)
            self.logger.info("K-mer PAV分析完成！| K-mer PAV analysis completed!")
            self.logger.info("=" * 90)
            self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
            self.logger.info(f"主要输出文件 | Main output files:")
            self.logger.info(f"  - {self.config.output_prefix}_counts.csv: 计数矩阵 | Count matrix")
            self.logger.info(f"  - {self.config.output_prefix}_presence_absence.csv: 存在缺失矩阵 | Presence/absence matrix")
            
        except Exception as e:
            self.logger.error(f"分析过程中发生错误 | Error occurred during analysis: {e}")
            sys.exit(1)
    
    def _final_cleanup(self, database_name: str):
        """最终清理工作 | Final cleanup"""
        if not self.config.keep_intermediate:
            self.logger.info("执行最终清理 | Performing final cleanup")
            
            # 清理数据库文件 | Clean up database files
            cleanup_files([f"{database_name}.kmc*"], self.config.output_path, self.logger)
            
            # 清理分割的序列文件 | Clean up split sequence files
            seq_dir = self.config.output_path / "individual_sequences"
            if seq_dir.exists():
                try:
                    import shutil
                    shutil.rmtree(seq_dir)
                    self.logger.info(f"删除临时序列目录 | Removed temporary sequence directory: {seq_dir}")
                except Exception as e:
                    self.logger.warning(f"删除临时序列目录失败 | Failed to remove temporary sequence directory: {e}")
            
            # 清理临时目录 | Clean up temporary directory
            if self.config.tmp_path.exists():
                try:
                    import shutil
                    shutil.rmtree(self.config.tmp_path)
                    self.logger.info(f"删除临时目录 | Removed temporary directory: {self.config.tmp_path}")
                except Exception as e:
                    self.logger.warning(f"删除临时目录失败 | Failed to remove temporary directory: {e}")

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="K-mer PAV (Presence/Absence Variation) 分析工具 (双阶段设计) | K-mer PAV Analysis Tool (Two-Stage Design)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
双阶段设计说明 | Two-Stage Design Description:
  阶段1 | Phase 1: 从数据库文件构建统一k-mer数据库
  阶段2 | Phase 2: 查询文件与k-mer数据库比较分析
  
样本处理逻辑 | Sample Processing Logic:
  - FASTQ文件: 整个文件作为一个样本
  - FASTA文件: 每条序列作为一个样本
  
依赖工具 | Required Tools:
  - KMC: K-mer计数工具 | K-mer counting tool
  - kmc_tools: KMC工具套件 | KMC tools suite
  - BioPython: FASTA序列处理 | FASTA sequence processing
  
安装方法 | Installation:
  Ubuntu/Debian: sudo apt-get install kmc
  Conda: conda install -c bioconda kmc
  Python: pip install biopython pandas numpy

示例 | Examples:
  %(prog)s --database-input db_files/ --query-input samples/ -o analysis
  %(prog)s --database-input db_files/ --query-input query.fasta -s 25
  %(prog)s --database-input ref.fasta --query-input samples/ --threads 16
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument("--database-input", required=True,
                       help="数据库输入文件/目录 (用于构建k-mer数据库) | Database input file/directory (for building k-mer database)")
    parser.add_argument("--query-input", required=True,
                       help="查询输入文件/目录 (用于样本比较) | Query input file/directory (for sample comparison)")
    
    # 输出参数 | Output arguments
    parser.add_argument("-o", "--output-prefix", default="kmer_analysis",
                       help="输出文件前缀 | Output file prefix")
    parser.add_argument("--output-dir", default="./kmer_output",
                       help="输出目录 | Output directory")
    
    # K-mer参数 | K-mer arguments
    parser.add_argument("-s", "--size", type=int, default=31,
                       help="K-mer大小 | K-mer size")
    parser.add_argument("-r", "--reverse-complement", action="store_true",
                       help="包含反向互补序列 | Include reverse complement sequences")
    parser.add_argument("--min-count", type=int, default=1,
                       help="最小k-mer计数阈值 | Minimum k-mer count threshold")
    parser.add_argument("--max-count", type=int, default=1000000,
                       help="最大k-mer计数阈值 | Maximum k-mer count threshold")
    
    # 文件处理参数 | File processing arguments
    parser.add_argument("--database-pattern",
                       help="数据库文件匹配模式 | Database file matching pattern")
    parser.add_argument("--query-pattern",
                       help="查询文件匹配模式 | Query file matching pattern")
    parser.add_argument("-t", "--threads", type=int, default=8,
                       help="线程数 | Number of threads")
    
    # KMC特定参数 | KMC specific arguments
    parser.add_argument("--kmc-memory", type=int, default=16,
                       help="KMC内存限制(GB) | KMC memory limit (GB)")
    parser.add_argument("--kmc-tmp-dir", default="kmc_tmp",
                       help="KMC临时目录 | KMC temporary directory")
    parser.add_argument("--keep-intermediate", action="store_true",
                       help="保留中间文件 | Keep intermediate files")
    
    # 工具路径 | Tool paths
    parser.add_argument("--kmc-path", default="kmc",
                       help="KMC可执行文件路径 | KMC executable path")
    parser.add_argument("--kmc-tools-path", default="kmc_tools",
                       help="kmc_tools可执行文件路径 | kmc_tools executable path")
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = KmerPAVAnalyzer(
        database_input=args.database_input,
        query_input=args.query_input,
        output_prefix=args.output_prefix,
        output_dir=args.output_dir,
        kmer_size=args.size,
        min_count=args.min_count,
        max_count=args.max_count,
        reverse_complement=args.reverse_complement,
        database_pattern=args.database_pattern,
        query_pattern=args.query_pattern,
        threads=args.threads,
        kmc_memory_gb=args.kmc_memory,
        kmc_tmp_dir=args.kmc_tmp_dir,
        keep_intermediate=args.keep_intermediate,
        kmc_path=args.kmc_path,
        kmc_tools_path=args.kmc_tools_path
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
