"""
MSA主程序模块 | MSA Main Module
"""

import argparse
import sys
from .config import MSAConfig
from .utils import MSALogger, CommandRunner, check_dependencies, count_sequences
from .aligner import MSAAligner
from .stats import AlignmentStats

class MSAAnalyzer:
    """MSA分析主类 | Main MSA Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = MSAConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = MSALogger(self.config.output_prefix)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)
        
        # 初始化各个处理器 | Initialize processors
        self.aligner = MSAAligner(self.config, self.logger, self.cmd_runner)
        self.stats_analyzer = AlignmentStats(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_alignment(self):
        """运行完整的比对流程 | Run complete alignment pipeline"""
        try:
            self.logger.info("="*60)
            self.logger.info("🧬 多序列比对分析开始 | MSA Analysis Started")
            self.logger.info("="*60)
            
            # 1. 检查依赖 | Check dependencies
            self.check_dependencies()
            
            # 2. 统计输入序列 | Count input sequences
            seq_count = count_sequences(self.config.input_file)
            self.logger.info(f"📊 输入序列数量 | Input sequences: {seq_count}")
            
            # 3. 执行比对 | Perform alignment
            if not self.aligner.align():
                raise RuntimeError("❌ 比对失败 | Alignment failed")
            
            # 4. 计算统计信息 | Calculate statistics
            self.stats_analyzer.calculate_stats()
            
            # 5. 完成 | Complete
            self.logger.info("="*60)
            self.logger.info("✅ 多序列比对分析完成 | MSA Analysis Completed")
            self.logger.info("="*60)
            self.logger.info(f"📁 结果保存在 | Results saved: {self.config.output_prefix}.*")
            
        except Exception as e:
            self.logger.error(f"❌ 分析流程异常终止 | Analysis terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 多序列比对工具 (模块化版本) | Multiple Sequence Alignment Tool (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s -i sequences.fasta -o alignment
  %(prog)s -i seqs.fa -o result -m clustalo -t 16
  %(prog)s -i input.fasta -o out -m mafft --mafft-strategy linsi
  %(prog)s -i seqs.fa -o align -m muscle --muscle-maxiters 32
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', required=True, 
                       help='📁 输入序列文件 (FASTA格式) | Input sequence file (FASTA format)')
    parser.add_argument('-o', '--output', required=True,
                       help='📝 输出文件前缀 | Output file prefix')
    
    # 比对方法 | Alignment method
    parser.add_argument('-m', '--method', default='mafft',
                       choices=['mafft', 'clustalo', 'muscle', 't_coffee'],
                       help='🔧 比对方法 | Alignment method')
    
    # 通用参数 | Common parameters
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='🧵 线程数 | Number of threads')
    parser.add_argument('-f', '--format', default='fasta',
                       choices=['fasta', 'clustal', 'phylip', 'nexus'],
                       help='💾 输出格式 | Output format')
    
    # MAFFT参数 | MAFFT parameters
    mafft_group = parser.add_argument_group('MAFFT参数 | MAFFT Parameters')
    mafft_group.add_argument('--mafft-strategy', default='auto',
                            choices=['auto', 'linsi', 'ginsi', 'einsi', 'fftns', 'fftnsi'],
                            help='🎯 MAFFT策略 | MAFFT strategy')
    mafft_group.add_argument('--mafft-maxiterate', type=int, default=1000,
                            help='🔄 MAFFT最大迭代次数 | MAFFT max iterations')
    
    # Clustal Omega参数 | Clustal Omega parameters
    clustalo_group = parser.add_argument_group('Clustal Omega参数 | Clustal Omega Parameters')
    clustalo_group.add_argument('--clustalo-iterations', type=int, default=0,
                               help='🔄 Clustal Omega迭代次数 | Clustal Omega iterations')
    
    # MUSCLE参数 | MUSCLE parameters
    muscle_group = parser.add_argument_group('MUSCLE参数 | MUSCLE Parameters')
    muscle_group.add_argument('--muscle-maxiters', type=int, default=16,
                             help='🔄 MUSCLE最大迭代次数 | MUSCLE max iterations')
    
    # 工具路径 | Tool paths
    tool_group = parser.add_argument_group('工具路径 | Tool Paths')
    tool_group.add_argument('--mafft-path', default='mafft',
                           help='📍 MAFFT路径 | MAFFT path')
    tool_group.add_argument('--clustalo-path', default='clustalo',
                           help='📍 Clustal Omega路径 | Clustal Omega path')
    tool_group.add_argument('--muscle-path', default='muscle',
                           help='📍 MUSCLE路径 | MUSCLE path')
    tool_group.add_argument('--tcoffee-path', default='t_coffee',
                           help='📍 T-Coffee路径 | T-Coffee path')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = MSAAnalyzer(
        input_file=args.input,
        output_prefix=args.output,
        method=args.method,
        threads=args.threads,
        output_format=args.format,
        mafft_strategy=args.mafft_strategy,
        mafft_maxiterate=args.mafft_maxiterate,
        clustalo_iterations=args.clustalo_iterations,
        muscle_maxiters=args.muscle_maxiters,
        mafft_path=args.mafft_path,
        clustalo_path=args.clustalo_path,
        muscle_path=args.muscle_path,
        tcoffee_path=args.tcoffee_path
    )
    
    analyzer.run_alignment()

if __name__ == "__main__":
    main()
