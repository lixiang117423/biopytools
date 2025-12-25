"""
🌳 RAxML系统发育分析主程序模块 | RAxML Phylogenetic Analysis Main Module
"""

import argparse
import sys
from .config import RAxMLConfig
from .utils import RAxMLLogger, CommandRunner, check_dependencies
from .preprocessing import SequenceProcessor
from .analysis import RAxMLAnalyzer
from .results import ResultsManager

class RAxMLPhylogeneticAnalyzer:
    """RAxML系统发育分析主类 | Main RAxML Phylogenetic Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = RAxMLConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = RAxMLLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.sequence_processor = SequenceProcessor(self.config, self.logger, self.cmd_runner)
        self.raxml_analyzer = RAxMLAnalyzer(self.config, self.logger, self.cmd_runner)
        self.results_manager = ResultsManager(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的RAxML分析流程 | Run complete RAxML analysis pipeline"""
        try:
            self.logger.info("🌳" + "=" * 77)
            self.logger.info("🌳 开始RAxML系统发育分析 | Starting RAxML phylogenetic analysis")
            self.logger.info("🌳" + "=" * 77)
            
            # 步骤0: 检查依赖和环境 | Step 0: Check dependencies and environment
            self.logger.info("-" * 50)
            self.logger.info("🔍 步骤0: 环境检查 | Step 0: Environment check")
            self.logger.info("-" * 50)
            
            if not self.check_dependencies():
                raise RuntimeError("依赖检查失败 | Dependency check failed")
            
            # 步骤1: 序列文件预处理 | Step 1: Sequence file preprocessing
            self.logger.info("-" * 50)
            self.logger.info("📊 步骤1: 序列文件预处理 | Step 1: Sequence file preprocessing")
            self.logger.info("-" * 50)
            
            if not self.sequence_processor.validate_sequence_format():
                raise RuntimeError("序列文件格式验证失败 | Sequence file format validation failed")
            
            # 获取序列统计信息 | Get sequence statistics
            seq_stats = self.sequence_processor.get_sequence_statistics()
            
            # 准备工作目录 | Prepare working directory
            self.sequence_processor.prepare_working_directory()
            
            # 步骤2: 运行RAxML分析 | Step 2: Run RAxML analysis
            self.logger.info("-" * 50)
            self.logger.info("🌳 步骤2: RAxML系统发育分析 | Step 2: RAxML phylogenetic analysis")
            self.logger.info("-" * 50)
            
            if not self.raxml_analyzer.run_raxml_analysis():
                raise RuntimeError("RAxML分析失败 | RAxML analysis failed")
            
            # 步骤3: 结果处理和验证 | Step 3: Results processing and validation
            self.logger.info("-" * 50)
            self.logger.info("📋 步骤3: 结果处理 | Step 3: Results processing")
            self.logger.info("-" * 50)
            
            if not self.results_manager.validate_results():
                self.logger.warning("⚠️ 结果验证发现问题 | Issues found during result validation")
            
            # 生成总结报告 | Generate summary report
            summary_report = self.results_manager.generate_summary_report()
            
            self.logger.info("🌳" + "=" * 77)
            self.logger.info("🎉 RAxML系统发育分析完成 | RAxML phylogenetic analysis completed successfully")
            self.logger.info("🌳" + "=" * 77)
            
            # 输出重要文件路径 | Output important file paths
            output_files = self.results_manager.collect_output_files()
            if 'best_tree' in output_files:
                self.logger.info(f"🌟 最佳树文件 | Best tree file: {output_files['best_tree']}")
            if 'bipartitions' in output_files:
                self.logger.info(f"🔢 支持值树文件 | Support values tree file: {output_files['bipartitions']}")
            
            self.logger.info(f"📊 详细报告 | Detailed report: {summary_report}")
            
        except Exception as e:
            self.logger.error(f"❌ 分析流程失败 | Analysis pipeline failed: {e}")
            raise

def create_parser():
    """创建命令行参数解析器 | Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='🌳 RAxML系统发育分析工具 v1.0 | RAxML Phylogenetic Analysis Tool v1.0',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
📖 使用示例 | Usage Examples:
  %(prog)s -s alignment.phy -n my_tree
  %(prog)s -s data.phy -n phylo -m GTRGAMMA -f a -x 12345 -# 1000
  %(prog)s -s proteins.phy -n prot_tree -m PROTGAMMAWAG -T 64
  %(prog)s -s alignment.phy -n constrained_tree -g constraint.tre -T 88
        """
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('📋 必需参数 | Required arguments')
    required.add_argument('-s', '--sequence-file', required=True, 
                         help='🧬 输入序列文件 (PHYLIP格式) | Input sequence file (PHYLIP format)')
    required.add_argument('-n', '--output-name', required=True,
                         help='📄 输出文件名称 | Output file name')
    
    # 模型参数 | Model parameters
    model_group = parser.add_argument_group('🧮 模型参数 | Model parameters')
    model_group.add_argument('-m', '--model', default='GTRGAMMA',
                            help='🔬 替换模型 | Substitution model (GTRGAMMA, PROTGAMMAWAG, etc.)')
    model_group.add_argument('-c', '--categories', type=int, default=25,
                            help='📊 速率异质性类别数 | Number of rate heterogeneity categories')
    model_group.add_argument('-e', '--likelihood-epsilon', type=float, default=0.1,
                            help='🎯 似然优化精度 | Likelihood optimization precision')
    
    # 算法参数 | Algorithm parameters
    algo_group = parser.add_argument_group('⚙️ 算法参数 | Algorithm parameters')
    algo_group.add_argument('-f', '--algorithm', default='d',
                           help='🔄 算法类型 | Algorithm type (d=rapid hill-climbing, a=rapid bootstrap, etc.)')
    algo_group.add_argument('-p', '--parsimony-seed', type=int,
                           help='🎲 简约法随机种子 | Parsimony random seed')
    algo_group.add_argument('-#', '-N', '--runs', default='1',
                           help='🔢 运行次数或bootstrap次数 | Number of runs or bootstrap replicates')
    
    # Bootstrap参数 | Bootstrap parameters
    bootstrap_group = parser.add_argument_group('🔄 Bootstrap参数 | Bootstrap parameters')
    bootstrap_group.add_argument('-b', '--bootstrap-seed', type=int,
                                help='🎲 Bootstrap随机种子 | Bootstrap random seed')
    bootstrap_group.add_argument('-x', '--rapid-bootstrap-seed', type=int,
                                help='⚡ 快速bootstrap随机种子 | Rapid bootstrap random seed')
    bootstrap_group.add_argument('-I', '--bootstrap-convergence',
                                choices=['autoFC', 'autoMR', 'autoMRE', 'autoMRE_IGN'],
                                help='🛑 Bootstrap收敛标准 | Bootstrap convergence criterion')
    bootstrap_group.add_argument('-B', '--bootstop-threshold', type=float, default=0.03,
                                help='🚦 Bootstrap停止阈值 | Bootstrap stop threshold')
    bootstrap_group.add_argument('--bootstop-perms', type=int, default=100,
                                help='🔢 Bootstrap停止检验次数 | Bootstrap stop test permutations')
    bootstrap_group.add_argument('-k', '--print-bootstrap-trees', action='store_true',
                                help='📊 输出带分支长度的bootstrap树 | Print bootstrap trees with branch lengths')
    
    # 树参数 | Tree parameters
    tree_group = parser.add_argument_group('🌲 树参数 | Tree parameters')
    tree_group.add_argument('-t', '--starting-tree',
                           help='🌱 起始树文件 | Starting tree file')
    tree_group.add_argument('-g', '--constraint-tree',
                           help='🔗 约束树文件 | Constraint tree file')
    tree_group.add_argument('-o', '--outgroup',
                           help='🎯 外群名称 (逗号分隔多个) | Outgroup name(s) (comma-separated)')
    
    # 性能参数 | Performance parameters
    perf_group = parser.add_argument_group('⚡ 性能参数 | Performance parameters')
    perf_group.add_argument('-T', '--threads', type=int, default=88,
                           help='🔥 线程数 | Number of threads')
    perf_group.add_argument('-U', '--memory-saving', action='store_true',
                           help='💾 启用内存节省模式 | Enable memory saving mode')
    perf_group.add_argument('-D', '--ml-search-convergence', action='store_true',
                           help='🎯 启用ML搜索收敛标准 | Enable ML search convergence criterion')
    
    # 高级参数 | Advanced parameters
    advanced_group = parser.add_argument_group('🔬 高级参数 | Advanced parameters')
    advanced_group.add_argument('-d', '--random-starting-tree', action='store_true',
                                help='🎲 使用随机起始树 | Use random starting tree')
    advanced_group.add_argument('-V', '--disable-rate-heterogeneity', action='store_true',
                                help='🚫 禁用速率异质性模型 | Disable rate heterogeneity model')
    advanced_group.add_argument('-u', '--gamma-median', action='store_true',
                                help='📊 使用GAMMA模型中位数 | Use median for GAMMA model')
    advanced_group.add_argument('-H', '--disable-pattern-compression', action='store_true',
                                help='🚫 禁用模式压缩 | Disable pattern compression')
    
    # 输出和工具参数 | Output and tool parameters
    output_group = parser.add_argument_group('📂 输出和工具参数 | Output and tool parameters')
    output_group.add_argument('-w', '--output-dir', default='./raxml_output',
                             help='📁 输出目录 | Output directory')
    output_group.add_argument('--raxml-path', default='raxmlHPC-PTHREADS',
                             help='🔧 RAxML程序路径 | RAxML program path')
    
    # 质量控制参数 | Quality control parameters
    qc_group = parser.add_argument_group('🔍 质量控制参数 | Quality control parameters')
    qc_group.add_argument('--no-seq-check', action='store_true',
                         help='🚫 跳过序列检查 | Skip sequence checking')
    qc_group.add_argument('--silent', action='store_true',
                         help='🔇 静默模式 | Silent mode')
    
    return parser

def main():
    """主函数 | Main function"""
    parser = create_parser()
    args = parser.parse_args()
    
    # 处理特殊参数 | Handle special parameters
    rate_het_model = not args.disable_rate_heterogeneity
    pattern_compression = not args.disable_pattern_compression
    
    # 创建分析器并运行 | Create analyzer and run
    try:
        analyzer = RAxMLPhylogeneticAnalyzer(
            sequence_file=args.sequence_file,
            output_name=args.output_name,
            model=args.model,
            output_dir=args.output_dir,
            algorithm=args.algorithm,
            parsimony_seed=args.parsimony_seed,
            bootstrap_seed=args.bootstrap_seed,
            rapid_bootstrap_seed=args.rapid_bootstrap_seed,
            runs=args.runs,
            starting_tree=args.starting_tree,
            constraint_tree=args.constraint_tree,
            outgroup=args.outgroup,
            categories=args.categories,
            rate_het_model=rate_het_model,
            gamma_median=args.gamma_median,
            threads=args.threads,
            likelihood_epsilon=args.likelihood_epsilon,
            ml_search_convergence=args.ml_search_convergence,
            random_starting_tree=args.random_starting_tree,
            memory_saving=args.memory_saving,
            pattern_compression=pattern_compression,
            bootstrap_convergence=args.bootstrap_convergence,
            bootstop_threshold=args.bootstop_threshold,
            bootstop_perms=args.bootstop_perms,
            print_bootstrap_trees=args.print_bootstrap_trees,
            raxml_path=args.raxml_path,
            no_seq_check=args.no_seq_check,
            silent=args.silent
        )
        
        analyzer.run_analysis()
        
    except KeyboardInterrupt:
        print("\n⚠️ 分析被用户中断 | Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 分析失败 | Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
