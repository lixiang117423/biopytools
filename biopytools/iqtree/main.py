"""
IQ-TREE 分析主程序模块|IQ-TREE Analysis Main Module
"""

import argparse
import sys
from ..common.paths import expand_path
from .config import TreeConfig
from .utils import TreeLogger, CommandRunner, check_dependencies
from .analysis import TreeBuilder, ConcordanceAnalyzer, AncestralReconstructor
from .ancestral_parser import AncestralStateParser


class IQTreeAnalyzer:
    """IQ-TREE 系统发育树分析主类|Main IQ-TREE Phylogenetic Tree Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = TreeConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = TreeLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        # 初始化各个分析器|Initialize analyzers
        self.tree_builder = TreeBuilder(self.config, self.logger, self.cmd_runner)
        self.concordance_analyzer = ConcordanceAnalyzer(self.config, self.logger, self.cmd_runner)
        self.ancestral_reconstructor = AncestralReconstructor(self.config, self.logger, self.cmd_runner)
        self.ancestral_parser = AncestralStateParser(self.config, self.logger)

    def check_dependencies(self):
        """检查依赖软件|Check dependencies"""
        return check_dependencies(self.config, self.logger)

    def run_analysis(self):
        """运行完整的系统发育分析流程|Run complete phylogenetic analysis pipeline"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("IQ-TREE 系统发育树分析流程启动|IQ-TREE Phylogenetic Analysis Pipeline Started")
            self.logger.info("=" * 60)

            # 1. 检查依赖|Check dependencies
            self.logger.info("步骤1: 检查软件依赖|Step 1: Checking software dependencies")
            if not self.check_dependencies():
                self.logger.error("依赖检查失败|Dependency check failed")
                sys.exit(1)

            # 2. 系统发育树构建|Phylogenetic tree construction
            self.logger.info("步骤2: 构建系统发育树|Step 2: Building phylogenetic tree")
            if not self.tree_builder.run_tree_inference():
                self.logger.error("系统发育树构建失败|Tree construction failed")
                sys.exit(1)

            # 3. 一致性因子分析（可选）|Concordance factor analysis (optional)
            if self.config.enable_concordance:
                self.logger.info("步骤3: 一致性因子分析|Step 3: Concordance factor analysis")
                if not self.concordance_analyzer.run_concordance_analysis():
                    self.logger.warning("一致性因子分析失败，继续其他分析|Concordance analysis failed, continuing")

            # 4. 祖先状态重建（可选）|Ancestral state reconstruction (optional)
            if self.config.enable_ancestral:
                self.logger.info("步骤4: 祖先状态重建|Step 4: Ancestral state reconstruction")
                if not self.ancestral_reconstructor.run_ancestral_reconstruction():
                    self.logger.warning("祖先状态重建失败，继续其他分析|Ancestral reconstruction failed, continuing")
                else:
                    self.logger.info("步骤4.1: 解析祖先状态结果|Step 4.1: Parsing ancestral state results")
                    self.ancestral_parser.parse_ancestral_results()

            # 5. 完成总结|Final summary
            self.logger.info("=" * 60)
            self.logger.info("IQ-TREE 分析流程完成|IQ-TREE Analysis Pipeline Completed")
            self.logger.info("=" * 60)
            self.logger.info(f"结果保存在|Results saved in: {self.config.output_dir}")

        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止|Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='IQ-TREE 系统发育树分析脚本|IQ-TREE Phylogenetic Tree Analysis Script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i alignment.fasta -o tree_results -p my_tree
        """
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-i', '--input', required=True,
                         help='输入比对文件 (FASTA/PHYLIP/NEXUS)|Input alignment file (FASTA/PHYLIP/NEXUS)')
    required.add_argument('-o', '--output-dir', required=True,
                         help='输出目录|Output directory')
    required.add_argument('-p', '--prefix', required=True,
                         help='输出文件前缀|Output file prefix')

    # 核心参数|Core parameters
    core = parser.add_argument_group('核心参数|Core parameters')
    core.add_argument('-m', '--model',
                     help='进化模型 (不指定则自动选择)|Evolutionary model (auto-select if not specified)')
    core.add_argument('-t', '--threads', type=int, default=12,
                     help='线程数|Number of threads (default: 12)')

    # Bootstrap参数|Bootstrap parameters
    boot = parser.add_argument_group('Bootstrap参数|Bootstrap parameters')
    boot.add_argument('-b', '--bootstrap', type=int, default=1000,
                     help='Bootstrap重复次数|Bootstrap replicates (default: 1000)')
    boot.add_argument('--boot-type', choices=['ufboot', 'standard'], default='ufboot',
                     help='Bootstrap类型|Bootstrap type (default: ufboot)')
    boot.add_argument('--save-boot-trees', action='store_true',
                     help='保存所有bootstrap树到文件|Save all bootstrap trees to file')

    # 外群和约束|Outgroup and constraints
    tree_opts = parser.add_argument_group('树选项|Tree options')
    tree_opts.add_argument('--outgroup',
                          help='外群名称 (多个用逗号分隔)|Outgroup taxon names (comma-separated)')
    tree_opts.add_argument('--constraint',
                          help='约束树文件|Constraint tree file')

    # 分区分析|Partition analysis
    partition = parser.add_argument_group('分区分析|Partition analysis')
    partition.add_argument('--partition',
                          help='分区文件|Partition file')
    partition.add_argument('--partition-mode', choices=['edge-linked', 'edge-equal', 'edge-unlinked'],
                          default='edge-linked',
                          help='分区模式|Partition mode (default: edge-linked)')

    # 高级功能|Advanced features
    advanced = parser.add_argument_group('高级功能|Advanced features')
    advanced.add_argument('--concordance',
                         help='一致性因子分析：基因树文件|Concordance factor: gene tree file')
    advanced.add_argument('--ancestral', action='store_true',
                         help='启用祖先状态重建|Enable ancestral state reconstruction')

    # 其他参数|Other parameters
    other = parser.add_argument_group('其他参数|Other parameters')
    other.add_argument('--seed', type=int,
                      help='随机种子|Random seed')
    other.add_argument('--runs', type=int, default=1,
                      help='独立运行次数|Number of independent runs (default: 1)')
    other.add_argument('--redo', action='store_true',
                      help='重新运行分析|Redo analysis')
    other.add_argument('--iqtree-path',
                      help='IQ-TREE程序路径|IQ-TREE program path')

    args = parser.parse_args()

    # 构建配置参数|Build config parameters
    kwargs = dict(
        input_file=args.input,
        output_dir=args.output_dir,
        prefix=args.prefix,
        model=args.model,
        threads=args.threads,
        bootstrap=args.bootstrap,
        bootstrap_type=args.boot_type,
        save_boot_trees=args.save_boot_trees,
        outgroup=args.outgroup,
        partition_file=args.partition,
        partition_mode=args.partition_mode,
        enable_concordance=bool(args.concordance),
        concordance_trees=args.concordance,
        enable_ancestral=args.ancestral,
        constraint_tree=args.constraint,
        seed=args.seed,
        runs=args.runs,
        redo=args.redo,
    )

    if args.iqtree_path:
        kwargs['iqtree_path'] = expand_path(args.iqtree_path)

    analyzer = IQTreeAnalyzer(**kwargs)
    analyzer.run_analysis()


if __name__ == "__main__":
    main()
