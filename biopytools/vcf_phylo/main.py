"""
VCF系统发育分析主程序模块 | VCF Phylogenetic Analysis Main Module
"""

import argparse
import sys
from .config import PhyloConfig
from .utils import PhyloLogger, CommandRunner, check_dependencies
from .distance_calc import DistanceCalculator
from .tree_builder import TreeBuilder
from .results import ResultsManager

class VCFPhyloAnalyzer:
    """VCF系统发育分析主类 | Main VCF Phylogenetic Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = PhyloConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = PhyloLogger(self.config.output_prefix)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.working_path)
        
        # 初始化各个处理器 | Initialize processors
        self.distance_calculator = DistanceCalculator(self.config, self.logger, self.cmd_runner)
        self.tree_builder = TreeBuilder(self.config, self.logger)
        self.results_manager = ResultsManager(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的系统发育分析流程 | Run complete phylogenetic analysis pipeline"""
        try:
            self.logger.info("="*80)
            self.logger.info("开始VCF系统发育分析 | Starting VCF phylogenetic analysis")
            self.logger.info("="*80)
            
            # 步骤0: 检查依赖 | Step 0: Check dependencies
            if not self.check_dependencies():
                raise RuntimeError("依赖检查失败 | Dependency check failed")
            
            # 步骤1: 计算距离矩阵 | Step 1: Calculate distance matrix
            self.logger.info("-" * 40)
            self.logger.info("步骤1: 计算距离矩阵 | Step 1: Calculate distance matrix")
            self.logger.info("-" * 40)
            
            if not self.distance_calculator.calculate_distance_matrix():
                raise RuntimeError("距离矩阵计算失败 | Distance matrix calculation failed")
            
            # 步骤2: 构建NJ系统发育树 | Step 2: Build NJ phylogenetic tree
            self.logger.info("-" * 40)
            self.logger.info("步骤2: 构建系统发育树 | Step 2: Build phylogenetic tree")
            self.logger.info("-" * 40)
            
            if not self.tree_builder.build_nj_tree():
                raise RuntimeError("NJ系统发育树构建失败 | NJ phylogenetic tree construction failed")
            
            # 步骤3: 验证结果和生成报告 | Step 3: Validate results and generate report
            self.logger.info("-" * 40)
            self.logger.info("步骤3: 生成结果报告 | Step 3: Generate results report")
            self.logger.info("-" * 40)
            
            if not self.results_manager.validate_results():
                self.logger.warning("结果验证发现问题 | Issues found during result validation")
            
            # 生成总结报告 | Generate summary report
            self.results_manager.generate_summary_report()
            
            self.logger.info("="*80)
            self.logger.info("VCF系统发育分析完成 | VCF phylogenetic analysis completed successfully")
            self.logger.info("="*80)
            
            # 输出结果文件信息 | Output result file information
            self.results_manager.log_output_files()
            
        except Exception as e:
            self.logger.error(f"分析失败 | Analysis failed: {e}")
            raise

def create_parser():
    """创建命令行参数解析器 | Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description="VCF系统发育分析工具 v2.0 | VCF Phylogenetic Analysis Tool v2.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  # 基本分析 | Basic analysis
  %(prog)s -i wild.snp.vcf -o wild_snp
  
  # 指定完整路径 | Specify full paths
  %(prog)s --input 01.data/wild.snp.new.record.vcf.recode.vcf \\
           --output 02.tree/wild_snp_dis_mat \\
           --tree-output 02.tree/wild_snp_dis.nwk
  
  # 从已有距离矩阵构建树 | Build tree from existing matrix
  %(prog)s --distance-matrix existing_matrix.txt \\
           --tree-output tree.nwk \\
           --skip-vcf2dis
        """
    )
    
    # 输入文件参数 | Input file arguments
    input_group = parser.add_argument_group("输入文件 | Input Files")
    input_group.add_argument(
        "-i", "--input",
        dest="vcf_file",
        help="输入VCF文件路径 | Input VCF file path"
    )
    input_group.add_argument(
        "-d", "--distance-matrix",
        help="已有的距离矩阵文件路径（用于跳过VCF2Dis步骤）| Existing distance matrix file path (for skipping VCF2Dis step)"
    )
    
    # 输出文件参数 | Output file arguments
    output_group = parser.add_argument_group("输出文件 | Output Files")
    output_group.add_argument(
        "-o", "--output",
        dest="output_prefix",
        default="phylo_analysis",
        help="输出文件前缀（默认: phylo_analysis）| Output file prefix (default: phylo_analysis)"
    )
    output_group.add_argument(
        "-t", "--tree-output",
        help="系统发育树输出文件路径（默认: OUTPUT.nwk）| Phylogenetic tree output file path (default: OUTPUT.nwk)"
    )
    
    # 工具参数 | Tool arguments
    tool_group = parser.add_argument_group("工具设置 | Tool Settings")
    tool_group.add_argument(
        "--vcf2dis-path",
        default="VCF2Dis",
        help="VCF2Dis程序路径（默认: VCF2Dis）| VCF2Dis program path (default: VCF2Dis)"
    )
    tool_group.add_argument(
        "-w", "--working-dir",
        default=".",
        help="工作目录（默认: 当前目录）| Working directory (default: current directory)"
    )
    
    # 行为控制参数 | Behavior control arguments
    behavior_group = parser.add_argument_group("行为控制 | Behavior Control")
    behavior_group.add_argument(
        "--skip-vcf2dis",
        action="store_true",
        help="跳过VCF2Dis步骤，直接从距离矩阵构建树 | Skip VCF2Dis step, build tree directly from distance matrix"
    )
    
    return parser

def main():
    """主函数 | Main function"""
    parser = create_parser()
    args = parser.parse_args()
    
    # 检查参数逻辑 | Check argument logic
    if args.skip_vcf2dis:
        if not args.distance_matrix:
            parser.error("使用 --skip-vcf2dis 时必须指定 --distance-matrix | --distance-matrix must be specified when using --skip-vcf2dis")
    else:
        if not args.vcf_file:
            parser.error("必须指定输入VCF文件 | Input VCF file must be specified")
    
    try:
        # 创建分析器 | Create analyzer
        analyzer = VCFPhyloAnalyzer(
            vcf_file=args.vcf_file,
            output_prefix=args.output_prefix,
            distance_matrix=args.distance_matrix,
            tree_output=args.tree_output,
            vcf2dis_path=args.vcf2dis_path,
            skip_vcf2dis=args.skip_vcf2dis,
            working_dir=args.working_dir
        )
        
        # 运行分析 | Run analysis
        analyzer.run_analysis()
        
    except Exception as e:
        print(f"错误 | Error: {e}", file=sys.stderr)
        sys.exit(1)
