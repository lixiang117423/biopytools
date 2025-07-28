# """
# 群体遗传分析主程序模块 | Population Genetics Analysis Main Module
# """

# import argparse
# import sys
# from .config import PopGenConfig
# from .utils import PopGenLogger, CommandRunner, check_dependencies
# from .preprocessing import VCFPreprocessor
# from .diversity import DiversityCalculator
# from .fst import FstCalculator
# from .ibd import IBDCalculator
# from .ld import LDCalculator
# from .ne import EffectivePopSizeCalculator
# from .results import ResultsProcessor, SummaryGenerator

# class PopulationGeneticsAnalyzer:
#     """群体遗传分析主类 | Main Population Genetics Analyzer Class"""
    
#     def __init__(self, **kwargs):
#         # 初始化配置 | Initialize configuration
#         self.config = PopGenConfig(**kwargs)
#         self.config.validate()
        
#         # 初始化日志 | Initialize logging
#         self.logger_manager = PopGenLogger(self.config.output_path)
#         self.logger = self.logger_manager.get_logger()
        
#         # 初始化命令执行器 | Initialize command runner
#         self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
#         # 初始化各个处理器 | Initialize processors
#         self.vcf_preprocessor = VCFPreprocessor(self.config, self.logger, self.cmd_runner)
#         self.diversity_calculator = DiversityCalculator(self.config, self.logger, self.cmd_runner)
#         self.fst_calculator = FstCalculator(self.config, self.logger, self.cmd_runner)
#         self.ibd_calculator = IBDCalculator(self.config, self.logger, self.cmd_runner)
#         self.ld_calculator = LDCalculator(self.config, self.logger, self.cmd_runner)
#         self.ne_calculator = EffectivePopSizeCalculator(self.config, self.logger, self.cmd_runner)
#         self.results_processor = ResultsProcessor(self.config, self.logger)
#         self.summary_generator = SummaryGenerator(self.config, self.logger)
    
#     def check_dependencies(self):
#         """检查依赖软件 | Check dependencies"""
#         return check_dependencies(self.config, self.logger)
    
#     def run_analysis(self):
#         """运行完整的群体遗传分析流程 | Run complete population genetics analysis pipeline"""
#         try:
#             self.logger.info("="*60)
#             self.logger.info("开始群体遗传分析流程 | Starting population genetics analysis pipeline")
#             self.logger.info("="*60)
            
#             # 检查依赖软件 | Check dependencies
#             if not self.check_dependencies():
#                 self.logger.error("依赖软件检查失败，退出分析 | Dependency check failed, exiting analysis")
#                 sys.exit(1)
            
#             # 步骤1: VCF预处理和质量控制 | Step 1: VCF preprocessing and quality control
#             self.logger.info("\n步骤1: VCF预处理和质量控制 | Step 1: VCF preprocessing and quality control")
#             filtered_vcf = self.vcf_preprocessor.quality_control()
#             if not filtered_vcf:
#                 self.logger.error("VCF质量控制失败 | VCF quality control failed")
#                 sys.exit(1)
            
#             # 步骤2: 加载分组信息 | Step 2: Load group information
#             self.logger.info("\n步骤2: 加载分组信息 | Step 2: Load group information")
#             group_dict = self.vcf_preprocessor.load_group_info()
            
#             # 步骤3: 计算多样性参数 | Step 3: Calculate diversity parameters
#             if any([self.config.calculate_pi, self.config.calculate_theta_w, self.config.calculate_tajima_d]):
#                 self.logger.info("\n步骤3: 计算多样性参数 | Step 3: Calculate diversity parameters")
#                 diversity_results = []
#                 for window_size in self.config.window_sizes:
#                     result = self.diversity_calculator.calculate_pi_theta_w_tajima_d(filtered_vcf, window_size)
#                     diversity_results.append(result)
#                 self.results_processor.process_diversity_results(diversity_results)
            
#             # 步骤4: 计算Fst | Step 4: Calculate Fst
#             if self.config.calculate_fst:
#                 self.logger.info("\n步骤4: 计算群体分化指数 (Fst) | Step 4: Calculate population differentiation (Fst)")
#                 fst_results = self.fst_calculator.calculate_fst(filtered_vcf, group_dict)
#                 self.results_processor.process_fst_results(fst_results)
            
#             # 步骤5: 计算IBD | Step 5: Calculate IBD
#             if self.config.calculate_ibd:
#                 self.logger.info("\n步骤5: 计算同源性分析 (IBD) | Step 5: Calculate identity-by-descent (IBD)")
#                 ibd_result = self.ibd_calculator.calculate_ibd(filtered_vcf)
            
#             # 步骤6: 计算LD | Step 6: Calculate LD
#             if self.config.calculate_ld:
#                 self.logger.info("\n步骤6: 计算连锁不平衡 (LD) | Step 6: Calculate linkage disequilibrium (LD)")
#                 ld_result = self.ld_calculator.calculate_ld(filtered_vcf)
            
#             # 步骤7: 计算有效群体大小 | Step 7: Calculate effective population size
#             if self.config.calculate_ne:
#                 self.logger.info("\n步骤7: 使用SMC++计算有效群体大小 | Step 7: Calculate effective population size using SMC++")
#                 ne_results = self.ne_calculator.calculate_ne_smcpp(filtered_vcf, group_dict)
            
#             # 步骤8: 生成总结报告 | Step 8: Generate summary report
#             self.logger.info("\n步骤8: 生成总结报告 | Step 8: Generate summary report")
#             self.summary_generator.generate_summary()
            
#             # 完成信息 | Completion information
#             self.logger.info("\n" + "="*60)
#             self.logger.info("群体遗传分析完成！| Population genetics analysis completed!")
#             self.logger.info("="*60)
#             self.logger.info(f"结果保存在 | Results saved in: {self.config.output_dir}")
            
#         except Exception as e:
#             self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
#             sys.exit(1)

# def main():
#     """主函数 | Main function"""
#     parser = argparse.ArgumentParser(
#         description='群体遗传分析工具 (模块化版本) | Population Genetics Analysis Tool (Modular Version)',
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#         epilog="""
# 示例 | Examples:
#   %(prog)s -v variants.vcf.gz -o popgen_results
#   %(prog)s -v data.vcf.gz -o results -g groups.txt --fst --pi --ld
#   %(prog)s -v variants.vcf -o results -g samples_groups.txt --all -t 8
#   %(prog)s -v filtered.vcf.gz -o results --ne --ibd -f csv
#         """
#     )
    
#     # 必需参数 | Required arguments
#     parser.add_argument('-v', '--vcf', required=True, 
#                        help='输入VCF文件路径 | Input VCF file path')
    
#     # 可选参数 | Optional arguments
#     parser.add_argument('-o', '--output', default='./popgen_output', 
#                        help='输出目录 | Output directory')
#     parser.add_argument('-g', '--groups', 
#                        help='分组信息文件 | Group information file')
    
#     # 分析选择参数 | Analysis selection parameters
#     parser.add_argument('--all', action='store_true', default=True,
#                        help='计算所有参数 (默认) | Calculate all parameters (default)')
#     parser.add_argument('--fst', action='store_true', 
#                        help='计算Fst | Calculate Fst')
#     parser.add_argument('--pi', action='store_true', 
#                        help='计算π | Calculate π')
#     parser.add_argument('--theta-w', action='store_true', 
#                        help='计算θw | Calculate θw')
#     parser.add_argument('--tajima-d', action='store_true', 
#                        help='计算Tajima\'s D | Calculate Tajima\'s D')
#     parser.add_argument('--ibd', action='store_true', 
#                        help='计算IBD | Calculate IBD')
#     parser.add_argument('--ld', action='store_true', 
#                        help='计算LD | Calculate LD')
#     parser.add_argument('--ne', action='store_true', 
#                        help='计算有效群体大小 | Calculate effective population size')
    
#     # 滑动窗口参数 | Sliding window parameters
#     parser.add_argument('-w', '--windows', nargs='+', type=int, default=[10000, 100000, 500000],
#                        help='滑动窗口大小 (bp) | Sliding window sizes (bp)')
#     parser.add_argument('--overlap', type=float, default=0.9, 
#                        help='窗口重叠率 | Window overlap rate')
    
#     # 质控参数 | Quality control parameters
#     parser.add_argument('-m', '--maf', type=float, default=0.01, 
#                        help='MAF阈值 | MAF threshold')
#     parser.add_argument('-M', '--missing', type=float, default=0.1, 
#                        help='缺失率阈值 | Missing rate threshold')
#     parser.add_argument('-H', '--hwe', type=float, default=1e-6, 
#                        help='HWE p值阈值 | HWE p-value threshold')
#     parser.add_argument('--min-dp', type=int, default=10, 
#                        help='最小测序深度 | Minimum depth')
#     parser.add_argument('--max-dp', type=int, default=100, 
#                        help='最大测序深度 | Maximum depth')
    
#     # 输出格式 | Output format
#     parser.add_argument('-f', '--format', choices=['txt', 'csv', 'tsv', 'json'], default='txt',
#                        help='输出文件格式 | Output file format')
    
#     # 计算资源 | Computing resources
#     parser.add_argument('-t', '--threads', type=int, default=4, 
#                        help='线程数 | Number of threads')
    
#     # 工具路径 | Tool paths
#     parser.add_argument('--vcftools-path', default='vcftools', 
#                        help='VCFtools路径 | VCFtools path')
#     parser.add_argument('--plink-path', default='plink', 
#                        help='PLINK路径 | PLINK path')
#     parser.add_argument('--bcftools-path', default='bcftools', 
#                        help='BCFtools路径 | BCFtools path')
#     parser.add_argument('--smcpp-path', default='smc++', 
#                        help='SMC++路径 | SMC++ path')
    
#     args = parser.parse_args()
    
#     # 设置分析参数 | Set analysis parameters
#     if not args.all:
#         specific_analyses = any([args.fst, args.pi, args.theta_w, args.tajima_d, 
#                                args.ibd, args.ld, args.ne])
#         if not specific_analyses:
#             args.all = True
    
#     # 设置分析参数 | Set analysis parameters
#     if args.all:
#         calculate_all = True
#         calculate_fst = True
#         calculate_pi = True
#         calculate_theta_w = True
#         calculate_tajima_d = True
#         calculate_ibd = True
#         calculate_ld = True
#         calculate_ne = True
#     else:
#         calculate_all = False
#         calculate_fst = args.fst
#         calculate_pi = args.pi
#         calculate_theta_w = args.theta_w
#         calculate_tajima_d = args.tajima_d
#         calculate_ibd = args.ibd
#         calculate_ld = args.ld
#         calculate_ne = args.ne
    
#     # 创建分析器并运行 | Create analyzer and run
#     analyzer = PopulationGeneticsAnalyzer(
#         vcf_file=args.vcf,
#         output_dir=args.output,
#         group_file=args.groups,
#         calculate_all=calculate_all,
#         calculate_fst=calculate_fst,
#         calculate_pi=calculate_pi,
#         calculate_theta_w=calculate_theta_w,
#         calculate_tajima_d=calculate_tajima_d,
#         calculate_ibd=calculate_ibd,
#         calculate_ld=calculate_ld,
#         calculate_ne=calculate_ne,
#         window_sizes=args.windows,
#         window_overlap=args.overlap,
#         maf=args.maf,
#         missing_rate=args.missing,
#         hwe_pvalue=args.hwe,
#         min_dp=args.min_dp,
#         max_dp=args.max_dp,
#         output_format=args.format,
#         threads=args.threads,
#         vcftools_path=args.vcftools_path,
#         plink_path=args.plink_path,
#         bcftools_path=args.bcftools_path,
#         smcpp_path=args.smcpp_path
#     )
    
#     analyzer.run_analysis()

"""
群体遗传分析主程序模块 | Population Genetics Analysis Main Module
"""

import argparse
import sys
from .config import PopGenConfig
from .utils import PopGenLogger, CommandRunner, check_dependencies
from .preprocessing import VCFPreprocessor
from .diversity import DiversityCalculator
from .fst import FstCalculator
from .ibd import IBDCalculator
from .ld import LDCalculator
from .ne import EffectivePopSizeCalculator
from .results import ResultsProcessor, SummaryGenerator

class PopulationGeneticsAnalyzer:
    """群体遗传分析主类 | Main Population Genetics Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = PopGenConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = PopGenLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.vcf_preprocessor = VCFPreprocessor(self.config, self.logger, self.cmd_runner)
        self.diversity_calculator = DiversityCalculator(self.config, self.logger, self.cmd_runner)
        self.fst_calculator = FstCalculator(self.config, self.logger, self.cmd_runner)
        self.ibd_calculator = IBDCalculator(self.config, self.logger, self.cmd_runner)
        self.ld_calculator = LDCalculator(self.config, self.logger, self.cmd_runner)
        self.ne_calculator = EffectivePopSizeCalculator(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的群体遗传分析流程 | Run complete population genetics analysis pipeline"""
        try:
            self.logger.info("="*60)
            self.logger.info("开始群体遗传分析流程 | Starting population genetics analysis pipeline")
            self.logger.info("="*60)
            
            # 检查依赖软件 | Check dependencies
            if not self.check_dependencies():
                self.logger.error("依赖软件检查失败，退出分析 | Dependency check failed, exiting analysis")
                sys.exit(1)
            
            # 步骤1: VCF预处理和质量控制 | Step 1: VCF preprocessing and quality control
            self.logger.info("\n步骤1: VCF预处理和质量控制 | Step 1: VCF preprocessing and quality control")
            
            # 检查VCF文件格式 | Check VCF file format
            vcf_info = self.vcf_preprocessor.check_vcf_format(self.config.vcf_file)
            
            # 执行质量控制（或跳过） | Execute quality control (or skip)
            filtered_vcf = self.vcf_preprocessor.quality_control()
            if not filtered_vcf:
                self.logger.error("VCF质量控制失败 | VCF quality control failed")
                sys.exit(1)
            
            # 步骤2: 加载分组信息 | Step 2: Load group information
            self.logger.info("\n步骤2: 加载分组信息 | Step 2: Load group information")
            group_dict = self.vcf_preprocessor.load_group_info()
            
            # 步骤3: 计算多样性参数 | Step 3: Calculate diversity parameters
            if any([self.config.calculate_pi, self.config.calculate_theta_w, self.config.calculate_tajima_d]):
                self.logger.info("\n步骤3: 计算多样性参数 | Step 3: Calculate diversity parameters")
                diversity_results = []
                for window_size in self.config.window_sizes:
                    result = self.diversity_calculator.calculate_pi_theta_w_tajima_d(filtered_vcf, window_size)
                    diversity_results.append(result)
            
            # 步骤4: 计算Fst | Step 4: Calculate Fst
            if self.config.calculate_fst and group_dict:
                self.logger.info("\n步骤4: 计算Fst | Step 4: Calculate Fst")
                fst_results = []
                for window_size in self.config.window_sizes:
                    result = self.fst_calculator.calculate_fst(filtered_vcf, group_dict, window_size)
                    fst_results.append(result)
            elif self.config.calculate_fst and not group_dict:
                self.logger.warning("跳过Fst计算：未提供分组信息 | Skipping Fst calculation: no group information provided")
            
            # 步骤5: 计算IBD | Step 5: Calculate IBD
            if self.config.calculate_ibd:
                self.logger.info("\n步骤5: 计算IBD | Step 5: Calculate IBD")
                ibd_result = self.ibd_calculator.calculate_ibd(filtered_vcf)
            
            # 步骤6: 计算LD | Step 6: Calculate LD
            if self.config.calculate_ld:
                self.logger.info("\n步骤6: 计算LD | Step 6: Calculate LD")
                ld_result = self.ld_calculator.calculate_ld(filtered_vcf)
            
            # 步骤7: 计算有效群体大小 | Step 7: Calculate effective population size
            if self.config.calculate_ne and group_dict:
                self.logger.info("\n步骤7: 计算有效群体大小 | Step 7: Calculate effective population size")
                ne_results = []
                for group in set(group_dict.values()):
                    result = self.ne_calculator.calculate_effective_population_size(filtered_vcf, group, group_dict)
                    ne_results.append(result)
            elif self.config.calculate_ne and not group_dict:
                self.logger.warning("跳过Ne计算：未提供分组信息 | Skipping Ne calculation: no group information provided")
            
            # 步骤8: 处理和总结结果 | Step 8: Process and summarize results
            self.logger.info("\n步骤8: 处理和总结结果 | Step 8: Process and summarize results")
            
            # 处理所有结果 | Process all results
            self.results_processor.process_all_results()
            
            # 生成总结报告 | Generate summary report
            summary_file = self.summary_generator.generate_summary()
            
            # 完成信息 | Completion information
            self.logger.info("\n" + "="*60)
            self.logger.info("群体遗传分析完成！| Population genetics analysis completed!")
            self.logger.info(f"结果保存在 | Results saved in: {self.config.output_dir}")
            self.logger.info(f"总结报告 | Summary report: {summary_file}")
            self.logger.info("="*60)
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='群体遗传分析工具包 (模块化版本) | Population Genetics Analysis Toolkit (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s -v variants.vcf.gz -o popgen_results --pi
  %(prog)s -v data.vcf.gz -o results -g groups.txt --all
  %(prog)s -v filtered.vcf.gz -o output --skip-qc --fst --pi
  %(prog)s -v data.vcf.gz -o results -w 10000,100000 -t 8
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-v', '--vcf', required=True,
                       help='输入VCF文件路径 (支持压缩和未压缩) | Input VCF file path (supports compressed and uncompressed)')
    
    # 基本参数 | Basic parameters
    parser.add_argument('-o', '--output', default='./popgen_output',
                       help='输出目录 | Output directory')
    parser.add_argument('-g', '--groups',
                       help='分组信息文件 (制表符分隔: 样本名\t分组名) | Group information file (tab-separated: sample_name\\tgroup_name)')
    
    # 分析参数选择 | Analysis parameter selection
    analysis_group = parser.add_argument_group('分析参数选择 | Analysis Parameter Selection')
    analysis_group.add_argument('--all', action='store_true',
                               help='计算所有群体遗传参数 | Calculate all population genetics parameters')
    analysis_group.add_argument('--fst', action='store_true',
                               help='计算Fst | Calculate Fst')
    analysis_group.add_argument('--pi', action='store_true',
                               help='计算核苷酸多样性π | Calculate nucleotide diversity π')
    analysis_group.add_argument('--theta-w', action='store_true',
                               help='计算Watterson\'s θw | Calculate Watterson\'s θw')
    analysis_group.add_argument('--tajima-d', action='store_true',
                               help='计算Tajima\'s D | Calculate Tajima\'s D')
    analysis_group.add_argument('--ibd', action='store_true',
                               help='计算同源性分析 | Calculate IBD analysis')
    analysis_group.add_argument('--ld', action='store_true',
                               help='计算连锁不平衡 | Calculate linkage disequilibrium')
    analysis_group.add_argument('--ne', action='store_true',
                               help='计算有效群体大小 | Calculate effective population size')
    
    # 滑动窗口参数 | Sliding window parameters
    window_group = parser.add_argument_group('滑动窗口参数 | Sliding Window Parameters')
    window_group.add_argument('-w', '--windows',
                             help='窗口大小列表，逗号分隔 (如: 10000,100000,500000) | Window sizes, comma-separated (e.g., 10000,100000,500000)')
    window_group.add_argument('--window-overlap', type=float, default=0.9,
                             help='窗口重叠比例 | Window overlap ratio')
    
    # 质控参数 | Quality control parameters
    qc_group = parser.add_argument_group('质量控制参数 | Quality Control Parameters')
    qc_group.add_argument('-m', '--maf', type=float, default=0.01,
                         help='最小等位基因频率阈值 | Minor allele frequency threshold')
    qc_group.add_argument('-M', '--missing', type=float, default=0.1,
                         help='最大缺失率阈值 | Maximum missing rate threshold')
    qc_group.add_argument('-H', '--hwe', type=float, default=1e-6,
                         help='Hardy-Weinberg平衡p值阈值 | Hardy-Weinberg equilibrium p-value threshold')
    qc_group.add_argument('--min-dp', type=int, default=10,
                         help='最小测序深度 | Minimum sequencing depth')
    qc_group.add_argument('--max-dp', type=int, default=100,
                         help='最大测序深度 | Maximum sequencing depth')
    qc_group.add_argument('--skip-qc', action='store_true',
                         help='跳过VCF质量控制过滤，直接使用输入文件 | Skip VCF quality control filtering, use input file directly')
    
    # 输出参数 | Output parameters
    output_group = parser.add_argument_group('输出参数 | Output Parameters')
    output_group.add_argument('-f', '--format', choices=['txt', 'csv', 'json'], default='txt',
                             help='输出格式 | Output format')
    
    # 软件路径 | Software paths
    software_group = parser.add_argument_group('软件路径 | Software Paths')
    software_group.add_argument('--vcftools-path', default='vcftools',
                               help='VCFtools软件路径 | VCFtools software path')
    software_group.add_argument('--plink-path', default='plink',
                               help='PLINK软件路径 | PLINK software path')
    software_group.add_argument('--bcftools-path', default='bcftools',
                               help='BCFtools软件路径 | BCFtools software path')
    software_group.add_argument('--smcpp-path', default='smc++',
                               help='SMC++软件路径 | SMC++ software path')
    
    # 系统参数 | System parameters
    system_group = parser.add_argument_group('系统参数 | System Parameters')
    system_group.add_argument('-t', '--threads', type=int, default=4,
                             help='线程数 | Number of threads')
    
    args = parser.parse_args()
    
    # 处理窗口大小参数 | Process window sizes parameter
    window_sizes = None
    if args.windows:
        try:
            window_sizes = [int(w.strip()) for w in args.windows.split(',')]
        except ValueError:
            parser.error("窗口大小必须是逗号分隔的整数 | Window sizes must be comma-separated integers")
    
    # 处理分析参数选择 | Process analysis parameter selection
    if args.all:
        calculate_fst = calculate_pi = calculate_theta_w = True
        calculate_tajima_d = calculate_ibd = calculate_ld = calculate_ne = True
    else:
        calculate_fst = args.fst
        calculate_pi = args.pi
        calculate_theta_w = args.theta_w
        calculate_tajima_d = args.tajima_d
        calculate_ibd = args.ibd
        calculate_ld = args.ld
        calculate_ne = args.ne
        
        # 如果没有选择任何分析，默认计算pi | If no analysis selected, default to pi
        if not any([calculate_fst, calculate_pi, calculate_theta_w, calculate_tajima_d, 
                   calculate_ibd, calculate_ld, calculate_ne]):
            calculate_pi = True
            print("未指定分析参数，默认计算核苷酸多样性π | No analysis specified, defaulting to nucleotide diversity π")
    
    # 创建分析器并运行 | Create analyzer and run
    try:
        analyzer = PopulationGeneticsAnalyzer(
            vcf_file=args.vcf,
            output_dir=args.output,
            group_file=args.groups,
            window_sizes=window_sizes,
            window_overlap=args.window_overlap,
            calculate_fst=calculate_fst,
            calculate_pi=calculate_pi,
            calculate_theta_w=calculate_theta_w,
            calculate_tajima_d=calculate_tajima_d,
            calculate_ibd=calculate_ibd,
            calculate_ld=calculate_ld,
            calculate_ne=calculate_ne,
            maf=args.maf,
            missing_rate=args.missing,
            hwe_pvalue=args.hwe,
            min_dp=args.min_dp,
            max_dp=args.max_dp,
            skip_qc=args.skip_qc,  # 新增参数 | New parameter
            output_format=args.format,
            vcftools_path=args.vcftools_path,
            plink_path=args.plink_path,
            bcftools_path=args.bcftools_path,
            smcpp_path=args.smcpp_path,
            threads=args.threads
        )
        
        analyzer.run_analysis()
        
    except Exception as e:
        print(f"错误 | Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()