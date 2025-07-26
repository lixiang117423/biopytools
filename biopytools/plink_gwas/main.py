# """
# PLINK GWAS分析主程序模块 | PLINK GWAS Analysis Main Module
# """

# import argparse
# import sys
# import os
# from .config import PlinkGWASConfig
# from .utils import PlinkGWASLogger, CommandRunner, FileManager
# from .data_processing import PhenotypeProcessor, VCFProcessor, QualityController
# from .population_analysis import PopulationAnalyzer
# from .association_analysis import AssociationAnalyzer
# from .results import ResultsProcessor, ReportGenerator, VisualizationGenerator

# class PlinkGWAS:
#     """PLINK GWAS分析主类 | Main PLINK GWAS Analysis Class"""
    
#     def __init__(self, **kwargs):
#         # 初始化配置 | Initialize configuration
#         self.config = PlinkGWASConfig(**kwargs)
#         self.config.validate()
        
#         # 切换到工作目录 | Change to working directory
#         os.chdir(self.config.output_path)
        
#         # 初始化日志 | Initialize logging
#         self.logger_manager = PlinkGWASLogger(self.config.output_path)
#         self.logger = self.logger_manager.get_logger()
        
#         # 初始化命令执行器 | Initialize command runner
#         self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
#         # 初始化文件管理器 | Initialize file manager
#         self.file_manager = FileManager(self.logger, self.config.output_path)
        
#         # 初始化各个处理器 | Initialize processors
#         self.pheno_processor = PhenotypeProcessor(self.config, self.logger, self.cmd_runner)
#         self.vcf_processor = VCFProcessor(self.config, self.logger, self.cmd_runner)
#         self.quality_controller = QualityController(self.config, self.logger, self.cmd_runner)
#         self.population_analyzer = PopulationAnalyzer(self.config, self.logger, self.cmd_runner)
#         self.association_analyzer = AssociationAnalyzer(self.config, self.logger, self.cmd_runner)
#         self.results_processor = ResultsProcessor(self.config, self.logger)
#         self.report_generator = ReportGenerator(self.config, self.logger)
#         self.viz_generator = VisualizationGenerator(self.config, self.logger, self.cmd_runner)
    
#     def run_analysis(self):
#         """运行完整的PLINK GWAS分析流程 | Run complete PLINK GWAS analysis pipeline"""
#         try:
#             self.logger.info("开始PLINK GWAS分析流程 | Starting PLINK GWAS analysis pipeline...")
#             self.logger.info(f"表型类型 | Trait type: {self.config.trait_type}")
#             self.logger.info(f"显著性校正方法 | Correction method: {self.config.correction_method}")
#             self.logger.info("="*50)
            
#             # 1. 检查并复制输入文件 | Check and copy input files
#             self.logger.info("步骤1: 检查输入文件 | Step 1: Checking input files")
#             self.file_manager.copy_input_files(self.config.vcf_file, self.config.phenotype_file)
            
#             # 2. 转换表型文件 | Convert phenotype file
#             self.logger.info("步骤2: 转换表型文件 | Step 2: Converting phenotype file")
#             self.pheno_processor.convert_phenotype()
            
#             # 3. 转换VCF文件 | Convert VCF file
#             self.logger.info("步骤3: 转换VCF文件 | Step 3: Converting VCF file")
#             self.vcf_processor.convert_vcf_to_plink()
            
#             # 4. 合并表型信息 | Merge phenotype information
#             self.logger.info("步骤4: 合并表型信息 | Step 4: Merging phenotype information")
#             self.vcf_processor.merge_phenotype()
            
#             # 5. 质量控制 | Quality control
#             self.logger.info("步骤5: 质量控制 | Step 5: Quality control")
#             self.quality_controller.quality_control()
            
#             # 6. 群体分层控制 | Population stratification control
#             self.logger.info("步骤6: 群体分层控制 | Step 6: Population stratification control")
#             use_pca = self.population_analyzer.population_stratification()
            
#             # 7. 关联分析 | Association analysis
#             self.logger.info("步骤7: 关联分析 | Step 7: Association analysis")
#             main_result = self.association_analyzer.association_analysis(use_pca)
            
#             # 8. 处理结果 | Process results
#             self.logger.info("步骤8: 处理结果 | Step 8: Processing results")
#             results_df = self.results_processor.process_results(main_result)
            
#             # 9. 生成报告 | Generate report
#             self.logger.info("步骤9: 生成报告 | Step 9: Generating report")
#             self.report_generator.generate_report(results_df, main_result)
            
#             # 10. 生成可视化 | Generate visualization
#             self.logger.info("步骤10: 生成可视化 | Step 10: Generating visualization")
#             self.viz_generator.generate_plots()
            
#             self.logger.info("="*50)
#             self.logger.info("PLINK GWAS分析流程完成！ | PLINK GWAS analysis pipeline completed!")
#             self.logger.info(f"表型类型 | Trait type: {self.config.trait_type}")
#             self.logger.info(f"显著性校正方法 | Correction method: {self.config.correction_method}")
#             self.logger.info(f"所有结果保存在 | All results saved in: {self.config.output_path.absolute()}")
#             self.logger.info("主要输出文件 | Main output files:")
#             self.logger.info("- analysis_report.txt: 分析报告 | Analysis report")
#             self.logger.info("- gwas_results_ADD.txt: 主要关联结果 | Main association results")
            
#             # 根据校正方法显示相应的输出文件 | Display output files based on correction method
#             if self.config.correction_method in ["bonferroni", "all"]:
#                 self.logger.info("- significant_bonferroni.txt: Bonferroni校正显著位点 | Bonferroni significant loci")
#                 self.logger.info("- bonferroni_info.txt: Bonferroni校正信息 | Bonferroni correction information")
            
#             if self.config.correction_method in ["suggestive", "all"]:
#                 self.logger.info("- significant_suggestive.txt: 提示性关联位点 | Suggestive association loci")
#                 self.logger.info("- suggestive_info.txt: 提示性关联信息 | Suggestive association information")
            
#             if self.config.correction_method in ["fdr", "all"]:
#                 self.logger.info("- significant_fdr.txt: FDR校正显著位点 | FDR significant loci")
#                 self.logger.info("- fdr_info.txt: FDR校正信息 | FDR correction information")
            
#             self.logger.info("- manhattan_plot.png: Manhattan图 | Manhattan plot")
#             self.logger.info("- qq_plot.png: QQ图 | QQ plot")
            
#         except Exception as e:
#             self.logger.error(f"分析流程失败 | Analysis pipeline failed: {e}")
#             sys.exit(1)

# def main():
#     """主函数 | Main function"""
#     parser = argparse.ArgumentParser(
#         description="完整的PLINK GWAS分析流程 (模块化版本) - 支持质量性状和数量性状，多种显著性校正方法 | Complete PLINK GWAS Analysis Pipeline (Modular Version) - Supporting both qualitative and quantitative traits, multiple significance correction methods",
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter
#     )
    
#     # 必需参数 | Required parameters
#     parser.add_argument("-v", "--vcf-file", required=True,
#                        help="输入VCF文件路径（支持.gz压缩） | Input VCF file path (supports .gz compression)")
#     parser.add_argument("-p", "--phenotype-file", required=True,
#                        help="表型文件路径（样本ID和表型值，以空格或制表符分隔） | Phenotype file path (sample ID and phenotype value, space or tab separated)")
    
#     # 表型类型参数 | Trait type parameter
#     parser.add_argument("-t", "--trait-type", choices=["qualitative", "quantitative"], default="qualitative",
#                        help="表型类型 | Trait type: 'qualitative' for binary traits (0/1 -> 1/2), 'quantitative' for continuous traits (keep original values)")
    
#     # 输出参数 | Output parameters
#     parser.add_argument("-o", "--output-dir", default="plink_results",
#                        help="输出目录 | Output directory")
    
#     # 质控参数 | Quality control parameters
#     parser.add_argument("--mind", type=float, default=0.05,
#                        help="个体缺失率阈值（移除缺失率大于此值的个体） | Individual missing rate threshold")
#     parser.add_argument("--geno", type=float, default=0.05,
#                        help="SNP缺失率阈值（移除缺失率大于此值的SNP） | SNP missing rate threshold")
#     parser.add_argument("--maf", type=float, default=0.01,
#                        help="最小等位基因频率阈值 | Minor allele frequency threshold")
#     parser.add_argument("--hwe", type=float, default=1e-6,
#                        help="Hardy-Weinberg平衡检验P值阈值 | Hardy-Weinberg equilibrium p-value threshold")
    
#     # LD剪枝参数 | LD pruning parameters
#     parser.add_argument("--ld-window-size", type=int, default=50,
#                        help="LD剪枝窗口大小（kb） | LD pruning window size (kb)")
#     parser.add_argument("--ld-step-size", type=int, default=5,
#                        help="LD剪枝步长（SNP数） | LD pruning step size (number of SNPs)")
#     parser.add_argument("--ld-r2-threshold", type=float, default=0.2,
#                        help="LD剪枝r²阈值 | LD pruning r² threshold")
    
#     # 主成分参数 | PCA parameters
#     parser.add_argument("--pca-components", type=int, default=10,
#                        help="计算的主成分数量 | Number of principal components to compute")
#     parser.add_argument("--pca-use", type=int, default=5,
#                        help="关联分析中使用的主成分数量 | Number of PCs to use in association analysis")
    
#     # 显著性校正参数 | Significance correction parameters
#     parser.add_argument("--correction-method", choices=["bonferroni", "suggestive", "fdr", "all"], default="all",
#                        help="显著性校正方法 | Significance correction method: 'bonferroni' for Bonferroni correction, 'suggestive' for suggestive threshold, 'fdr' for false discovery rate, 'all' for all methods")
#     parser.add_argument("--bonferroni-alpha", type=float, default=0.05,
#                        help="Bonferroni校正的alpha水平 | Alpha level for Bonferroni correction")
#     parser.add_argument("--suggestive-threshold", type=float, default=1e-5,
#                        help="提示性关联阈值 | Suggestive association threshold")
#     parser.add_argument("--fdr-alpha", type=float, default=0.05,
#                        help="FDR校正的q值阈值 | q-value threshold for FDR correction")
    
#     # 其他参数 | Other parameters
#     parser.add_argument("--threads", type=int, default=1,
#                        help="使用的线程数 | Number of threads to use")
    
#     args = parser.parse_args()
    
#     # 创建分析器并运行 | Create analyzer and run
#     plink = PlinkGWAS(**vars(args))
#     plink.run_analysis()

# if __name__ == "__main__":
#     main()


# 20250726添加显性和隐性表型类型支持 | Added support for dominant and recessive trait types on 20250726
"""
PLINK GWAS分析主程序模块 | PLINK GWAS Analysis Main Module
"""

import argparse
import sys
import os
from pathlib import Path
from .config import PlinkGWASConfig
from .utils import PlinkGWASLogger, CommandRunner, FileManager
from .data_processing import PhenotypeProcessor, VCFProcessor, QualityController
from .population_analysis import PopulationAnalyzer
from .association_analysis import AssociationAnalyzer
from .results import ResultsProcessor, ReportGenerator, VisualizationGenerator

class PlinkGWAS:
    """PLINK GWAS分析主类 | Main PLINK GWAS Analysis Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = PlinkGWASConfig(**kwargs)
        self.config.validate()
        
        # 确保输出目录存在 | Ensure output directory exists
        self.config.output_path.mkdir(parents=True, exist_ok=True)
        
        # 保存原始工作目录 | Save original working directory
        self.original_cwd = Path.cwd()
        
        # 初始化日志（在切换目录之前）| Initialize logging (before changing directory)
        self.logger_manager = PlinkGWASLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 切换到工作目录 | Change to working directory
        os.chdir(self.config.output_path)
        self.logger.info(f"工作目录 | Working directory: {self.config.output_path}")
        
        # 初始化命令执行器和文件管理器（使用当前目录"."）| Initialize command runner and file manager (use current directory ".")
        current_dir = Path(".")
        self.cmd_runner = CommandRunner(self.logger, current_dir)
        self.file_manager = FileManager(self.logger, current_dir)
        
        # 初始化各个处理器 | Initialize processors
        self.pheno_processor = PhenotypeProcessor(self.config, self.logger, self.cmd_runner)
        self.vcf_processor = VCFProcessor(self.config, self.logger, self.cmd_runner)
        self.quality_controller = QualityController(self.config, self.logger, self.cmd_runner)
        self.population_analyzer = PopulationAnalyzer(self.config, self.logger, self.cmd_runner)
        self.association_analyzer = AssociationAnalyzer(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.report_generator = ReportGenerator(self.config, self.logger)
        self.viz_generator = VisualizationGenerator(self.config, self.logger, self.cmd_runner)
    
    def __del__(self):
        """析构函数，恢复原始工作目录 | Destructor, restore original working directory"""
        try:
            if hasattr(self, 'original_cwd'):
                os.chdir(self.original_cwd)
        except:
            pass
    
    def run_analysis(self):
        """运行完整的PLINK GWAS分析流程 | Run complete PLINK GWAS analysis pipeline"""
        try:
            self.logger.info("开始PLINK GWAS分析流程 | Starting PLINK GWAS analysis pipeline...")
            self.logger.info(f"表型类型 | Trait type: {self.config.trait_type}")
            self.logger.info(f"遗传模型 | Genetic model: {self.config.genetic_model}")
            self.logger.info(f"显著性校正方法 | Correction method: {self.config.correction_method}")
            self.logger.info("="*50)
            
            # 1. 检查并复制输入文件 | Check and copy input files
            self.logger.info("步骤1: 检查输入文件 | Step 1: Checking input files")
            # 构建绝对路径 | Build absolute paths
            vcf_path = self.original_cwd / self.config.vcf_file if not Path(self.config.vcf_file).is_absolute() else Path(self.config.vcf_file)
            pheno_path = self.original_cwd / self.config.phenotype_file if not Path(self.config.phenotype_file).is_absolute() else Path(self.config.phenotype_file)
            
            self.file_manager.copy_input_files(str(vcf_path), str(pheno_path))
            
            # 2. 转换表型文件 | Convert phenotype file
            self.logger.info("步骤2: 转换表型文件 | Step 2: Converting phenotype file")
            self.pheno_processor.convert_phenotype()
            
            # 3. 转换VCF文件 | Convert VCF file
            self.logger.info("步骤3: 转换VCF文件 | Step 3: Converting VCF file")
            self.vcf_processor.convert_vcf_to_plink()
            
            # 4. 合并表型信息 | Merge phenotype information
            self.logger.info("步骤4: 合并表型信息 | Step 4: Merging phenotype information")
            self.vcf_processor.merge_phenotype()
            
            # 5. 数据质量控制 | Data quality control
            self.logger.info("步骤5: 数据质量控制 | Step 5: Data quality control")
            self.quality_controller.quality_control()
            
            # 6. 群体结构分析 | Population structure analysis
            self.logger.info("步骤6: 群体结构分析 | Step 6: Population structure analysis")
            self.population_analyzer.population_analysis()
            
            # 7. 关联分析 | Association analysis
            self.logger.info("步骤7: 关联分析 | Step 7: Association analysis")
            main_result = self.association_analyzer.association_analysis(use_pca=True)
            
            # 7. 结果处理 | Results processing
            self.logger.info("步骤7: 结果处理 | Step 7: Results processing")
            model_results = self.results_processor.process_results(main_result)
            
            # 8. 生成报告 | Generate reports
            self.logger.info("步骤8: 生成报告 | Step 8: Generating reports")
            stats = self._collect_statistics()
            self.report_generator.generate_summary_report(model_results, stats)
            
            # 9. 生成可视化 | Generate visualizations
            self.logger.info("步骤9: 生成可视化 | Step 9: Generating visualizations")
            self.viz_generator.generate_plots(model_results)
            
            self.logger.info("="*50)
            self.logger.info("PLINK GWAS分析完成! | PLINK GWAS analysis completed!")
            
            # 输出结果摘要 | Output results summary
            self._print_results_summary(model_results, stats)
            
        except Exception as e:
            self.logger.error(f"分析过程中发生错误 | Error during analysis: {e}")
            raise
        finally:
            # 确保最后恢复原始工作目录 | Ensure original working directory is restored
            try:
                os.chdir(self.original_cwd)
            except:
                pass
    
    def _collect_statistics(self) -> dict:
        """收集分析统计信息 | Collect analysis statistics"""
        stats = {}
        
        try:
            # 读取样本信息 | Read sample information
            if os.path.exists("data_qc1.fam"):
                with open("data_qc1.fam", 'r') as f:
                    lines = f.readlines()
                    stats['total_samples'] = len(lines)
                    
                    # 如果是质量性状，统计病例对照数 | Count cases/controls for qualitative traits
                    if self.config.trait_type == "qualitative":
                        cases = sum(1 for line in lines if line.strip().split()[-1] == '2')
                        controls = sum(1 for line in lines if line.strip().split()[-1] == '1')
                        stats['cases'] = cases
                        stats['controls'] = controls
            
            # 读取SNP信息 | Read SNP information
            if os.path.exists("data_qc1.bim"):
                with open("data_qc1.bim", 'r') as f:
                    lines = f.readlines()
                    stats['total_snps'] = len(lines)
                    
                    # 统计染色体数 | Count chromosomes
                    chromosomes = set(line.strip().split()[0] for line in lines)
                    stats['chromosomes'] = len(chromosomes)
            
        except Exception as e:
            self.logger.warning(f"收集统计信息时出错 | Error collecting statistics: {e}")
        
        return stats
    
    def _print_results_summary(self, model_results: dict, stats: dict):
        """打印结果摘要 | Print results summary"""
        print("\n" + "="*60)
        print("PLINK GWAS 分析结果摘要 | Analysis Results Summary")
        print("="*60)
        print(f"遗传模型 | Genetic model: {self.config.genetic_model}")
        print(f"表型类型 | Trait type: {self.config.trait_type}")
        print(f"质控后样本数 | Samples after QC: {stats.get('total_samples', 'N/A')}")
        print(f"质控后SNP数 | SNPs after QC: {stats.get('total_snps', 'N/A')}")
        
        if self.config.trait_type == "qualitative":
            print(f"病例数 | Cases: {stats.get('cases', 'N/A')}")
            print(f"对照数 | Controls: {stats.get('controls', 'N/A')}")
        
        print("\n模型结果 | Model Results:")
        for model_name, results_df in model_results.items():
            if len(results_df) > 0:
                min_p = results_df['P'].min()
                print(f"  {model_name.upper()}: 最小P值 | Min P-value = {min_p:.2e}")
                
                # 显示显著位点数 | Show significant loci count
                if self.config.correction_method in ["bonferroni", "all"]:
                    bonferroni_threshold = self.config.bonferroni_alpha / len(results_df)
                    bonferroni_hits = len(results_df[results_df['P'] < bonferroni_threshold])
                    print(f"            Bonferroni显著位点 | Bonferroni significant: {bonferroni_hits}")
                
                if self.config.correction_method in ["suggestive", "all"]:
                    suggestive_hits = len(results_df[results_df['P'] < self.config.suggestive_threshold])
                    print(f"            提示性关联位点 | Suggestive loci: {suggestive_hits}")
        
        print(f"\n主要输出文件 | Main output files:")
        print(f"  - gwas_summary_report.txt: 总结报告 | Summary report")
        if self.config.genetic_model == "all":
            print(f"  - model_comparison_report.txt: 模型比较报告 | Model comparison report")
        print(f"  - gwas_results_*.txt: 详细关联分析结果 | Detailed association results")
        print("="*60)


def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='PLINK GWAS分析工具 | PLINK GWAS Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例 | Usage Examples:

# 基本质量性状分析（默认加性模型）| Basic qualitative trait analysis (default additive model)
python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t qualitative -o results

# 使用显性模型分析 | Analysis with dominant model
python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t qualitative -m dominant -o results

# 测试所有遗传模型 | Test all genetic models
python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t qualitative -m all -o results

# 数量性状分析 | Quantitative trait analysis
python run_plink_gwas.py -v data.vcf.gz -p pheno.txt -t quantitative -m additive -o results
        """
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('必需参数 | Required arguments')
    required.add_argument('-v', '--vcf-file', required=True,
                         help='VCF文件路径 | Path to VCF file')
    required.add_argument('-p', '--phenotype-file', required=True,
                         help='表型文件路径 | Path to phenotype file')
    
    # 分析参数 | Analysis parameters
    analysis = parser.add_argument_group('分析参数 | Analysis parameters')
    analysis.add_argument('-t', '--trait-type', choices=['qualitative', 'quantitative'],
                         default='qualitative',
                         help='表型类型 | Trait type (default: qualitative)')
    analysis.add_argument('-m', '--genetic-model', 
                         choices=['additive', 'dominant', 'recessive', 'all'],
                         default='additive',
                         help='遗传模型 | Genetic model (default: additive)')
    analysis.add_argument('-o', '--output-dir', default='gwas_results',
                         help='输出目录 | Output directory (default: gwas_results)')
    
    # 质量控制参数 | Quality control parameters
    qc = parser.add_argument_group('质量控制参数 | Quality control parameters')
    qc.add_argument('--mind', type=float, default=0.05,
                    help='个体缺失率阈值 | Individual missing rate threshold (default: 0.05)')
    qc.add_argument('--geno', type=float, default=0.05,
                    help='SNP缺失率阈值 | SNP missing rate threshold (default: 0.05)')
    qc.add_argument('--maf', type=float, default=0.01,
                    help='最小等位基因频率 | Minor allele frequency (default: 0.01)')
    qc.add_argument('--hwe', type=float, default=1e-6,
                    help='Hardy-Weinberg平衡P值阈值 | HWE p-value threshold (default: 1e-6)')
    
    # LD剪枝参数 | LD pruning parameters
    ld = parser.add_argument_group('LD剪枝参数 | LD pruning parameters')
    ld.add_argument('--ld-window-size', type=int, default=50,
                    help='LD剪枝窗口大小(kb) | LD window size in kb (default: 50)')
    ld.add_argument('--ld-step-size', type=int, default=5,
                    help='LD剪枝步长(SNP数) | LD step size in SNPs (default: 5)')
    ld.add_argument('--ld-r2-threshold', type=float, default=0.2,
                    help='LD剪枝r²阈值 | LD r² threshold (default: 0.2)')
    
    # 主成分分析参数 | PCA parameters
    pca = parser.add_argument_group('主成分分析参数 | PCA parameters')
    pca.add_argument('--pca-components', type=int, default=10,
                     help='计算的主成分数量 | Number of PCA components (default: 10)')
    pca.add_argument('--pca-use', type=int, default=5,
                     help='关联分析中使用的主成分数量 | Number of PCs to use in association (default: 5)')
    
    # 显著性校正参数 | Significance correction parameters
    correction = parser.add_argument_group('显著性校正参数 | Significance correction parameters')
    correction.add_argument('--correction-method', 
                           choices=['bonferroni', 'suggestive', 'fdr', 'all'],
                           default='all',
                           help='显著性校正方法 | Significance correction method (default: all)')
    correction.add_argument('--bonferroni-alpha', type=float, default=0.05,
                           help='Bonferroni校正alpha水平 | Bonferroni alpha level (default: 0.05)')
    correction.add_argument('--suggestive-threshold', type=float, default=1e-5,
                           help='提示性关联阈值 | Suggestive threshold (default: 1e-5)')
    correction.add_argument('--fdr-alpha', type=float, default=0.05,
                           help='FDR校正q值阈值 | FDR q-value threshold (default: 0.05)')
    
    # 计算参数 | Computing parameters
    compute = parser.add_argument_group('计算参数 | Computing parameters')
    compute.add_argument('--threads', type=int, default=1,
                        help='使用的线程数 | Number of threads (default: 1)')
    
    args = parser.parse_args()
    
    try:
        # 创建PLINK GWAS分析器 | Create PLINK GWAS analyzer
        gwas = PlinkGWAS(
            vcf_file=args.vcf_file,
            phenotype_file=args.phenotype_file,
            output_dir=args.output_dir,
            trait_type=args.trait_type,
            genetic_model=args.genetic_model,
            mind=args.mind,
            geno=args.geno,
            maf=args.maf,
            hwe=args.hwe,
            ld_window_size=args.ld_window_size,
            ld_step_size=args.ld_step_size,
            ld_r2_threshold=args.ld_r2_threshold,
            pca_components=args.pca_components,
            pca_use=args.pca_use,
            correction_method=args.correction_method,
            bonferroni_alpha=args.bonferroni_alpha,
            suggestive_threshold=args.suggestive_threshold,
            fdr_alpha=args.fdr_alpha,
            threads=args.threads
        )
        
        # 运行分析 | Run analysis
        gwas.run_analysis()
        
        print("\n✅ 分析成功完成! | Analysis completed successfully!")
        print(f"📁 结果保存在: {args.output_dir} | Results saved in: {args.output_dir}")
        
    except KeyboardInterrupt:
        print("\n⚠️  用户中断分析 | Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 分析失败 | Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()