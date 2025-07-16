"""
PLINK GWAS分析主程序模块 | PLINK GWAS Analysis Main Module
"""

import argparse
import sys
import os
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
        
        # 切换到工作目录 | Change to working directory
        os.chdir(self.config.output_path)
        
        # 初始化日志 | Initialize logging
        self.logger_manager = PlinkGWASLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化文件管理器 | Initialize file manager
        self.file_manager = FileManager(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.pheno_processor = PhenotypeProcessor(self.config, self.logger, self.cmd_runner)
        self.vcf_processor = VCFProcessor(self.config, self.logger, self.cmd_runner)
        self.quality_controller = QualityController(self.config, self.logger, self.cmd_runner)
        self.population_analyzer = PopulationAnalyzer(self.config, self.logger, self.cmd_runner)
        self.association_analyzer = AssociationAnalyzer(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.report_generator = ReportGenerator(self.config, self.logger)
        self.viz_generator = VisualizationGenerator(self.config, self.logger, self.cmd_runner)
    
    def run_analysis(self):
        """运行完整的PLINK GWAS分析流程 | Run complete PLINK GWAS analysis pipeline"""
        try:
            self.logger.info("开始PLINK GWAS分析流程 | Starting PLINK GWAS analysis pipeline...")
            self.logger.info("="*50)
            
            # 1. 检查并复制输入文件 | Check and copy input files
            self.logger.info("步骤1: 检查输入文件 | Step 1: Checking input files")
            self.file_manager.copy_input_files(self.config.vcf_file, self.config.phenotype_file)
            
            # 2. 转换表型文件 | Convert phenotype file
            self.logger.info("步骤2: 转换表型文件 | Step 2: Converting phenotype file")
            self.pheno_processor.convert_phenotype()
            
            # 3. 转换VCF文件 | Convert VCF file
            self.logger.info("步骤3: 转换VCF文件 | Step 3: Converting VCF file")
            self.vcf_processor.convert_vcf_to_plink()
            
            # 4. 合并表型信息 | Merge phenotype information
            self.logger.info("步骤4: 合并表型信息 | Step 4: Merging phenotype information")
            self.vcf_processor.merge_phenotype()
            
            # 5. 质量控制 | Quality control
            self.logger.info("步骤5: 质量控制 | Step 5: Quality control")
            self.quality_controller.quality_control()
            
            # 6. 群体分层控制 | Population stratification control
            self.logger.info("步骤6: 群体分层控制 | Step 6: Population stratification control")
            use_pca = self.population_analyzer.population_stratification()
            
            # 7. 关联分析 | Association analysis
            self.logger.info("步骤7: 关联分析 | Step 7: Association analysis")
            main_result = self.association_analyzer.association_analysis(use_pca)
            
            # 8. 处理结果 | Process results
            self.logger.info("步骤8: 处理结果 | Step 8: Processing results")
            results_df = self.results_processor.process_results(main_result)
            
            # 9. 生成报告 | Generate report
            self.logger.info("步骤9: 生成报告 | Step 9: Generating report")
            self.report_generator.generate_report(results_df, main_result)
            
            # 10. 生成可视化 | Generate visualization
            self.logger.info("步骤10: 生成可视化 | Step 10: Generating visualization")
            self.viz_generator.generate_plots()
            
            self.logger.info("="*50)
            self.logger.info("PLINK GWAS分析流程完成！ | PLINK GWAS analysis pipeline completed!")
            self.logger.info(f"所有结果保存在 | All results saved in: {self.config.output_path.absolute()}")
            self.logger.info("主要输出文件 | Main output files:")
            self.logger.info("- analysis_report.txt: 分析报告 | Analysis report")
            self.logger.info("- gwas_results_ADD.txt: 主要关联结果 | Main association results")
            self.logger.info("- significant_hits.txt: 显著关联位点 | Significant association loci")
            self.logger.info("- manhattan_plot.png: Manhattan图 | Manhattan plot")
            self.logger.info("- qq_plot.png: QQ图 | QQ plot")
            
        except Exception as e:
            self.logger.error(f"分析流程失败 | Analysis pipeline failed: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="完整的PLINK质量性状GWAS分析流程 (模块化版本) | Complete PLINK Qualitative Trait GWAS Analysis Pipeline (Modular Version)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 必需参数 | Required parameters
    parser.add_argument("-v", "--vcf-file", required=True,
                       help="输入VCF文件路径（支持.gz压缩） | Input VCF file path (supports .gz compression)")
    parser.add_argument("-p", "--phenotype-file", required=True,
                       help="表型文件路径（样本ID和表型值，以空格或制表符分隔） | Phenotype file path (sample ID and phenotype value, space or tab separated)")
    
    # 输出参数 | Output parameters
    parser.add_argument("-o", "--output-dir", default="plink_results",
                       help="输出目录 | Output directory")
    
    # 质控参数 | Quality control parameters
    parser.add_argument("--mind", type=float, default=0.05,
                       help="个体缺失率阈值（移除缺失率大于此值的个体） | Individual missing rate threshold")
    parser.add_argument("--geno", type=float, default=0.05,
                       help="SNP缺失率阈值（移除缺失率大于此值的SNP） | SNP missing rate threshold")
    parser.add_argument("--maf", type=float, default=0.01,
                       help="最小等位基因频率阈值 | Minor allele frequency threshold")
    parser.add_argument("--hwe", type=float, default=1e-6,
                       help="Hardy-Weinberg平衡检验P值阈值 | Hardy-Weinberg equilibrium p-value threshold")
    
    # LD剪枝参数 | LD pruning parameters
    parser.add_argument("--ld-window-size", type=int, default=50,
                       help="LD剪枝窗口大小（kb） | LD pruning window size (kb)")
    parser.add_argument("--ld-step-size", type=int, default=5,
                       help="LD剪枝步长（SNP数） | LD pruning step size (number of SNPs)")
    parser.add_argument("--ld-r2-threshold", type=float, default=0.2,
                       help="LD剪枝r²阈值 | LD pruning r² threshold")
    
    # 主成分参数 | PCA parameters
    parser.add_argument("--pca-components", type=int, default=10,
                       help="计算的主成分数量 | Number of principal components to compute")
    parser.add_argument("--pca-use", type=int, default=5,
                       help="关联分析中使用的主成分数量 | Number of PCs to use in association analysis")
    
    # 显著性阈值 | Significance thresholds
    parser.add_argument("--significance-threshold", default="5e-8",
                       help="全基因组显著性阈值 | Genome-wide significance threshold")
    parser.add_argument("--suggestive-threshold", default="1e-5",
                       help="提示性关联阈值 | Suggestive association threshold")
    
    # 其他参数 | Other parameters
    parser.add_argument("--threads", type=int, default=1,
                       help="使用的线程数 | Number of threads to use")
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    plink = PlinkGWAS(**vars(args))
    plink.run_analysis()

if __name__ == "__main__":
    main()
