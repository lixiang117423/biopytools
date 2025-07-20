"""
ADMIXTURE分析主程序模块 | ADMIXTURE Analysis Main Module
"""

import argparse
import sys
from .config import AdmixtureConfig
from .utils import AdmixtureLogger, CommandRunner, SoftwareChecker
from .data_processing import VCFProcessor, PlinkProcessor
# 修正 #1: 导入 AdmixtureAnalyzer 并使用别名 CoreAnalyzer 避免命名冲突
from .analysis import AdmixtureAnalyzer as CoreAnalyzer, ResultsProcessor 
from .results import CovariateGenerator, PlotGenerator, SummaryGenerator

class AdmixtureAnalyzer:
    """ADMIXTURE分析主类 | Main ADMIXTURE Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = AdmixtureConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = AdmixtureLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 检查软件环境 | Check software environment
        self.software_checker = SoftwareChecker(self.logger)
        if not self.software_checker.check_dependencies():
            sys.exit(1)
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.vcf_processor = VCFProcessor(self.config, self.logger, self.cmd_runner)
        self.plink_processor = PlinkProcessor(self.config, self.logger, self.cmd_runner)
        # 修正 #1: 使用别名 CoreAnalyzer 进行实例化
        self.admixture_analyzer = CoreAnalyzer(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.covariate_generator = CovariateGenerator(self.config, self.logger)
        self.plot_generator = PlotGenerator(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的ADMIXTURE分析流程 | Run complete ADMIXTURE analysis pipeline"""
        try:
            self.logger.info("="*60)
            self.logger.info("开始ADMIXTURE群体结构分析 | Starting ADMIXTURE Population Structure Analysis")
            self.logger.info("="*60)
            
            # 步骤1: VCF预处理 | Step 1: VCF preprocessing
            self.logger.info("\n步骤1: VCF文件预处理 | Step 1: VCF file preprocessing")
            processed_vcf = self.vcf_processor.preprocess_vcf()
            if processed_vcf is None:
                processed_vcf = self.config.vcf_file
            
            # 步骤2: 转换为PLINK格式 | Step 2: Convert to PLINK format
            self.logger.info("\n步骤2: 转换为PLINK格式 | Step 2: Convert to PLINK format")
            raw_prefix = self.plink_processor.convert_vcf_to_plink(processed_vcf)
            
            # 步骤3: 质量控制 | Step 3: Quality control
            self.logger.info("\n步骤3: 质量控制 | Step 3: Quality control")
            qc_prefix = self.plink_processor.quality_control(raw_prefix)
            
            # 步骤4: ADMIXTURE分析 | Step 4: ADMIXTURE analysis
            self.logger.info("\n步骤4: ADMIXTURE分析 | Step 4: ADMIXTURE analysis")
            best_k = self.admixture_analyzer.run_admixture_analysis(qc_prefix)
            
            # 步骤5: 结果处理 | Step 5: Results processing
            self.logger.info("\n步骤5: 结果处理 | Step 5: Results processing")
            q_data, p_data, stats = self.results_processor.process_results(best_k)
            
            # 步骤6: 生成GWAS协变量 | Step 6: Generate GWAS covariates
            self.logger.info("\n步骤6: 生成GWAS协变量 | Step 6: Generate GWAS covariates")
            covar_file = self.covariate_generator.generate_gwas_covariates(best_k)
            
            # 步骤7: 生成可视化图表 | Step 7: Generate visualization plots
            self.logger.info("\n步骤7: 生成可视化图表 | Step 7: Generate visualization plots")
            plot_file = self.plot_generator.generate_plots(q_data, best_k)
            
            # 步骤8: 生成总结报告 | Step 8: Generate summary report
            self.logger.info("\n步骤8: 生成总结报告 | Step 8: Generate summary report")
            summary_file = self.summary_generator.generate_summary(best_k, stats)
            
            # 完成信息 | Completion information
            self.logger.info("\n" + "="*60)
            self.logger.info("ADMIXTURE分析完成！| ADMIXTURE analysis completed!")
            self.logger.info("="*60)
            self.logger.info(f"最优K值 | Best K value: {best_k}")
            self.logger.info(f"结果目录 | Results directory: {self.config.output_dir}")
            self.logger.info("主要输出文件 | Main output files:")
            self.logger.info("  - admixture_proportions.csv: 个体祖先成分 | Individual ancestry proportions")
            self.logger.info("  - gwas_covariates.txt: GWAS协变量文件 | GWAS covariate file")
            self.logger.info("  - cv_results.csv: 交叉验证结果 | Cross-validation results")
            self.logger.info("  - *.pdf: 可视化图表 | Visualization plots")
            self.logger.info("  - analysis_summary.txt: 分析总结报告 | Analysis summary report")
            
        except Exception as e:
            self.logger.error(f"分析过程中发生错误 | Error occurred during analysis: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="ADMIXTURE群体结构分析工具 (模块化版本) | ADMIXTURE Population Structure Analysis Tool (Modular Version)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 必需参数 | Required arguments
    parser.add_argument("-v", "--vcf", required=True,
                       help="输入VCF文件路径 | Input VCF file path")
    
    # 可选参数 | Optional arguments
    parser.add_argument("-o", "--output", default="admixture_results",
                       help="输出目录 | Output directory")
    parser.add_argument("-k", "--min-k", type=int, default=2,
                       help="最小K值 | Minimum K value")
    parser.add_argument("-K", "--max-k", type=int, default=10,
                       help="最大K值 | Maximum K value")
    parser.add_argument("-c", "--cv-folds", type=int, default=5,
                       help="交叉验证折数 | Cross-validation folds")
    parser.add_argument("-t", "--threads", type=int, default=4,
                       help="线程数 | Number of threads")
    
    # 质控参数 | Quality control parameters
    parser.add_argument("-m", "--maf", type=float, default=0.01,
                       help="MAF阈值 | MAF threshold")
    parser.add_argument("-M", "--missing", type=float, default=0.1,
                       help="缺失率阈值 | Missing rate threshold")
    parser.add_argument("-H", "--hwe", type=float, default=1e-6,
                       help="HWE p值阈值 | HWE p-value threshold")
    
    # 处理选项 | Processing options
    parser.add_argument("-s", "--skip-preprocessing", action="store_true",
                       help="跳过VCF预处理 | Skip VCF preprocessing")
    parser.add_argument("-i", "--keep-intermediate", action="store_true",
                       help="保留中间文件 | Keep intermediate files")
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = AdmixtureAnalyzer(
        vcf_file=args.vcf,
        output_dir=args.output,
        min_k=args.min_k,
        max_k=args.max_k,
        cv_folds=args.cv_folds,
        threads=args.threads,
        maf=args.maf,
        missing_rate=args.missing,
        hwe_pvalue=args.hwe,
        skip_preprocessing=args.skip_preprocessing,
        keep_intermediate=args.keep_intermediate
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
