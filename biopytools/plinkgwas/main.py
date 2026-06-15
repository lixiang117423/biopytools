
# 20250726添加显性和隐性表型类型支持|Added support for dominant and recessive trait types on 20250726
"""
PLINK GWAS分析主程序模块|PLINK GWAS Analysis Main Module
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
    """PLINK GWAS分析主类|Main PLINK GWAS Analysis Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = PlinkGWASConfig(**kwargs)
        self.config.validate()

        # 确保输出目录存在|Ensure output directory exists
        self.config.output_path.mkdir(parents=True, exist_ok=True)

        # 保存原始工作目录|Save original working directory
        self.original_cwd = Path.cwd()

        # 初始化日志（在切换目录之前）| Initialize logging (before changing directory)
        # 根据verbose/quiet参数确定日志级别|Determine log level based on verbose/quiet
        import logging
        if self.config.quiet:
            log_level = logging.ERROR
        elif self.config.verbose >= 2:
            log_level = logging.DEBUG
        else:
            # 默认显示INFO级别|Default to INFO level
            log_level = logging.INFO

        self.logger_manager = PlinkGWASLogger(
            self.config.output_path,
            log_name="plink_analysis.log",
            log_level=log_level
        )
        self.logger = self.logger_manager.get_logger()
        
        # 切换到工作目录|Change to working directory
        os.chdir(self.config.output_path)
        self.logger.info(f"工作目录|Working directory: {self.config.output_path}")
        
        # 初始化命令执行器和文件管理器（使用当前目录"."）| Initialize command runner and file manager (use current directory ".")
        current_dir = Path(".")
        self.cmd_runner = CommandRunner(self.logger, current_dir)
        self.file_manager = FileManager(self.logger, current_dir)
        
        # 初始化各个处理器|Initialize processors
        self.pheno_processor = PhenotypeProcessor(self.config, self.logger, self.cmd_runner)
        self.vcf_processor = VCFProcessor(self.config, self.logger, self.cmd_runner)
        self.quality_controller = QualityController(self.config, self.logger, self.cmd_runner)
        self.population_analyzer = PopulationAnalyzer(self.config, self.logger, self.cmd_runner)
        self.association_analyzer = AssociationAnalyzer(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.report_generator = ReportGenerator(self.config, self.logger)
        self.viz_generator = VisualizationGenerator(self.config, self.logger, self.cmd_runner)
    
    def __del__(self):
        """析构函数，恢复原始工作目录|Destructor, restore original working directory"""
        try:
            if hasattr(self, 'original_cwd'):
                os.chdir(self.original_cwd)
        except:
            pass
    
    def run_analysis(self):
        """运行完整的PLINK GWAS分析流程|Run complete PLINK GWAS analysis pipeline"""
        try:
            # 程序启动信息|Program startup information
            self.logger.info("=" * 60)
            self.logger.info("PLINK GWAS Analysis Pipeline|PLINK GWAS分析流程")
            self.logger.info("Version: 1.0.0")
            self.logger.info("=" * 60)
            self.logger.info(f"Input VCF|VCF输入文件: {self.config.vcf_file}")
            self.logger.info(f"Input phenotype|表型输入文件: {self.config.phenotype_file}")
            self.logger.info(f"Output directory|输出目录: {self.config.output_dir}")
            self.logger.info(f"Trait type|表型类型: {self.config.trait_type}")
            self.logger.info(f"Genetic model|遗传模型: {self.config.genetic_model}")
            self.logger.info(f"Threads|线程数: {self.config.threads}")
            self.logger.info(f"Correction method|显著性校正方法: {self.config.correction_method}")
            self.logger.info("=" * 60)

            self.logger.info("开始PLINK GWAS分析流程|Starting PLINK GWAS analysis pipeline...")
            self.logger.info("=" * 60)

            # 1. 检查并复制输入文件|Check and copy input files
            self.logger.info("=" * 60)
            self.logger.info("STEP 1: 检查输入文件|Checking input files")
            self.logger.info("=" * 60)
            # 构建绝对路径|Build absolute paths
            vcf_path = self.original_cwd / self.config.vcf_file if not Path(self.config.vcf_file).is_absolute() else Path(self.config.vcf_file)
            pheno_path = self.original_cwd / self.config.phenotype_file if not Path(self.config.phenotype_file).is_absolute() else Path(self.config.phenotype_file)

            self.file_manager.copy_input_files(str(vcf_path), str(pheno_path))

            # 2. 转换表型文件|Convert phenotype file
            self.logger.info("=" * 60)
            self.logger.info("STEP 2: 转换表型文件|Converting phenotype file")
            self.logger.info("=" * 60)
            self.pheno_processor.convert_phenotype()

            # 3. 转换VCF文件|Convert VCF file
            self.logger.info("=" * 60)
            self.logger.info("STEP 3: 转换VCF文件|Converting VCF file")
            self.logger.info("=" * 60)
            self.vcf_processor.convert_vcf_to_plink()

            # 4. 合并表型信息|Merge phenotype information
            self.logger.info("=" * 60)
            self.logger.info("STEP 4: 合并表型信息|Merging phenotype information")
            self.logger.info("=" * 60)
            self.vcf_processor.merge_phenotype()

            # 5. 数据质量控制|Data quality control
            self.logger.info("=" * 60)
            self.logger.info("STEP 5: 数据质量控制|Data quality control")
            self.logger.info("=" * 60)
            self.quality_controller.quality_control()

            # 6. 群体结构分析|Population structure analysis
            self.logger.info("=" * 60)
            self.logger.info("STEP 6: 群体结构分析|Population structure analysis")
            self.logger.info("=" * 60)
            self.population_analyzer.population_analysis()

            # 7. 关联分析|Association analysis
            self.logger.info("=" * 60)
            self.logger.info("STEP 7: 关联分析|Association analysis")
            self.logger.info("=" * 60)
            main_result = self.association_analyzer.association_analysis(use_pca=True)

            # 8. 结果处理|Results processing
            self.logger.info("=" * 60)
            self.logger.info("STEP 8: 结果处理|Results processing")
            self.logger.info("=" * 60)
            model_results = self.results_processor.process_results(main_result)

            # 9. 生成报告|Generate reports
            self.logger.info("=" * 60)
            self.logger.info("STEP 9: 生成报告|Generating reports")
            self.logger.info("=" * 60)
            stats = self._collect_statistics()
            self.report_generator.generate_summary_report(model_results, stats)

            # 10. 生成可视化|Generate visualizations
            self.logger.info("=" * 60)
            self.logger.info("STEP 10: 生成可视化|Generating visualizations")
            self.logger.info("=" * 60)
            self.viz_generator.generate_plots(model_results)

            self.logger.info("=" * 60)
            self.logger.info("PLINK GWAS分析完成!|PLINK GWAS analysis completed!")
            self.logger.info("=" * 60)
            
            # 输出结果摘要|Output results summary
            self._print_results_summary(model_results, stats)
            
        except Exception as e:
            self.logger.error(f"分析过程中发生错误|Error during analysis: {e}")
            raise
        finally:
            # 确保最后恢复原始工作目录|Ensure original working directory is restored
            try:
                os.chdir(self.original_cwd)
            except:
                pass
    
    def _collect_statistics(self) -> dict:
        """收集分析统计信息|Collect analysis statistics"""
        stats = {}
        
        try:
            # 读取样本信息|Read sample information
            if os.path.exists("data_qc1.fam"):
                with open("data_qc1.fam", 'r') as f:
                    lines = f.readlines()
                    stats['total_samples'] = len(lines)
                    
                    # 如果是质量性状，统计病例对照数|Count cases/controls for qualitative traits
                    if self.config.trait_type == "qualitative":
                        cases = sum(1 for line in lines if line.strip().split()[-1] == '2')
                        controls = sum(1 for line in lines if line.strip().split()[-1] == '1')
                        stats['cases'] = cases
                        stats['controls'] = controls
            
            # 读取SNP信息|Read SNP information
            if os.path.exists("data_qc1.bim"):
                with open("data_qc1.bim", 'r') as f:
                    lines = f.readlines()
                    stats['total_snps'] = len(lines)
                    
                    # 统计染色体数|Count chromosomes
                    chromosomes = set(line.strip().split()[0] for line in lines)
                    stats['chromosomes'] = len(chromosomes)
            
        except Exception as e:
            self.logger.warning(f"收集统计信息时出错|Error collecting statistics: {e}")
        
        return stats
    
    def _print_results_summary(self, model_results: dict, stats: dict):
        """打印结果摘要|Print results summary"""
        print("="*60)
        print("PLINK GWAS 分析结果摘要|Analysis Results Summary")
        print("="*60)
        print(f"遗传模型|Genetic model: {self.config.genetic_model}")
        print(f"表型类型|Trait type: {self.config.trait_type}")
        print(f"质控后样本数|Samples after QC: {stats.get('total_samples', 'N/A')}")
        print(f"质控后SNP数|SNPs after QC: {stats.get('total_snps', 'N/A')}")
        
        if self.config.trait_type == "qualitative":
            print(f"病例数|Cases: {stats.get('cases', 'N/A')}")
            print(f"对照数|Controls: {stats.get('controls', 'N/A')}")

        print("模型结果|Model Results:")
        for model_name, results_df in model_results.items():
            if len(results_df) > 0:
                min_p = results_df['P'].min()
                print(f"  {model_name.upper()}: 最小P值|Min P-value = {min_p:.2e}")
                
                # 显示显著位点数|Show significant loci count
                if self.config.correction_method in ["bonferroni", "all"]:
                    bonferroni_threshold = self.config.bonferroni_alpha / len(results_df)
                    bonferroni_hits = len(results_df[results_df['P'] < bonferroni_threshold])
                    print(f"            Bonferroni显著位点|Bonferroni significant: {bonferroni_hits}")
                
                if self.config.correction_method in ["suggestive", "all"]:
                    suggestive_hits = len(results_df[results_df['P'] < self.config.suggestive_threshold])
                    print(f"            提示性关联位点|Suggestive loci: {suggestive_hits}")
        
        print(f"\n主要输出文件|Main output files:")
        print(f"  - gwas_summary_report.txt: 总结报告|Summary report")
        if self.config.genetic_model == "all":
            print(f"  - model_comparison_report.txt: 模型比较报告|Model comparison report")
        print(f"  - gwas_results_*.txt: 详细关联分析结果|Detailed association results")
        print("="*60)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='PLINK GWAS分析工具|PLINK GWAS Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  python run_plink_gwas.py -i data.vcf.gz -p pheno.txt -o results
        """
    )
    
    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-i', '--vcf', required=True,
                         help='VCF文件路径|Path to VCF file')
    required.add_argument('-p', '--phenotype', required=True,
                         help='表型文件路径|Path to phenotype file')
    
    # 分析参数|Analysis parameters
    analysis = parser.add_argument_group('分析参数|Analysis parameters')
    analysis.add_argument('-T', '--trait-type', choices=['qualitative', 'quantitative'],
                         default='qualitative',
                         help='表型类型|Trait type ')
    analysis.add_argument('-m', '--genetic-model', 
                         choices=['additive', 'dominant', 'recessive', 'all'],
                         default='additive',
                         help='遗传模型|Genetic model ')
    analysis.add_argument('-o', '--output-dir', default='gwas_results',
                         help='输出目录|Output directory ')
    # --- 新增代码开始 ---
    analysis.add_argument('--no-strat-corr', 
                         action='store_false', 
                         dest='population_strat_correction',
                         help='禁用群体结构校正（跳过LD剪枝和PCA）| Disable population stratification correction (skip LD pruning and PCA)')
    # --- 新增代码结束 ---
    
    # 质量控制参数|Quality control parameters
    qc = parser.add_argument_group('质量控制参数|Quality control parameters')
    qc.add_argument('--mind', type=float, default=None,
                    help='个体缺失率阈值|Individual missing rate threshold ')
    qc.add_argument('--geno', type=float, default=0.05,
                    help='SNP缺失率阈值|SNP missing rate threshold ')
    qc.add_argument('--maf', type=float, default=0.05,
                    help='最小等位基因频率|Minor allele frequency ')
    qc.add_argument('--hwe', type=float, default=None,
                    help='Hardy-Weinberg平衡P值阈值|HWE p-value threshold ')
    
    # LD剪枝参数|LD pruning parameters
    ld = parser.add_argument_group('LD剪枝参数|LD pruning parameters')
    ld.add_argument('--ld-window-size', type=int, default=50,
                    help='LD剪枝窗口大小(kb)|LD window size in kb ')
    ld.add_argument('--ld-step-size', type=int, default=5,
                    help='LD剪枝步长(SNP数)|LD step size in SNPs ')
    ld.add_argument('--ld-r2-threshold', type=float, default=0.2,
                    help='LD剪枝r²阈值|LD r² threshold ')
    
    # 主成分分析参数|PCA parameters
    pca = parser.add_argument_group('主成分分析参数|PCA parameters')
    pca.add_argument('--pca-components', type=int, default=10,
                     help='计算的主成分数量|Number of PCA components ')
    pca.add_argument('--pca-use', type=int, default=5,
                     help='关联分析中使用的主成分数量|Number of PCs to use in association ')
    
    # 显著性校正参数|Significance correction parameters
    correction = parser.add_argument_group('显著性校正参数|Significance correction parameters')
    correction.add_argument('--correction-method', 
                           choices=['bonferroni', 'suggestive', 'fdr', 'all'],
                           default='all',
                           help='显著性校正方法|Significance correction method ')
    correction.add_argument('--bonferroni-alpha', type=float, default=0.05,
                           help='Bonferroni校正alpha水平|Bonferroni alpha level ')
    correction.add_argument('--suggestive-threshold', type=float, default=1e-5,
                           help='提示性关联阈值|Suggestive threshold ')
    correction.add_argument('--fdr-alpha', type=float, default=0.05,
                           help='FDR校正q值阈值|FDR q-value threshold ')
    
    # 计算参数|Computing parameters
    compute = parser.add_argument_group('计算参数|Computing parameters')
    compute.add_argument('-t', '--threads', type=int, default=1,
                        help='使用的线程数|Number of threads ')

    # 日志控制参数|Logging control parameters
    logging_group = parser.add_argument_group('日志控制参数|Logging control parameters')
    logging_group.add_argument('-v', '--verbose', action='count', default=0,
                              help='详细输出模式|Verbose mode (-v: INFO, -vv: DEBUG)')
    logging_group.add_argument('--quiet', action='store_true',
                              help='静默模式，只输出ERROR|Quiet mode (only ERROR)')
    logging_group.add_argument('--log-file',
                              help='日志文件路径|Log file path')

    # 执行控制参数|Execution control parameters
    exec_group = parser.add_argument_group('执行控制参数|Execution control parameters')
    exec_group.add_argument('-f', '--force', action='store_true',
                           help='强制覆盖已存在的输出目录|Force overwrite existing output directory')
    exec_group.add_argument('--dry-run', action='store_true',
                           help='模拟运行，不实际执行分析|Dry run without actual analysis')

    # 版本信息|Version information
    parser.add_argument('-V', '--version', action='version',
                       version='%(prog)s 1.0.0')

    args = parser.parse_args()

    # --- 新增代码开始 ---
    # 确保在创建实例前设置好默认值
    # 这一步可以省略，因为 dataclass 已经有默认值，但这样做更明确
    if 'population_strat_correction' not in args:
        args.population_strat_correction = True
    # --- 新增代码结束 ---
    
    try:
        # 创建PLINK GWAS分析器|Create PLINK GWAS analyzer
        gwas = PlinkGWAS(
            vcf_file=args.vcf,
            phenotype_file=args.phenotype,
            output_dir=args.output_dir,
            trait_type=args.trait_type,
            genetic_model=args.genetic_model,
            population_strat_correction=args.population_strat_correction,
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
            threads=args.threads,
            verbose=args.verbose,
            quiet=args.quiet,
            log_file=args.log_file if hasattr(args, 'log_file') else None
        )

        # 运行分析|Run analysis
        gwas.run_analysis()

        print("分析成功完成!|Analysis completed successfully!")
        print(f"结果保存在: {args.output_dir}|Results saved in: {args.output_dir}")

    except KeyboardInterrupt:
        print("用户中断分析|Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"分析失败|Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()