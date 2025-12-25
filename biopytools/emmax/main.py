"""
🧬 GWAS分析主程序模块 | GWAS Analysis Main Module
"""

import argparse
import sys
import traceback

try:
    from .config import GWASConfig
    from .gwas_engine import GWASEngine
    from .population_analysis import PopulationAnalyzer
    from .report_generator import ReportGenerator
    from .utils import CommandRunner, GWASLogger, check_dependencies, cleanup_temp_files
    from .vcf_processor import VCFProcessor
    from .visualization import GWASVisualizer
except ImportError:
    # 如果是直接运行单文件，使用绝对导入
    from config import GWASConfig
    from gwas_engine import GWASEngine
    from population_analysis import PopulationAnalyzer
    from report_generator import ReportGenerator
    from utils import CommandRunner, GWASLogger, check_dependencies, cleanup_temp_files
    from vcf_processor import VCFProcessor
    from visualization import GWASVisualizer


class GWASAnalyzer:
    """GWAS分析主类 | Main GWAS Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = GWASConfig(**kwargs)
        self.config.validate()

        # 初始化日志 | Initialize logging
        self.logger_manager = GWASLogger(self.config.output_prefix)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.working_path)

        # 初始化各个处理器 | Initialize processors
        self.vcf_processor = VCFProcessor(self.config, self.logger, self.cmd_runner)
        self.population_analyzer = PopulationAnalyzer(
            self.config, self.logger, self.cmd_runner
        )
        self.gwas_engine = GWASEngine(self.config, self.logger, self.cmd_runner)
        self.visualizer = GWASVisualizer(self.config, self.logger)
        self.report_generator = ReportGenerator(self.config, self.logger)

        # 临时文件列表 | Temporary files list
        self.temp_dirs = []

    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)

    def run_analysis(self):
        """运行完整的GWAS分析流程 | Run complete GWAS analysis pipeline"""
        try:
            self.logger.info("=" * 80)
            self.logger.info("🧬 开始AutoGWAS分析 | Starting AutoGWAS Analysis")
            self.logger.info("=" * 80)

            # 步骤0: 检查依赖 | Step 0: Check dependencies
            if not self.check_dependencies():
                raise RuntimeError("依赖检查失败 | Dependency check failed")

            # 步骤1: VCF文件处理 | Step 1: VCF file processing
            self.logger.info("-" * 40)
            self.logger.info("步骤1: VCF文件处理 | Step 1: VCF File Processing")
            self.logger.info("-" * 40)

            if not self.vcf_processor.process_vcf():
                raise RuntimeError("VCF文件处理失败 | VCF file processing failed")

            processed_files = self.vcf_processor.get_processed_files()
            self.temp_dirs.append(processed_files["temp_dir"])

            # 步骤2: 群体结构分析 | Step 2: Population structure analysis
            self.logger.info("-" * 40)
            self.logger.info(
                "步骤2: 群体结构分析 | Step 2: Population Structure Analysis"
            )
            self.logger.info("-" * 40)

            if not self.population_analyzer.analyze_population_structure(
                processed_files["plink_prefix"]
            ):
                raise RuntimeError(
                    "群体结构分析失败 | Population structure analysis failed"
                )

            analysis_files = self.population_analyzer.get_analysis_files()
            self.temp_dirs.append(analysis_files["temp_dir"])

            # 步骤3: GWAS关联分析 | Step 3: GWAS association analysis
            self.logger.info("-" * 40)
            self.logger.info("步骤3: GWAS关联分析 | Step 3: GWAS Association Analysis")
            self.logger.info("-" * 40)

            if not self.gwas_engine.run_gwas_analysis(
                processed_files["plink_prefix"],
                processed_files["phenotype_file"],
                analysis_files["kinship_prefix"],
            ):
                raise RuntimeError(
                    "GWAS关联分析失败 | GWAS association analysis failed"
                )

            results_file = self.gwas_engine.get_results_file()
            self.temp_dirs.append(self.gwas_engine.temp_dir)

            # 步骤4: 结果可视化 | Step 4: Results visualization
            self.logger.info("-" * 40)
            self.logger.info("步骤4: 结果可视化 | Step 4: Results Visualization")
            self.logger.info("-" * 40)

            if not self.visualizer.create_manhattan_plot(results_file):
                raise RuntimeError("曼哈顿图创建失败 | Manhattan plot creation failed")

            if not self.visualizer.create_qq_plot(results_file):
                raise RuntimeError("QQ图创建失败 | QQ plot creation failed")

            if not self.visualizer.extract_significant_snps(results_file):
                raise RuntimeError(
                    "显著SNP提取失败 | Significant SNP extraction failed"
                )

            # 步骤5: 生成报告 | Step 5: Generate reports
            self.logger.info("-" * 40)
            self.logger.info("步骤5: 生成报告 | Step 5: Generate Reports")
            self.logger.info("-" * 40)

            if not self.report_generator.generate_excel_report(
                results_file, self.config.significant_snps
            ):
                raise RuntimeError("Excel报告生成失败 | Excel report generation failed")

            if not self.report_generator.generate_html_report(
                results_file, self.config.significant_snps
            ):
                raise RuntimeError("HTML报告生成失败 | HTML report generation failed")

            self.logger.info("=" * 80)
            self.logger.info(
                "✅ AutoGWAS分析完成 | AutoGWAS Analysis Completed Successfully"
            )
            self.logger.info("=" * 80)

            # 输出结果文件信息 | Output result file information
            self._log_output_files()

        except Exception as e:
            self.logger.error(f"❌ 分析失败 | Analysis failed: {e}")
            self.logger.error(
                f"详细错误信息 | Detailed error: {traceback.format_exc()}"
            )
            raise

        finally:
            # 清理临时文件 | Cleanup temporary files
            cleanup_temp_files(self.temp_dirs, self.logger)

    def _log_output_files(self):
        """记录输出文件信息 | Log output file information"""
        self.logger.info("📁 输出文件 | Output files:")

        files_to_check = [
            (self.config.manhattan_plot, "曼哈顿图 | Manhattan plot"),
            (self.config.qq_plot, "QQ图 | QQ plot"),
            (self.config.significant_snps, "显著SNP | Significant SNPs"),
            (self.config.excel_report, "Excel报告 | Excel report"),
            (self.config.html_report, "HTML报告 | HTML report"),
            (f"{self.config.output_prefix}.log", "日志文件 | Log file"),
        ]

        for file_path, description in files_to_check:
            if os.path.exists(file_path):
                size = os.path.getsize(file_path)
                self.logger.info(f"  - {description}: {file_path} ({size:,} bytes)")


def create_parser():
    """创建命令行参数解析器 | Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description="🧬 AutoGWAS: 自动化GWAS分析工具 v1.0.0 | AutoGWAS: Automated GWAS Analysis Tool v1.0.0",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  # 基本分析 | Basic analysis
  %(prog)s -v input.vcf -p phenotypes.txt -o gwas_results
  
  # 指定完整路径 | Specify full paths
  %(prog)s --vcf 01.data/variants.vcf.gz \\
           --phenotype 02.pheno/traits.txt \\
           --output 03.results/gwas_analysis
  
  # 自定义分析参数 | Custom analysis parameters
  %(prog)s -v data.vcf -p phenotype.txt -o results \\
           --maf-threshold 0.01 \\
           --missing-threshold 0.1 \\
           --p-value-threshold 1e-8
  
  # 自定义工具路径 | Custom tool paths
  %(prog)s -v input.vcf -p traits.txt -o results \\
           --plink-path /path/to/plink \\
           --bcftools-path /path/to/bcftools \\
           --emmax-path /path/to/emmax
        """,
    )

    # 输入文件参数 | Input file arguments
    input_group = parser.add_argument_group("输入文件 | Input Files")
    input_group.add_argument(
        "-v", "--vcf", dest="vcf_file", help="输入VCF文件路径 | Input VCF file path"
    )
    input_group.add_argument(
        "-p",
        "--phenotype",
        dest="phenotype_file",
        help="输入表型文件路径 | Input phenotype file path",
    )

    # 输出文件参数 | Output file arguments
    output_group = parser.add_argument_group("输出文件 | Output Files")
    output_group.add_argument(
        "-o",
        "--output",
        dest="output_prefix",
        default="gwas_analysis",
        help="输出文件前缀（默认: gwas_analysis）| Output file prefix (default: gwas_analysis)",
    )
    output_group.add_argument(
        "--manhattan-plot",
        help="曼哈顿图输出文件路径（默认: OUTPUT_manhattan.png）| Manhattan plot output file path (default: OUTPUT_manhattan.png)",
    )
    output_group.add_argument(
        "--qq-plot",
        help="QQ图输出文件路径（默认: OUTPUT_qq.png）| QQ plot output file path (default: OUTPUT_qq.png)",
    )
    output_group.add_argument(
        "--significant-snps",
        help="显著SNP输出文件路径（默认: OUTPUT_significant_snps.txt）| Significant SNPs output file path (default: OUTPUT_significant_snps.txt)",
    )
    output_group.add_argument(
        "--excel-report",
        help="Excel报告输出文件路径（默认: OUTPUT_report.xlsx）| Excel report output file path (default: OUTPUT_report.xlsx)",
    )
    output_group.add_argument(
        "--html-report",
        help="HTML报告输出文件路径（默认: OUTPUT_report.html）| HTML report output file path (default: OUTPUT_report.html)",
    )

    # 工具路径参数 | Tool path arguments
    tool_group = parser.add_argument_group("工具路径 | Tool Paths")
    tool_group.add_argument(
        "--plink-path",
        default="plink",
        help="PLINK程序路径（默认: plink）| PLINK program path (default: plink)",
    )
    tool_group.add_argument(
        "--bcftools-path",
        default="bcftools",
        help="BCFtools程序路径（默认: bcftools）| BCFtools program path (default: bcftools)",
    )
    tool_group.add_argument(
        "--admixture-path",
        default="admixture",
        help="ADMIXTURE程序路径（默认: admixture）| ADMIXTURE program path (default: admixture)",
    )
    tool_group.add_argument(
        "--emmax-path",
        default="emmax",
        help="EMMAX程序路径（默认: emmax）| EMMAX program path (default: emmax)",
    )

    # VCF过滤参数 | VCF filtering parameters
    filter_group = parser.add_argument_group("VCF过滤参数 | VCF Filtering Parameters")
    filter_group.add_argument(
        "--maf-threshold",
        type=float,
        default=0.05,
        help="最小等位基因频率阈值（默认: 0.05）| Minor allele frequency threshold (default: 0.05)",
    )
    filter_group.add_argument(
        "--missing-threshold",
        type=float,
        default=0.2,
        help="缺失率阈值（默认: 0.2）| Missing rate threshold (default: 0.2)",
    )
    filter_group.add_argument(
        "--depth-min",
        type=int,
        default=3,
        help="最小测序深度（默认: 3）| Minimum sequencing depth (default: 3)",
    )
    filter_group.add_argument(
        "--depth-max",
        type=int,
        default=50,
        help="最大测序深度（默认: 50）| Maximum sequencing depth (default: 50)",
    )
    filter_group.add_argument(
        "--qual-min",
        type=float,
        default=20.0,
        help="最小质量值（默认: 20.0）| Minimum quality value (default: 20.0)",
    )

    # 群体结构分析参数 | Population structure analysis parameters
    pop_group = parser.add_argument_group(
        "群体结构分析参数 | Population Structure Analysis Parameters"
    )
    pop_group.add_argument(
        "--admixture-k-range",
        nargs=2,
        type=int,
        default=[1, 20],
        help="Admixture K值范围（默认: 1 20）| Admixture K value range (default: 1 20)",
    )
    pop_group.add_argument(
        "--pca-components",
        type=int,
        default=10,
        help="PCA主成分数量（默认: 10）| PCA components count (default: 10)",
    )
    pop_group.add_argument(
        "--ld-window",
        type=int,
        default=50,
        help="LD窗口大小（默认: 50）| LD window size (default: 50)",
    )
    pop_group.add_argument(
        "--ld-step",
        type=int,
        default=10,
        help="LD步长（默认: 10）| LD step size (default: 10)",
    )
    pop_group.add_argument(
        "--ld-r2",
        type=float,
        default=0.2,
        help="LD r²阈值（默认: 0.2）| LD r² threshold (default: 0.2)",
    )

    # EMMAX分析参数 | EMMAX analysis parameters
    emmax_group = parser.add_argument_group("EMMAX分析参数 | EMMAX Analysis Parameters")
    emmax_group.add_argument(
        "--emmax-precision",
        type=int,
        default=5,
        help="EMMAX输出精度（默认: 5）| EMMAX output precision (default: 5)",
    )
    emmax_group.add_argument(
        "--p-value-threshold",
        type=float,
        default=5e-8,
        help="显著性p值阈值（默认: 5e-8）| Significance p-value threshold (default: 5e-8)",
    )

    # 工作目录参数 | Working directory parameters
    work_group = parser.add_argument_group("工作目录 | Working Directory")
    work_group.add_argument(
        "-w",
        "--working-dir",
        default=".",
        help="工作目录（默认: 当前目录）| Working directory (default: current directory)",
    )

    return parser


def main():
    """主函数 | Main function"""
    parser = create_parser()
    args = parser.parse_args()

    # 检查必需参数 | Check required arguments
    if not args.vcf_file:
        parser.error("❌ 必须指定VCF文件 | VCF file must be specified")

    if not args.phenotype_file:
        parser.error("❌ 必须指定表型文件 | Phenotype file must be specified")

    try:
        # 创建分析器 | Create analyzer
        analyzer = GWASAnalyzer(
            vcf_file=args.vcf_file,
            phenotype_file=args.phenotype_file,
            output_prefix=args.output_prefix,
            manhattan_plot=args.manhattan_plot,
            qq_plot=args.qq_plot,
            significant_snps=args.significant_snps,
            excel_report=args.excel_report,
            html_report=args.html_report,
            plink_path=args.plink_path,
            bcftools_path=args.bcftools_path,
            admixture_path=args.admixture_path,
            emmax_path=args.emmax_path,
            maf_threshold=args.maf_threshold,
            missing_threshold=args.missing_threshold,
            depth_min=args.depth_min,
            depth_max=args.depth_max,
            qual_min=args.qual_min,
            admixture_k_range=tuple(args.admixture_k_range),
            pca_components=args.pca_components,
            ld_window=args.ld_window,
            ld_step=args.ld_step,
            ld_r2=args.ld_r2,
            emmax_precision=args.emmax_precision,
            p_value_threshold=args.p_value_threshold,
            working_dir=args.working_dir,
        )

        # 运行分析 | Run analysis
        analyzer.run_analysis()

    except Exception as e:
        print(f"❌ 错误 | Error: {e}", file=sys.stderr)
        sys.exit(1)
