"""
VCF PCA分析主程序模块|VCF PCA Analysis Main Module
"""

import argparse
import sys
from .config import PCAConfig
from .utils import PCALogger, CommandRunner, check_dependencies
from .data_processing import VCFProcessor, QualityController
from .pca_analysis import PCAAnalyzer
from .visualization import PCAVisualizer
from .results import SummaryGenerator

class VCFPCAAnalyzer:
    """VCF PCA分析主类|Main VCF PCA Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = PCAConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志|Initialize logging
        self.logger_manager = PCALogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器|Initialize processors
        self.vcf_processor = VCFProcessor(self.config, self.logger, self.cmd_runner)
        self.quality_controller = QualityController(self.config, self.logger, self.cmd_runner)
        self.pca_analyzer = PCAAnalyzer(self.config, self.logger, self.cmd_runner)
        self.visualizer = PCAVisualizer(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件|Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的PCA分析流程|Run complete PCA analysis pipeline"""
        try:
            self.logger.info("开始VCF PCA分析流程|Starting VCF PCA analysis pipeline")
            self.logger.info(f"VCF文件|VCF file: {self.config.vcf_file}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            
            # 检查依赖|Check dependencies
            self.check_dependencies()
            
            # 检查可视化库|Check visualization libraries
            if self.config.plot:
                try:
                    import matplotlib.pyplot as plt
                    import seaborn as sns
                    self.logger.info("可视化库可用|Visualization libraries available")
                except ImportError:
                    self.logger.warning("可视化库不可用，将跳过图表生成|Visualization libraries unavailable, skipping plots")
                    self.config.plot = False
            
            # 步骤1: VCF转PLINK|Step 1: VCF to PLINK conversion
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤1: VCF转PLINK格式|Step 1: VCF to PLINK conversion")
            self.logger.info(f"{'=' * 60}")
            
            if not self.vcf_processor.vcf_to_plink():
                raise RuntimeError("VCF转换失败|VCF conversion failed")
            
            # 步骤2: 质量控制|Step 2: Quality control
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤2: 质量控制|Step 2: Quality control")
            self.logger.info(f"{'=' * 60}")
            
            if not self.quality_controller.apply_quality_filters():
                raise RuntimeError("质量控制失败|Quality control failed")
            
            # 步骤3: PCA分析|Step 3: PCA analysis
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤3: PCA分析|Step 3: PCA analysis")
            self.logger.info(f"{'=' * 60}")
            
            if not self.pca_analyzer.run_pca_analysis():
                raise RuntimeError("PCA分析失败|PCA analysis failed")
            
            # 步骤4: 可视化|Step 4: Visualization
            if self.config.plot:
                self.logger.info(f"{'=' * 60}")
                self.logger.info("步骤4: 生成可视化图表|Step 4: Generate visualization plots")
                self.logger.info(f"{'=' * 60}")
                
                self.visualizer.create_plots()
            
            # 生成总结报告|Generate summary report
            self.summary_generator.generate_summary_report()
            
            self.logger.info(f"{'=' * 60}")
            self.logger.info("PCA分析完成！| PCA analysis completed!")
            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"结果保存在|Results saved in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止|Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='VCF PCA分析脚本|VCF PCA Analysis Script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i variants.vcf -o pca_results
        """
    )

    # 必需参数|Required arguments
    # New standard: -i/--input, Legacy: -v/--vcf (for backward compatibility)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input', dest='vcf_file',
                            help='输入VCF文件路径|Input VCF file path')
    input_group.add_argument('-v', '--vcf', dest='vcf_file',
                            help=argparse.SUPPRESS)

    # 可选参数|Optional arguments
    parser.add_argument('-o', '--output', default='./pca_output',
                       help='输出目录|Output directory')
    parser.add_argument('-s', '--sample-info',
                       help='样本信息文件|Sample information file')

    # PCA参数|PCA parameters
    parser.add_argument('-c', '--components', type=int, default=10,
                       help='主成分数量|Number of principal components')

    # 质控参数|Quality control parameters
    # Note: -m is reserved for --memory in standards, using --maf without short option
    # Legacy -m is still supported but hidden
    parser.add_argument('-m', dest='maf_legacy', type=float, help=argparse.SUPPRESS)
    parser.add_argument('--maf', dest='maf', type=float, default=0.05,
                       help='最小等位基因频率阈值|Minor allele frequency threshold')
    parser.add_argument('--missing', type=float, default=0.1,
                       help='最大缺失率阈值|Maximum missing rate threshold')
    parser.add_argument('--hwe', type=float, default=1e-6,
                       help='Hardy-Weinberg平衡p值阈值|Hardy-Weinberg equilibrium p-value threshold')
    parser.add_argument('--skip-qc', action='store_true',
                       help='跳过质量控制过滤|Skip quality control filtering')

    # 可视化参数|Visualization parameters
    # Note: -p is reserved for --prefix in standards, using -P for --plot
    # Legacy -p is still supported but hidden
    parser.add_argument('-p', dest='plot_legacy', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('-P', '--plot', dest='plot', action='store_true',
                       help='生成PCA可视化图表|Generate PCA visualization plots')
    parser.add_argument('-g', '--group-column',
                       help='分组列名|Column name for grouping')

    # 工具路径|Tool paths
    parser.add_argument('--plink-path', default='plink',
                       help='PLINK软件路径|PLINK software path')
    parser.add_argument('--bcftools-path', default='bcftools',
                       help='BCFtools软件路径|BCFtools software path')

    args = parser.parse_args()

    # Handle backward compatibility for maf and plot
    maf_value = args.maf_legacy if args.maf_legacy is not None else args.maf
    plot_value = args.plot or args.plot_legacy

    # 创建分析器并运行|Create analyzer and run
    analyzer = VCFPCAAnalyzer(
        vcf_file=args.vcf_file,
        output_dir=args.output,
        sample_info_file=args.sample_info,
        components=args.components,
        maf=maf_value,
        missing_rate=args.missing,
        hwe_pvalue=args.hwe,
        skip_qc=args.skip_qc,
        plot=plot_value,
        group_column=args.group_column,
        plink_path=args.plink_path,
        bcftools_path=args.bcftools_path
    )

    analyzer.run_analysis()

if __name__ == "__main__":
    main()
