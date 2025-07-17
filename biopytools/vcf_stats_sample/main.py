"""
VCF基因型统计主程序模块 | VCF Genotype Statistics Main Module
"""

import argparse
import sys
from .config import VCFStatsConfig
from .utils import VCFStatsLogger
from .data_processing import VCFProcessor, StatisticsCalculator
from .results import ResultsExporter, SummaryGenerator

class VCFStatsAnalyzer:
    """VCF基因型统计主类 | Main VCF Genotype Statistics Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = VCFStatsConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = VCFStatsLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化各个处理器 | Initialize processors
        self.vcf_processor = VCFProcessor(self.config, self.logger)
        self.stats_calculator = StatisticsCalculator(self.config, self.logger)
        self.results_exporter = ResultsExporter(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的VCF基因型统计分析流程 | Run complete VCF genotype statistics analysis pipeline"""
        try:
            self.logger.info("="*20 + " 开始VCF基因型统计分析 | Starting VCF Genotype Statistics Analysis " + "="*20)
            self.logger.info(f"VCF文件 | VCF file: {self.config.vcf_file}")
            self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
            
            if self.config.min_depth > 0:
                self.logger.info(f"最小深度过滤 | Minimum depth filter: {self.config.min_depth}")
            if self.config.min_qual > 0:
                self.logger.info(f"最小质量过滤 | Minimum quality filter: {self.config.min_qual}")
            if self.config.exclude_missing:
                self.logger.info("将排除缺失基因型 | Will exclude missing genotypes")
            
            # 步骤1: 处理VCF文件 | Step 1: Process VCF file
            self.logger.info("\n步骤 1/3: 处理VCF文件 | Step 1/3: Processing VCF file")
            if not self.vcf_processor.process_vcf():
                self.logger.error("VCF文件处理失败 | VCF file processing failed")
                sys.exit(1)
            
            # 步骤2: 计算统计结果 | Step 2: Calculate statistics
            self.logger.info("\n步骤 2/3: 计算统计结果 | Step 2/3: Calculating statistics")
            stats_results = self.stats_calculator.calculate_rates(self.vcf_processor.sample_stats)
            
            if not stats_results:
                self.logger.error("统计计算失败 | Statistics calculation failed")
                sys.exit(1)
            
            self.logger.info(f"成功计算了 {len(stats_results)} 个样本的统计结果 | Successfully calculated statistics for {len(stats_results)} samples")
            
            # 步骤3: 导出结果 | Step 3: Export results
            self.logger.info("\n步骤 3/3: 导出结果 | Step 3/3: Exporting results")
            self.results_exporter.export_summary_statistics(stats_results, self.vcf_processor.detailed_stats)
            self.results_exporter.export_detailed_statistics(self.vcf_processor.detailed_stats)
            self.results_exporter.export_per_sample_files(stats_results, self.vcf_processor.detailed_stats)
            
            # 生成分析总结 | Generate analysis summary
            self.summary_generator.generate_analysis_summary(stats_results)
            
            self.logger.info("\n" + "="*20 + " VCF基因型统计分析完成 | VCF Genotype Statistics Analysis Completed " + "="*20)
            self.logger.info(f"结果文件已保存到 | Results saved to: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="VCF基因型统计分析脚本 | VCF Genotype Statistics Analysis Script\n"
                   "支持短参数和长参数格式 | Supports both short and long parameter formats",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
        参数对照表 | Parameter Reference:
        -v/--vcf          输入VCF文件 | Input VCF file
        -o/--output       输出目录 | Output directory  
        -d/--min-depth    最小深度 | Minimum depth
        -q/--min-qual     最小质量 | Minimum quality
        -e/--exclude-missing  排除缺失 | Exclude missing
        -D/--no-detailed  无详细输出 | No detailed output
        -S/--no-summary   无汇总输出 | No summary output
        
        示例用法 | Example Usage:
        
        # 基本分析 | Basic analysis
        python run_vcf_stats.py -v variants.vcf -o vcf_stats_output
        
        # 应用质量和深度过滤 | Apply quality and depth filters
        python run_vcf_stats.py -v variants.vcf.gz -o results -d 10 -q 30.0
        
        # 排除缺失基因型，仅输出汇总统计 | Exclude missing genotypes, summary only
        python run_vcf_stats.py -v variants.vcf -o results -e -D
        
        # 长参数格式 | Long parameter format
        python run_vcf_stats.py --vcf variants.vcf --output results \\
            --min-depth 10 --min-qual 30.0 --exclude-missing --no-detailed
            
        支持的基因型格式 | Supported genotype formats:
        - 未定相: 0/0, 0/1, 1/1, ./. | Unphased: 0/0, 0/1, 1/1, ./.
        - 已定相: 0|0, 0|1, 1|1, .|. | Phased: 0|0, 0|1, 1|1, .|.
        - 多等位基因: 0/2, 1/2等 | Multi-allelic: 0/2, 1/2, etc.
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument("-v", "--vcf", required=True, 
                       help="输入VCF文件路径 (支持.gz压缩格式) | Input VCF file path (supports .gz compressed format)")
    
    # 可选参数 | Optional arguments
    parser.add_argument("-o", "--output", default="vcf_stats_output", 
                       help="输出目录 | Output directory")
    parser.add_argument("-d", "--min-depth", type=int, default=0, 
                       help="最小深度过滤阈值 (0表示不过滤) | Minimum depth filter threshold (0 = no filter)")
    parser.add_argument("-q", "--min-qual", type=float, default=0.0, 
                       help="最小质量分数过滤阈值 (0.0表示不过滤) | Minimum quality score filter threshold (0.0 = no filter)")
    parser.add_argument("-e", "--exclude-missing", action="store_true", 
                       help="排除缺失基因型 (./..) 的统计 | Exclude missing genotypes (./..) from statistics")
    parser.add_argument("-D", "--no-detailed", action="store_true", 
                       help="不输出详细统计结果 | Do not output detailed statistics")
    parser.add_argument("-S", "--no-summary", action="store_true", 
                       help="不输出汇总统计结果 | Do not output summary statistics")
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = VCFStatsAnalyzer(
        vcf_file=args.vcf,
        output_dir=args.output,
        min_depth=args.min_depth,
        min_qual=args.min_qual,
        exclude_missing=args.exclude_missing,
        output_detailed=not args.no_detailed,
        output_summary=not args.no_summary
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
