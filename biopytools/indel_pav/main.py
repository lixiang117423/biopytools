"""
INDEL PAV分析主程序模块 | INDEL PAV Analysis Main Module
"""

import argparse
import sys
from .config import PAVConfig
from .utils import PAVLogger, check_dependencies
from .pav_analysis import PAVAnalyzer, PAVFilter
from .results import ResultsWriter

class IndelPAVAnalyzer:
    """INDEL PAV分析主类 | Main INDEL PAV Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = PAVConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = PAVLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化各个处理器 | Initialize processors
        self.pav_analyzer = PAVAnalyzer(self.config, self.logger)
        self.pav_filter = PAVFilter(self.config, self.logger)
        self.results_writer = ResultsWriter(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的PAV分析流程 | Run complete PAV analysis pipeline"""
        try:
            self.logger.info("🚀 开始INDEL PAV分析流程 | Starting INDEL PAV analysis pipeline")
            self.logger.info(f"📁 VCF文件 | VCF file: {self.config.vcf_file}")
            self.logger.info(f"📁 输出文件 | Output file: {self.config.output_file}")
            
            # 检查依赖 | Check dependencies
            self.check_dependencies()
            
            # 步骤1: PAV分析 | Step 1: PAV analysis
            self.logger.info(f"{'=' * 60}")
            self.logger.info("🧬 步骤1: INDEL PAV分析 | Step 1: INDEL PAV analysis")
            self.logger.info(f"{'=' * 60}")
            
            samples, results = self.pav_analyzer.run_pav_analysis()
            
            if not results:
                self.logger.warning("⚠️  没有找到合格的INDEL，分析结束 | No qualified INDELs found, analysis terminated")
                return
            
            # 步骤2: 结果过滤 | Step 2: Result filtering
            self.logger.info(f"{'=' * 60}")
            self.logger.info("🔍 步骤2: 结果过滤 | Step 2: Result filtering")
            self.logger.info(f"{'=' * 60}")
            
            filtered_results = self.pav_filter.filter_results(results)
            
            # 步骤3: 结果输出 | Step 3: Result output
            self.logger.info(f"{'=' * 60}")
            self.logger.info("💾 步骤3: 结果输出 | Step 3: Result output")
            self.logger.info(f"{'=' * 60}")
            
            self.results_writer.write_pav_results(samples, filtered_results)
            self.results_writer.write_summary_report(samples, filtered_results)
            
            self.logger.info(f"{'=' * 60}")
            self.logger.info("🎉 INDEL PAV分析完成！| INDEL PAV analysis completed!")
            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"📊 结果保存在 | Results saved in: {self.config.output_file}")
            
        except Exception as e:
            self.logger.error(f"❌ 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 INDEL PAV分析脚本 (模块化版本) | INDEL PAV Analysis Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s -v variants.vcf -o indel_pav.txt
  %(prog)s -v data.vcf.gz -o pav_results.txt -t 16 --min-length 5
  %(prog)s -v variants.vcf -o pav.txt --max-length 100 --compress
  %(prog)s -v filtered.vcf.gz -o results.txt -q 30 -d 10 --max-missing 0.5
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-v', '--vcf', required=True, 
                       help='📁 输入VCF文件路径 (支持压缩和未压缩) | Input VCF file path (supports compressed and uncompressed)')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./indel_pav.txt', 
                       help='📁 输出文件路径 | Output file path')
    
    # 分析参数 | Analysis parameters
    parser.add_argument('-t', '--threads', type=int, default=88, 
                       help='🚀 线程数 | Number of threads')
    parser.add_argument('--min-length', type=int, default=1, 
                       help='📏 最小INDEL长度 (bp) | Minimum INDEL length (bp)')
    parser.add_argument('--max-length', type=int, 
                       help='📏 最大INDEL长度 (bp) | Maximum INDEL length (bp)')
    
    # 过滤参数 | Filtering parameters
    parser.add_argument('-q', '--min-quality', type=float, default=20.0, 
                       help='⭐ 最小质量分数 | Minimum quality score')
    parser.add_argument('-d', '--min-depth', type=int, default=5, 
                       help='📊 最小深度 | Minimum depth')
    parser.add_argument('--max-missing', type=float, default=0.8, 
                       help='❓ 最大缺失率 (0-1) | Maximum missing rate (0-1)')
    
    # 输出参数 | Output parameters
    parser.add_argument('--include-complex', action='store_true',
                       help='🔀 包含复杂变异 | Include complex variants')
    parser.add_argument('--compress', action='store_true',
                       help='🗜️  压缩输出文件 | Compress output file')
    
    # 工具路径 | Tool paths
    parser.add_argument('--bcftools-path', default='bcftools', 
                       help='🔧 BCFtools软件路径 | BCFtools software path')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = IndelPAVAnalyzer(
        vcf_file=args.vcf,
        output_file=args.output,
        threads=args.threads,
        min_length=args.min_length,
        max_length=args.max_length,
        min_quality=args.min_quality,
        min_depth=args.min_depth,
        max_missing_rate=args.max_missing,
        include_complex=args.include_complex,
        compress_output=args.compress,
        bcftools_path=args.bcftools_path
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
