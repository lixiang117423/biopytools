"""
VCF LD热图分析主程序模块 | VCF LD Heatmap Analysis Main Module
"""

import argparse
import sys
import os
from .config import LDHeatmapConfig
from .utils import LDLogger, check_dependencies
from .analyzer import LDHeatmapAnalyzer
from .results import SummaryGenerator

class LDHeatmapMain:
    """LD热图分析主类 | Main LD Heatmap Analysis Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = LDHeatmapConfig(**kwargs)
        
        # 初始化日志 | Initialize logging
        self.logger_manager = LDLogger(self.config.output_dir)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化各个处理器 | Initialize processors
        self.analyzer = LDHeatmapAnalyzer(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.logger)
    
    def run_analysis(self):
        """运行完整的LD分析流程 | Run complete LD analysis pipeline"""
        try:
            self.logger.info("开始VCF LD热图分析流程 | Starting VCF LD heatmap analysis pipeline")
            self.logger.info(f"VCF文件 | VCF file: {self.config.vcf_file}")
            self.logger.info(f"输出文件 | Output file: {self.config.output_file}")
            
            # 检查依赖 | Check dependencies
            self.check_dependencies()
            
            # 运行分析 | Run analysis
            self.analyzer.run_analysis()
            
            # 生成总结报告 | Generate summary report
            self.summary_generator.generate_summary_report()
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="VCF连锁不平衡热图生成器 (模块化版本) | VCF Linkage Disequilibrium Heatmap Generator (Modular Version)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
使用示例 | Examples:
  %(prog)s -i input.vcf -o output.png
  %(prog)s -i input.vcf -o output.pdf --region chr1:1000000-2000000
  %(prog)s -i input.vcf --maf 0.05 --max-snps 500 --figsize 12 10
  %(prog)s -i input.vcf -o heatmap.png --save-matrix ld_matrix.csv --triangle-only
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='输入VCF文件路径 | Input VCF file path')
    
    # 输出参数 | Output parameters
    parser.add_argument('-o', '--output', default='ld_heatmap.png',
                       help='输出图片文件路径 | Output image file path')
    
    # 过滤参数 | Filtering parameters
    parser.add_argument('--region', 
                       help='指定基因组区域 (格式: chr:start-end, 例如: chr1:1000000-2000000) | '
                            'Specify genomic region (format: chr:start-end, e.g.: chr1:1000000-2000000)')
    parser.add_argument('--maf', type=float, default=0.01,
                       help='最小等位基因频率阈值 | Minor allele frequency threshold')
    parser.add_argument('--max-snps', type=int, default=1000,
                       help='最大SNP数量限制 | Maximum SNP count limit')
    parser.add_argument('--samples', nargs='+',
                       help='指定样本名称列表 | Specify sample name list')
    
    # 图形参数 | Graphics parameters
    parser.add_argument('--figsize', nargs=2, type=int, default=[10, 8],
                       help='图形尺寸 (宽度 高度) | Figure size (width height)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='图像分辨率 | Image resolution')
    parser.add_argument('--colormap', default='RdYlGn_r',
                       help='颜色映射 | Color map')
    parser.add_argument('--title', 
                       help='图表标题 | Chart title')
    
    # 数据输出参数 | Data output parameters
    parser.add_argument('--save-matrix', 
                       help='保存LD矩阵到CSV文件 | Save LD matrix to CSV file')
    parser.add_argument('--ld-threshold', type=float, default=0.0,
                       help='LD阈值，低于此值显示为白色 | LD threshold, values below shown as white')
    
    # 其他参数 | Other parameters
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='显示详细信息 | Show verbose information')
    parser.add_argument('--triangle-only', action='store_true',
                       help='只显示上三角矩阵 | Show upper triangle only')
    
    args = parser.parse_args()
    
    # 检查输入文件 | Check input file
    if not os.path.exists(args.input):
        print(f"错误: 输入文件不存在 | Error: Input file does not exist - {args.input}")
        sys.exit(1)
    
    # 创建分析器并运行 | Create analyzer and run
    try:
        analyzer = LDHeatmapMain(
            vcf_file=args.input,
            output_file=args.output,
            region=args.region,
            maf=args.maf,
            max_snps=args.max_snps,
            samples=args.samples,
            figsize=args.figsize,
            dpi=args.dpi,
            colormap=args.colormap,
            title=args.title,
            save_matrix=args.save_matrix,
            ld_threshold=args.ld_threshold,
            triangle_only=args.triangle_only,
            verbose=args.verbose
        )
        
        analyzer.run_analysis()
        
    except ValueError as e:
        print(f"配置错误 | Configuration error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
