"""
主程序入口 | Main Program Entry
"""

import argparse
import sys
import os
from pathlib import Path

# 添加模块路径 | Add module path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from config import LDHeatmapConfig
from utils import LDLogger, check_dependencies
from analyzer import LDHeatmapAnalyzer

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
    
    # 创建配置 | Create configuration
    try:
        config = LDHeatmapConfig(
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
    except ValueError as e:
        print(f"配置错误 | Configuration error: {e}")
        sys.exit(1)
    
    # 设置日志 | Setup logging
    logger_manager = LDLogger(config.output_dir)
    logger = logger_manager.get_logger()
    
    # 检查依赖 | Check dependencies
    try:
        check_dependencies(logger)
    except RuntimeError as e:
        logger.error(f"依赖检查失败 | Dependency check failed: {e}")
        sys.exit(1)
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = LDHeatmapAnalyzer(config, logger)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
