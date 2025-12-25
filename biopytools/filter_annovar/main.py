"""
主程序模块 | Main Module
"""

import argparse
import sys
from .config import FilterConfig
from .logger import FilterLogger
from .processor import GeneProcessor

class FilterAnnovarAnalyzer:
    """Filter ANNOVAR分析主类 | Main Filter ANNOVAR Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = FilterConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = FilterLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化处理器 | Initialize processor
        self.processor = GeneProcessor(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的分析流程 | Run complete analysis pipeline"""
        try:
            self.logger.info("=" * 100)
            self.logger.info("🧬 Filter ANNOVAR - 基因区域变异提取工具")
            self.logger.info("=" * 100)
            
            # 处理所有基因 | Process all genes
            self.processor.process_all_genes()
            
            self.logger.info(f"\n{'=' * 100}")
            self.logger.info(f"✅ 分析完成 | Analysis completed successfully")
            self.logger.info(f"{'=' * 100}")
            
        except Exception as e:
            self.logger.error(f"❌ 分析过程中发生错误 | Error during analysis: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 Filter ANNOVAR - 基因区域变异提取工具 | Gene Region Variant Extraction Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:
  # 单个基因分析 | Single gene analysis
  %(prog)s -g gene-PHYSODRAFT_288440 --gff genome.gff \\
      --exonic variants.exonic_variant_function \\
      --all variants.variant_function -o output_dir
  
  # 批量基因分析 | Batch gene analysis  
  %(prog)s --gene-list genes.txt --gff genome.gff \\
      --exonic variants.exonic_variant_function \\
      --all variants.variant_function -o results -f txt
  
  # 自定义扩展范围和线程数 | Custom extension and threads
  %(prog)s -g gene-ID --gff genome.gff \\
      --exonic exonic.txt --all all.txt -e 10000 -t 44
        """
    )
    
    # 必需参数 | Required arguments
    gene_group = parser.add_mutually_exclusive_group(required=True)
    gene_group.add_argument('-g', '--gene-id',
                           help='🎯 单个基因ID | Single gene ID')
    gene_group.add_argument('-G', '--gene-list',
                           help='📋 基因ID列表文件(每行一个基因) | Gene ID list file (one gene per line)')
    
    parser.add_argument('--gff', required=True,
                       help='📄 GFF注释文件路径 | GFF annotation file path')
    parser.add_argument('--exonic', required=True,
                       help='📄 ANNOVAR外显子注释文件(.exonic_variant_function) | ANNOVAR exonic annotation file')
    parser.add_argument('--all', required=True, dest='all_variant',
                       help='📄 ANNOVAR所有变异注释文件(.variant_function) | ANNOVAR all variants annotation file')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./filter_output',
                       help='📁 输出目录 | Output directory (default: ./filter_output)')
    parser.add_argument('-f', '--format', choices=['excel', 'txt'], default='excel',
                       help='📊 输出格式 | Output format (default: excel)')
    parser.add_argument('-e', '--extend', type=int, default=5000,
                       help='📏 上下游扩展范围(bp) | Upstream/downstream extension range in bp (default: 5000)')
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='🧵 线程数 | Number of threads (default: 88)')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = FilterAnnovarAnalyzer(
        gff_file=args.gff,
        exonic_file=args.exonic,
        all_variant_file=args.all_variant,
        gene_id=args.gene_id,
        gene_list_file=args.gene_list,
        output_dir=args.output,
        output_format=args.format,
        extend_bp=args.extend,
        threads=args.threads
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
