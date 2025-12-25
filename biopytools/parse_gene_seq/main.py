"""
基因序列提取主程序模块 🧬 | Gene Sequence Extraction Main Module
"""

import argparse
import sys
from pathlib import Path
from .config import ExtractionConfig
from .utils import ExtractionLogger
from .sequence_extractor import GenomeLoader, GFFParser, GeneExtractor, FASTAWriter

class GeneSequenceExtractor:
    """基因序列提取主类 🧬 | Main Gene Sequence Extractor Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = ExtractionConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        output_dir = Path(self.config.output_file).parent
        self.logger_manager = ExtractionLogger(output_dir)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化各个处理器 | Initialize processors
        self.genome_loader = GenomeLoader(self.logger)
        self.gff_parser = GFFParser(self.logger)
        self.gene_extractor = GeneExtractor(self.logger)
        self.fasta_writer = FASTAWriter(self.logger)
    
    def run_extraction(self):
        """运行完整的基因序列提取流程 🚀 | Run complete gene sequence extraction pipeline"""
        try:
            self.logger.info("🧬 基因序列提取工具启动 | Gene sequence extraction tool started")
            self.logger.info(f"{'=' * 60}")
            
            if self.config.verbose:
                self.logger.info(f"📋 参数信息 | Parameter info:")
                self.logger.info(f"   基因组文件 | Genome file: {self.config.genome_file}")
                self.logger.info(f"   GFF文件 | GFF file: {self.config.gff_file}")
                self.logger.info(f"   输出文件 | Output file: {self.config.output_file}")
                self.logger.info(f"   特征类型 | Feature type: {self.config.feature_type}")
                self.logger.info(f"   最小长度 | Minimum length: {self.config.min_length}")
                self.logger.info(f"   线程数 | Threads: {self.config.threads}")
                self.logger.info("")
            
            # 步骤1: 加载基因组 | Step 1: Load genome
            genome_seqs = self.genome_loader.load_genome(self.config.genome_file)
            
            # 步骤2: 解析GFF文件 | Step 2: Parse GFF file
            features = self.gff_parser.parse_gff(self.config.gff_file, self.config.feature_type)
            
            if not features:
                self.logger.error(f"❌ 没有找到任何 '{self.config.feature_type}' 类型的特征 | No '{self.config.feature_type}' features found")
                sys.exit(1)
            
            # 步骤3: 提取基因序列 | Step 3: Extract gene sequences
            extracted_genes = self.gene_extractor.extract_genes(
                genome_seqs, features, 
                self.config.min_length, self.config.verbose
            )
            
            if not extracted_genes:
                self.logger.error("❌ 没有提取到任何基因序列 | No gene sequences extracted")
                sys.exit(1)
            
            # 步骤4: 写入输出文件 | Step 4: Write output file
            self.fasta_writer.write_fasta(extracted_genes, self.config.output_file, self.config.line_width)
            
            self.logger.info(f"{'=' * 60}")
            self.logger.info("🎉 基因序列提取完成 | Gene sequence extraction completed!")
            self.logger.info(f"📊 总计提取 {len(extracted_genes)} 个基因序列 | Total extracted {len(extracted_genes)} gene sequences")
            
        except KeyboardInterrupt:
            self.logger.warning("\n⚠️ 用户中断操作 | User interrupted operation")
            sys.exit(1)
        except Exception as e:
            self.logger.error(f"\n❌ 错误 | Error: {e}")
            sys.exit(1)

def main():
    """主函数 🚀 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 基因序列提取脚本 (模块化版本) | Gene Sequence Extraction Script (Modular Version)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
📚 使用示例 | Usage Examples:
  %(prog)s -g genome.fasta -f annotation.gff -o genes.fasta
  %(prog)s --genome genome.fa --gff genes.gff3 --output extracted_genes.fa
  %(prog)s -g genome.fasta -f annotation.gff -o genes.fasta --feature-type CDS --min-length 300 -v
  %(prog)s -g genome.fasta -f annotation.gff -o genes.fasta --threads 16 --line-width 80

📋 支持的文件格式 | Supported File Formats:
  - 基因组文件 | Genome files: FASTA格式 (.fasta, .fa, .fas)
  - GFF文件 | GFF files: GFF3格式 (.gff, .gff3)
  - 输出文件 | Output files: FASTA格式
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-g', '--genome', required=True,
                       help='🧬 基因组FASTA文件路径 | Genome FASTA file path')
    parser.add_argument('-f', '--gff', required=True,
                       help='📄 GFF注释文件路径 | GFF annotation file path')
    parser.add_argument('-o', '--output', required=True,
                       help='💾 输出FASTA文件路径 | Output FASTA file path')
    
    # 可选参数 | Optional arguments
    parser.add_argument('--feature-type', default='gene',
                       help='🎯 要提取的特征类型 | Feature type to extract (default: gene)')
    parser.add_argument('--min-length', type=int, default=0,
                       help='📏 最小基因长度过滤 | Minimum gene length filter (default: 0)')
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='🧵 线程数 | Number of threads (default: 88)')
    parser.add_argument('--line-width', type=int, default=60,
                       help='📐 FASTA序列行宽度 | FASTA sequence line width (default: 60)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='🔍 显示详细信息 | Show verbose output')
    
    args = parser.parse_args()
    
    # 创建提取器并运行 | Create extractor and run
    extractor = GeneSequenceExtractor(
        genome_file=args.genome,
        gff_file=args.gff,
        output_file=args.output,
        feature_type=args.feature_type,
        min_length=args.min_length,
        threads=args.threads,
        line_width=args.line_width,
        verbose=args.verbose
    )
    
    extractor.run_extraction()

if __name__ == "__main__":
    main()
