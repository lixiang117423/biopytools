"""
GenBank序列提取主程序模块 🚀 | GenBank Sequence Extraction Main Module
"""

import argparse
import sys
# 处理相对导入问题 | Handle relative import issues
try:
    from .config import ExtractorConfig
    from .utils import ExtractorLogger, FileManager, parallel_process_files
    from .sequence_processor import SequenceExtractor
    from .output_writer import SequenceWriter
    from .reporter import StatisticsReporter
    from .phylogenetic import PhylogeneticMatrixBuilder
except ImportError:
    # 如果相对导入失败，尝试绝对导入 | If relative import fails, try absolute import
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from config import ExtractorConfig
    from utils import ExtractorLogger, FileManager, parallel_process_files
    from sequence_processor import SequenceExtractor
    from output_writer import SequenceWriter
    from reporter import StatisticsReporter
    from phylogenetic import PhylogeneticMatrixBuilder

class GenBankExtractor:
    """GenBank序列提取主类 🧬 | Main GenBank Sequence Extractor Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = ExtractorConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = ExtractorLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化各个处理器 | Initialize processors
        self.file_manager = FileManager(self.config, self.logger)
        self.sequence_extractor = SequenceExtractor(self.config, self.logger)
        self.sequence_writer = SequenceWriter(self.config, self.logger)
        self.reporter = StatisticsReporter(self.config, self.logger)
        self.phylo_builder = PhylogeneticMatrixBuilder(self.config, self.logger)
    
    def run_extraction(self):
        """运行完整的序列提取流程 🚀 | Run complete sequence extraction pipeline"""
        try:
            self.logger.info("🚀 开始GenBank序列提取流程")
            self.logger.info(f"📂 输入目录: {self.config.input_dir}")
            self.logger.info(f"📁 输出目录: {self.config.output_dir}")
            self.logger.info(f"🔧 线程数: {self.config.threads}")
            self.logger.info("=" * 80)
            
            # 步骤1: 查找GenBank文件 | Step 1: Find GenBank files
            self.logger.info("📋 步骤1: 查找GenBank文件")
            genbank_files = self.file_manager.find_genbank_files()
            if not genbank_files:
                self.logger.error("❌ 未找到GenBank文件，退出程序")
                return False
            
            # 步骤2: 创建输出目录 | Step 2: Create output directories
            self.logger.info("🏗️ 步骤2: 创建输出目录结构")
            self.file_manager.create_output_directories()
            
            # 步骤3: 并行提取序列 | Step 3: Extract sequences in parallel
            self.logger.info(f"⚡ 步骤3: 并行提取序列 (使用{self.config.threads}个线程)")
            all_sample_data = parallel_process_files(
                genbank_files,
                self.sequence_extractor.process_single_file,
                self.config.threads,
                self.logger
            )
            
            if not all_sample_data:
                self.logger.error("❌ 序列提取失败，退出程序")
                return False
            
            # 步骤4: 写入输出文件 | Step 4: Write output files
            self.logger.info("💾 步骤4: 写入序列文件")
            self.sequence_writer.write_sequences(all_sample_data)
            
            # 步骤5: 生成统计报告 | Step 5: Generate statistics report
            self.logger.info("📊 步骤5: 生成统计报告")
            self.reporter.generate_comprehensive_report(all_sample_data, self.config.output_dir)
            
            # 步骤6: 创建系统发育矩阵 (可选) | Step 6: Create phylogenetic matrix (optional)
            if self.config.create_phylogenetic_matrix:
                self.logger.info("🌲 步骤6: 创建系统发育分析矩阵")
                self.phylo_builder.create_phylogenetic_matrix()
            
            # 输出总结 | Output summary
            total_genes_extracted = sum(data['genes_extracted'] for data in all_sample_data)
            unique_genes = len(set(gene for data in all_sample_data for gene in data['genes_list']))
            
            self.logger.info("=" * 80)
            self.logger.info("🎉 序列提取完成!")
            self.logger.info(f"📊 处理样本数: {len(all_sample_data)}")
            self.logger.info(f"🧬 总提取基因数: {total_genes_extracted}")
            self.logger.info(f"💎 不重复基因数: {unique_genes}")
            self.logger.info(f"📁 CDS文件输出到: {self.config.cds_dir}")
            self.logger.info(f"🧪 蛋白质文件输出到: {self.config.pep_dir}")
            self.logger.info("=" * 80)
            
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 序列提取流程在执行过程中意外终止: {e}")
            sys.exit(1)

def main():
    """主函数 🎯 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 GenBank序列提取工具 | GenBank Sequence Extraction Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
🌟 示例 | Examples:
  %(prog)s -i /path/to/genbank/files -o ./output
  %(prog)s --input ./gb_files --output ./results --threads 64
  %(prog)s -i ./genbank -o ./extract_results -t 32 --phylo
  %(prog)s -i /data/gb_files -o /results --min-length 20 --threads 88
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', '--input-dir', required=True, 
                       help='📂 输入GenBank文件目录 | Input GenBank files directory')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', '--output-dir', default='./genbank_output', 
                       help='📁 输出目录 | Output directory')
    
    # 处理参数 | Processing parameters
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='⚡ 并行线程数 | Number of parallel threads')
    parser.add_argument('--min-length', type=int, default=10,
                       help='🔬 最小蛋白质长度 (氨基酸) | Minimum protein length (amino acids)')
    
    # 输出选项 | Output options
    parser.add_argument('--phylo', '--create-phylogenetic-matrix', action='store_true',
                       help='🌲 创建系统发育分析矩阵 | Create phylogenetic analysis matrix')
    parser.add_argument('--no-sample-sep', action='store_true',
                       help='🚫 不按样品分离输出 | Do not separate output by sample')
    parser.add_argument('--no-gene-sep', action='store_true',
                       help='🚫 不按基因分离输出 | Do not separate output by gene')
    parser.add_argument('--keep-unknown', action='store_true',
                       help='📝 保留unknown基因 | Keep unknown genes')
    
    args = parser.parse_args()
    
    # 创建分析器 | Create analyzer
    try:
        extractor = GenBankExtractor(
            input_dir=args.input,
            output_dir=args.output,
            threads=args.threads,
            min_protein_length=args.min_length,
            create_phylogenetic_matrix=args.phylo,
            separate_by_sample=not args.no_sample_sep,
            separate_by_gene=not args.no_gene_sep,
            skip_unknown_genes=not args.keep_unknown
        )
        
        # 运行提取 | Run extraction
        success = extractor.run_extraction()
        sys.exit(0 if success else 1)
        
    except Exception as e:
        print(f"❌ 程序初始化失败: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
