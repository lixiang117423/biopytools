"""
序列提取工具主程序模块 🎯 | Sequence Extraction Tool Main Module
"""

import argparse
import sys
from .config import ExtractorConfig
from .utils import ExtractorLogger, CommandRunner, check_dependencies
from .sequence_processor import SequenceProcessor

class SequenceExtractor:
    """序列提取器主类 🧬 | Main Sequence Extractor Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = ExtractorConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = ExtractorLogger(self.config.output_file, verbose=self.config.verbose)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)
        
        # 初始化序列处理器 | Initialize sequence processor
        self.sequence_processor = SequenceProcessor(self.config, self.logger, self.cmd_runner)
    
    def check_dependencies(self):
        """检查依赖软件 🔍 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_extraction(self):
        """运行完整的序列提取流程 🚀 | Run complete sequence extraction pipeline"""
        try:
            self.logger.info("🧬 开始序列提取流程 | Starting sequence extraction pipeline")
            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"序列文件 | Sequence file: {self.config.sequence_file}")
            self.logger.info(f"区域文件 | Regions file: {self.config.regions_file}")
            self.logger.info(f"输出文件 | Output file: {self.config.output_file}")
            self.logger.info(f"序列类型 | Sequence type: {self.config.sequence_type}")
            self.logger.info(f"线程数 | Threads: {self.config.threads}")
            self.logger.info(f"{'=' * 60}")
            
            # 检查依赖 | Check dependencies
            if not self.check_dependencies():
                self.logger.error("❌ 依赖检查失败 | Dependency check failed")
                return False
            
            # 提取序列 | Extract sequences
            self.logger.info("🎯 开始序列提取 | Starting sequence extraction")
            if not self.sequence_processor.extract_sequences():
                self.logger.error("❌ 序列提取失败 | Sequence extraction failed")
                return False
            
            self.logger.info("🎉 序列提取流程完成 | Sequence extraction pipeline completed")
            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"结果保存在 | Results saved in: {self.config.output_file}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 提取流程在执行过程中意外终止 | Extraction pipeline terminated unexpectedly: {e}")
            return False

def main():
    """主函数 🎯 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 序列提取工具 (模块化版本) | Sequence Extraction Tool (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
🌟 示例 | Examples:
  %(prog)s -s genome.fasta -r regions.bed -o extracted.fasta --type dna
  %(prog)s -s proteins.fasta -r regions.bed -o extracted.fasta --type protein
  %(prog)s -s genome.fasta -r regions.bed -o extracted.fasta --type dna --threads 16
  %(prog)s -s genome.fasta -r regions.bed -o extracted.fasta --type dna --reverse-complement
  %(prog)s -s genome.fasta -r regions.bed -o extracted.fasta --type dna --translate
  
🔄 区域文件格式 | Regions file format:
  染色体名  起始位置  终止位置  [链信息]
  chr1      1000      2000      +
  chr2      5000      6000      -
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-s', '--sequence-file', required=True, 
                       help='🧬 输入序列文件 (FASTA格式) | Input sequence file (FASTA format)')
    parser.add_argument('-r', '--regions-file', required=True,
                       help='📊 区域文件 (类似BED格式: 染色体 起始 终止 [链]) | Regions file (BED-like format: chromosome start end [strand])')
    parser.add_argument('-o', '--output-file', required=True,
                       help='💾 输出序列文件 | Output sequence file')
    
    # 序列参数 | Sequence parameters
    parser.add_argument('--type', '--sequence-type', choices=['dna', 'protein'], default='dna',
                       help='🔬 序列类型 | Sequence type')
    
    # 处理参数 | Processing parameters
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='🚀 线程数 | Number of threads')
    parser.add_argument('--merge-output', action='store_true', default=True,
                       help='📦 合并输出到一个文件 | Merge output to one file')
    parser.add_argument('--separate-output', action='store_true',
                       help='📂 分别输出到多个文件 | Output to separate files')
    parser.add_argument('--no-headers', action='store_true',
                       help='🏷️ 不包含区域信息在序列名中 | Do not include region info in sequence names')
    
    # DNA特有参数 | DNA-specific parameters
    parser.add_argument('--reverse-complement', action='store_true',
                       help='🔄 反向互补DNA序列 (会覆盖第4列的链信息) | Reverse complement DNA sequences (overrides 4th column strand info)')
    parser.add_argument('--translate', action='store_true',
                       help='🔄 将DNA翻译为蛋白质 | Translate DNA to protein')
    
    # 输出参数 | Output parameters
    parser.add_argument('--line-width', type=int, default=80,
                       help='📏 FASTA序列每行字符数 | Characters per line in FASTA sequence')
    parser.add_argument('-v', '--verbose', action='store_true', default=True,
                       help='📝 详细输出 | Verbose output')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='🔇 安静模式 | Quiet mode')
    
    # 工具路径 | Tool paths
    parser.add_argument('--samtools-path', default='samtools',
                       help='🛠️ samtools软件路径 | samtools software path')
    
    args = parser.parse_args()
    
    # 处理互斥参数 | Handle mutually exclusive parameters
    if args.separate_output:
        args.merge_output = False
    
    if args.quiet:
        args.verbose = False
    
    # 创建提取器并运行 | Create extractor and run
    extractor = SequenceExtractor(
        sequence_file=args.sequence_file,
        regions_file=args.regions_file,
        output_file=args.output_file,
        sequence_type=args.type,
        threads=args.threads,
        merge_output=args.merge_output,
        include_headers=not args.no_headers,
        reverse_complement=args.reverse_complement,
        translate_dna=args.translate,
        line_width=args.line_width,
        verbose=args.verbose,
        samtools_path=args.samtools_path
    )
    
    success = extractor.run_extraction()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
