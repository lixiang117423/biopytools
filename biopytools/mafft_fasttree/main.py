"""
🌳 系统发育树构建主程序模块 | Phylogenetic Tree Builder Main Module
"""

import argparse
import sys
from .config import PhyloConfig
from .utils import PhyloLogger, CommandRunner, check_dependencies
from .sequence_processor import SequenceProcessor
from .alignment import MAFFTAligner
from .tree_builder import FastTreeBuilder

class PhyloTreeBuilder:
    """系统发育树构建主类 | Main Phylogenetic Tree Builder Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = PhyloConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = PhyloLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.seq_processor = SequenceProcessor(self.config, self.logger)
        self.aligner = MAFFTAligner(self.config, self.logger, self.cmd_runner)
        self.tree_builder = FastTreeBuilder(self.config, self.logger, self.cmd_runner)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_pipeline(self):
        """运行完整的系统发育树构建流程 | Run complete phylogenetic tree building pipeline"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("🌳 开始系统发育树构建流程 | Starting phylogenetic tree building pipeline")
            self.logger.info("=" * 60)
            
            # 1. 检查依赖
            self.logger.info("📋 步骤 1/5: 检查依赖软件 | Step 1/5: Checking dependencies")
            self.check_dependencies()
            
            # 2. 检测序列类型
            self.logger.info("📋 步骤 2/5: 检测序列类型 | Step 2/5: Detecting sequence type")
            if self.config.seq_type:
                seq_type = self.config.seq_type
                self.logger.info(f"✅ 使用用户指定的序列类型 | Using user-specified sequence type: {seq_type}")
            else:
                seq_type = self.seq_processor.detect_sequence_type(self.config.input_file)
            
            # 3. 清理序列和处理ID
            self.logger.info("📋 步骤 3/5: 清理序列和处理ID | Step 3/5: Cleaning sequences and processing IDs")
            id_mapping = self.seq_processor.clean_sequences(
                self.config.input_file,
                self.config.cleaned_file
            )
            
            # 4. 运行MAFFT多序列比对
            self.logger.info("📋 步骤 4/5: MAFFT多序列比对 | Step 4/5: MAFFT multiple sequence alignment")
            if not self.aligner.run_alignment(self.config.cleaned_file, self.config.mafft_file):
                raise RuntimeError("MAFFT比对失败 | MAFFT alignment failed")
            
            # 5. 运行FastTree构建系统发育树
            self.logger.info("📋 步骤 5/5: FastTree构建系统发育树 | Step 5/5: FastTree phylogenetic tree construction")
            if not self.tree_builder.build_tree(self.config.mafft_file, self.config.tree_file, seq_type):
                raise RuntimeError("FastTree构建失败 | FastTree construction failed")
            
            # 完成
            self.logger.info("=" * 60)
            self.logger.info("✅ 系统发育树构建流程完成 | Phylogenetic tree building pipeline completed")
            self.logger.info("=" * 60)
            
            # 输出结果文件信息
            self.logger.info("📁 输出文件 | Output files:")
            self.logger.info(f"  - 清理后序列 | Cleaned sequences: {self.config.cleaned_file}")
            self.logger.info(f"  - MAFFT比对结果 | MAFFT alignment: {self.config.mafft_file}")
            self.logger.info(f"  - FastTree系统发育树 | FastTree phylogenetic tree: {self.config.tree_file}")
            if id_mapping:
                self.logger.info(f"  - ID映射文件 | ID mapping file: {self.config.id_mapping_file}")
            self.logger.info(f"  - 日志文件 | Log file: {self.logger_manager.log_file}")
            self.logger.info(f"\n📂 输出目录 | Output directory: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"❌ 流程执行失败 | Pipeline execution failed: {e}")
            import traceback
            self.logger.debug(traceback.format_exc())
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🌳 系统发育树构建脚本 (模块化版本) | Phylogenetic Tree Builder Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s -i sequences.fa -o phylo_results
  %(prog)s -i proteins.fa -o results -t 64
  %(prog)s -i sequences.fa -o results --mafft-params "--maxiterate 1000" --fasttree-params "-gamma"
  %(prog)s -i nucleotides.fa -o results --seq-type nucleotide
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', required=True, 
                       help='🧬 输入序列文件 (FASTA格式) | Input sequence file (FASTA format)')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./phylo_output', 
                       help='📁 输出目录 | Output directory')
    
    # 运行参数 | Run parameters
    parser.add_argument('-t', '--threads', type=int, default=88, 
                       help='🧵 线程数 | Number of threads')
    
    parser.add_argument('--seq-type', choices=['protein', 'nucleotide'],
                       help='🔬 序列类型 (不指定则自动检测) | Sequence type (auto-detect if not specified)')
    
    # MAFFT参数 | MAFFT parameters
    parser.add_argument('--mafft-params', default='--auto',
                       help='🧩 MAFFT额外参数 | Additional MAFFT parameters')
    
    # FastTree参数 | FastTree parameters  
    parser.add_argument('--fasttree-params', default='',
                       help='🌳 FastTree额外参数 | Additional FastTree parameters')
    
    # 工具路径 | Tool paths
    parser.add_argument('--mafft-path', default='mafft', 
                       help='🔧 MAFFT软件路径 | MAFFT software path')
    parser.add_argument('--fasttree-path', default='fasttree', 
                       help='🔧 FastTree软件路径 | FastTree software path')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    builder = PhyloTreeBuilder(
        input_file=args.input,
        output_dir=args.output,
        seq_type=args.seq_type,
        threads=args.threads,
        mafft_params=args.mafft_params,
        fasttree_params=args.fasttree_params,
        mafft_path=args.mafft_path,
        fasttree_path=args.fasttree_path
    )
    
    builder.run_pipeline()

if __name__ == "__main__":
    main()
