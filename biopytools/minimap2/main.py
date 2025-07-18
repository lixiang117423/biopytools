"""
Minimap2分析主程序模块 | Minimap2 Analysis Main Module
"""

import argparse
import sys
import os
from .config import Minimap2Config
from .utils import Minimap2Logger, CommandRunner
from .alignment import Minimap2Aligner
from .data_processing import PAFProcessor
from .sequence_extraction import SequenceExtractor
from .results import SummaryGenerator

class Minimap2Analyzer:
    """Minimap2分析主类 | Main Minimap2 Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = Minimap2Config(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = Minimap2Logger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.aligner = Minimap2Aligner(self.config, self.logger, self.cmd_runner)
        self.paf_processor = PAFProcessor(self.config, self.logger)
        self.sequence_extractor = SequenceExtractor(self.config, self.logger, self.cmd_runner)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的Minimap2分析流程 | Run complete Minimap2 analysis pipeline"""
        try:
            self.logger.info("开始Minimap2分析流程 | Starting Minimap2 analysis pipeline")
            self.logger.info(f"目标基因组 | Target genome: {self.config.target_genome}")
            self.logger.info(f"查询基因组 | Query genome: {self.config.query_genome}")
            self.logger.info("开始Minimap2分析流程 | Starting Minimap2 analysis pipeline")
            self.logger.info(f"目标基因组 | Target genome: {self.config.target_genome}")
            self.logger.info(f"查询基因组 | Query genome: {self.config.query_genome}")
            self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
            self.logger.info(f"预设参数 | Preset: {self.config.preset}")
            self.logger.info(f"tp类型过滤 | tp type filter: {self.config.tp_type}")
            self.logger.info(f"线程数 | Threads: {self.config.threads}")
            self.logger.info(f"最小匹配长度 | Min match length: {self.config.min_match_length}")
            self.logger.info(f"最小未比对长度 | Min unmapped length: {self.config.min_unmapped_length}")
            
            # 步骤1: 运行Minimap2比对 | Step 1: Run Minimap2 alignment
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤1: 运行Minimap2比对 | Step 1: Running Minimap2 alignment")
            self.logger.info("="*60)
            
            if not self.aligner.run_alignment():
                self.logger.error("Minimap2比对失败 | Minimap2 alignment failed")
                sys.exit(1)
            
            # 步骤2: 处理PAF文件 | Step 2: Process PAF file
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤2: 处理PAF文件 | Step 2: Processing PAF file")
            self.logger.info("="*60)
            
            # 加载PAF文件 | Load PAF file
            df_paf = self.paf_processor.load_paf_file()
            if df_paf is None:
                self.logger.error("加载PAF文件失败 | Failed to load PAF file")
                sys.exit(1)
            
            # 筛选PAF数据 | Filter PAF data
            df_filtered = self.paf_processor.filter_paf_data(df_paf)
            if df_filtered.empty:
                self.logger.warning("筛选后没有找到有效的比对结果 | No valid alignments found after filtering")
                unmapped_regions = []
            else:
                # 查找未比对区间 | Find unmapped regions
                unmapped_regions = self.paf_processor.find_unmapped_regions(df_filtered)
            
            # 步骤3: 保存BED文件 | Step 3: Save BED file
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤3: 保存BED文件 | Step 3: Saving BED file")
            self.logger.info("="*60)
            
            if not self.paf_processor.save_bed_file(unmapped_regions):
                self.logger.error("保存BED文件失败 | Failed to save BED file")
                sys.exit(1)
            
            # 步骤4: 提取未比对序列 | Step 4: Extract unmapped sequences
            self.logger.info("\n" + "="*60)
            self.logger.info("步骤4: 提取未比对序列 | Step 4: Extracting unmapped sequences")
            self.logger.info("="*60)
            
            if not self.sequence_extractor.extract_unmapped_sequences(unmapped_regions):
                self.logger.error("提取未比对序列失败 | Failed to extract unmapped sequences")
                sys.exit(1)
            
            # 生成总结报告 | Generate summary report
            self.summary_generator.generate_summary_report(unmapped_regions)
            
            self.logger.info("\n" + "="*60)
            self.logger.info("分析完成！| Analysis completed!")
            self.logger.info("="*60)
            self.logger.info(f"输出文件位于 | Output files in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='Minimap2全基因组比对和未比对区间提取工具 | Minimap2 Whole Genome Alignment and Unmapped Region Extraction Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-t', '--target', required=True, 
                       help='目标基因组文件路径 | Target genome file path')
    parser.add_argument('-q', '--query', required=True, 
                       help='查询基因组文件路径 | Query genome file path')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output-dir', default='./minimap2_output',
                       help='输出目录 | Output directory')
    parser.add_argument('-x', '--preset', default='asm5', 
                       choices=['asm5', 'asm10', 'asm20', 'map-ont', 'map-pb'],
                       help='Minimap2预设参数 | Minimap2 preset parameters')
    parser.add_argument('-p', '--threads', type=int, default=8,
                       help='线程数 | Number of threads')
    parser.add_argument('-m', '--min-match', type=int, default=1000,
                       help='最小匹配长度阈值 | Minimum match length threshold')
    parser.add_argument('-u', '--min-unmapped', type=int, default=1000,
                       help='最小未比对区间长度阈值 | Minimum unmapped region length threshold')
    parser.add_argument('--tp-type', default='P', choices=['S', 'P', 'SP'],
                       help='保留的tp类型 | tp type to keep: S(secondary), P(primary), SP(both) - 默认P | default P')
    parser.add_argument('-M', '--minimap2-path', default='minimap2',
                       help='minimap2可执行文件路径 | minimap2 executable path')
    parser.add_argument('-S', '--seqkit-path', default='seqkit',
                       help='seqkit可执行文件路径 | seqkit executable path')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = Minimap2Analyzer(
        target_genome=args.target,
        query_genome=args.query,
        output_dir=args.output_dir,
        preset=args.preset,
        threads=args.threads,
        min_match_length=args.min_match,
        min_unmapped_length=args.min_unmapped,
        tp_type=args.tp_type,
        minimap2_path=args.minimap2_path,
        seqkit_path=args.seqkit_path
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
