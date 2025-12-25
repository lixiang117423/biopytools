"""
GenomeThreader 分析主程序模块 | GenomeThreader Analysis Main Module
"""

import argparse
import sys
from .config import GTHConfig
from .utils import GTHLogger, CommandRunner, check_dependencies
from .data_processing import SequenceProcessor, FileManager
from .alignment import GTHAligner, ConsensusProcessor
from .output_processing import OutputParser, SummaryGenerator

class GenomeThreaderAnalyzer:
    """GenomeThreader 分析主类 | Main GenomeThreader Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = GTHConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = GTHLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.sequence_processor = SequenceProcessor(self.config, self.logger, self.cmd_runner)
        self.file_manager = FileManager(self.config, self.logger)
        self.aligner = GTHAligner(self.config, self.logger, self.cmd_runner)
        self.consensus_processor = ConsensusProcessor(self.config, self.logger, self.cmd_runner)
        self.output_parser = OutputParser(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的GenomeThreader分析流程 | Run complete GenomeThreader analysis pipeline"""
        try:
            self.logger.info("🧬 开始GenomeThreader基因预测分析流程 | Starting GenomeThreader gene prediction analysis pipeline")
            self.logger.info(f"📂 基因组文件 | Genomic file: {self.config.genomic_file}")
            self.logger.info(f"📁 输出目录 | Output directory: {self.config.output_dir}")
            
            if self.config.cdna_file:
                self.logger.info(f"🧵 cDNA文件 | cDNA file: {self.config.cdna_file}")
            if self.config.protein_file:
                self.logger.info(f"🔤 蛋白质文件 | Protein file: {self.config.protein_file}")
            if self.config.est_file:
                self.logger.info(f"📝 EST文件 | EST file: {self.config.est_file}")
            
            self.logger.info(f"⚡ 线程数 | Threads: {self.config.threads}")
            if self.config.species:
                self.logger.info(f"🐾 物种模型 | Species model: {self.config.species}")

            # 步骤1: 检查依赖软件
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("步骤1: 检查依赖软件 | Step 1: Checking dependencies")
            self.logger.info(f"{'=' * 60}")
            self.check_dependencies()

            # 步骤2: 设置输出目录结构
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("步骤2: 设置输出目录结构 | Step 2: Setting up output directory structure")
            self.logger.info(f"{'=' * 60}")
            self.file_manager.setup_output_structure()

            # 步骤3: 验证输入文件
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("步骤3: 验证输入文件 | Step 3: Validating input files")
            self.logger.info(f"{'=' * 60}")
            if not self.sequence_processor.validate_fasta_files():
                raise RuntimeError("输入文件验证失败 | Input file validation failed")

            # 步骤4: 收集序列统计信息
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("步骤4: 收集序列统计信息 | Step 4: Collecting sequence statistics")
            self.logger.info(f"{'=' * 60}")
            sequence_stats = self.sequence_processor.get_sequence_statistics()

            # 步骤5: 运行GenomeThreader比对分析
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("步骤5: 运行GenomeThreader比对分析 | Step 5: Running GenomeThreader alignment analysis")
            self.logger.info(f"{'=' * 60}")
            if not self.aligner.run_alignment():
                raise RuntimeError("GenomeThreader比对分析失败 | GenomeThreader alignment analysis failed")

            # 步骤6: 处理中间结果（如果启用）
            if self.config.intermediate:
                self.logger.info(f"\n{'=' * 60}")
                self.logger.info("步骤6: 处理中间结果 | Step 6: Processing intermediate results")
                self.logger.info(f"{'=' * 60}")
                if not self.consensus_processor.process_intermediate_results():
                    self.logger.warning("⚠️ 中间结果处理失败，但基本分析已完成 | Intermediate results processing failed, but basic analysis completed")

            # 步骤7: 解析和验证输出结果
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("步骤7: 解析和验证输出结果 | Step 7: Parsing and validating output results")
            self.logger.info(f"{'=' * 60}")
            output_stats = self.output_parser.parse_gff3_output()
            
            if not self.output_parser.validate_output_quality():
                self.logger.warning("⚠️ 输出质量检查发现问题，但分析已完成 | Output quality check found issues, but analysis completed")

            # 步骤8: 生成总结报告
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("步骤8: 生成总结报告 | Step 8: Generating summary report")
            self.logger.info(f"{'=' * 60}")
            self.summary_generator.generate_analysis_summary(sequence_stats, output_stats)
            self.summary_generator.generate_json_summary(sequence_stats, output_stats)

            # 步骤9: 清理临时文件
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info("步骤9: 清理临时文件 | Step 9: Cleaning up temporary files")
            self.logger.info(f"{'=' * 60}")
            self.file_manager.cleanup_temp_files()

            # 分析完成
            self.logger.info(f"\n{'=' * 80}")
            self.logger.info("🎉 GenomeThreader基因预测分析完成！| GenomeThreader gene prediction analysis completed!")
            self.logger.info(f"{'=' * 80}")
            self.logger.info("📁 输出文件:")
            self.logger.info(f"  📊 分析总结报告: analysis_summary.txt")
            self.logger.info(f"  📊 JSON格式总结: analysis_summary.json")
            
            if output_stats:
                self.logger.info(f"  🧬 预测基因数: {output_stats.get('genes', 0)}")
                self.logger.info(f"  🧵 预测转录本数: {output_stats.get('mrnas', 0)}")
            
            self.logger.info(f"  📁 详细结果目录: {self.config.output_dir}")
            self.logger.info(f"  📝 日志文件: {self.logger_manager.log_file}")
            self.logger.info(f"{'=' * 80}")
            self.logger.info(f"✅ 结果保存在: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"❌ 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 GenomeThreader基因预测分析脚本 (模块化版本) | GenomeThreader Gene Prediction Analysis Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
🧬 示例 | Examples:
  %(prog)s -g genome.fa -c cdna.fa -o gth_results
  %(prog)s -g genome.fa -p proteins.fa -s arabidopsis -o results
  %(prog)s -g genome.fa -c cdna.fa -p proteins.fa -s human --paralogs -o results
  %(prog)s -g genome.fa -c cdna.fa --fastdp --intron-cutout -o fast_results
  %(prog)s -g genome.fa -c est.fa --intermediate --min-score 0.8 -o high_quality
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-g', '--genomic', required=True, 
                       help='🧬 输入基因组序列文件 (FASTA格式) | Input genomic sequence file (FASTA format)')
    
    # 序列文件参数 | Sequence file arguments (至少需要一个 | at least one required)
    sequence_group = parser.add_argument_group('序列文件参数 | Sequence file arguments (至少需要一个 | at least one required)')
    sequence_group.add_argument('-c', '--cdna',
                               help='🧵 cDNA/转录本序列文件 | cDNA/transcript sequence file')
    sequence_group.add_argument('-p', '--protein',
                               help='🔤 蛋白质序列文件 | Protein sequence file')
    sequence_group.add_argument('-e', '--est',
                               help='📝 EST序列文件 | EST sequence file')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./gth_output', 
                       help='📁 输出目录 | Output directory')
    parser.add_argument('-s', '--species',
                       choices=['human', 'mouse', 'rat', 'chicken', 'drosophila', 'nematode', 
                               'fission_yeast', 'aspergillus', 'arabidopsis', 'maize', 'rice', 'medicago'],
                       help='🐾 物种模型 | Species model')
    
    # 分析参数 | Analysis parameters
    analysis_group = parser.add_argument_group('分析参数 | Analysis parameters')
    analysis_group.add_argument('-t', '--threads', type=int, default=88, 
                               help='⚡ 线程数 | Number of threads')
    analysis_group.add_argument('--min-score', type=float, default=0.65, 
                               help='🎯 最小比对得分阈值 | Minimum alignment score threshold')
    analysis_group.add_argument('--gc-min-coverage', type=int, default=50, 
                               help='📊 全局链最小覆盖度百分比 | Global chain minimum coverage percentage')
    
    # 高级分析选项 | Advanced analysis options
    advanced_group = parser.add_argument_group('高级分析选项 | Advanced analysis options')
    advanced_group.add_argument('--paralogs', action='store_true',
                               help='🔄 计算旁系同源基因 | Compute paralogous genes')
    advanced_group.add_argument('--intron-cutout', action='store_true',
                               help='✂️ 启用内含子切除技术 | Enable intron cutout technique')
    advanced_group.add_argument('--auto-intron-cutout', type=int, default=0,
                               help='🤖 自动内含子切除矩阵大小(MB) | Auto intron cutout matrix size (MB)')
    advanced_group.add_argument('--fastdp', action='store_true',
                               help='🚀 使用快速DP算法 | Use fast DP algorithm')
    
    # 链方向参数 | Strand direction parameters
    strand_group = parser.add_argument_group('链方向参数 | Strand direction parameters')
    strand_group.add_argument('--forward-only', action='store_true',
                             help='➡️ 只分析正链 | Analyze forward strand only')
    strand_group.add_argument('--reverse-only', action='store_true',
                             help='⬅️ 只分析负链 | Analyze reverse strand only')
    strand_group.add_argument('--cdna-forward', action='store_true',
                             help='🧵➡️ 只比对cDNA正链 | Align cDNA forward strand only')
    
    # 输出格式参数 | Output format parameters
    output_group = parser.add_argument_group('输出格式参数 | Output format parameters')
    output_group.add_argument('--xml-output', action='store_true',
                             help='📄 输出XML格式 | Output XML format')
    output_group.add_argument('--intermediate', action='store_true',
                             help='🔄 输出中间结果用于后续处理 | Output intermediate results for further processing')
    output_group.add_argument('--skip-alignment-out', action='store_true',
                             help='⏩ 跳过比对输出 | Skip alignment output')
    
    # 高级参数 | Advanced parameters
    advanced_params = parser.add_argument_group('高级参数 | Advanced parameters')
    advanced_params.add_argument('--score-matrix', default='BLOSUM62',
                                help='🧮 氨基酸替换评分矩阵 | Amino acid substitution scoring matrix')
    advanced_params.add_argument('--translation-table', type=int, default=1,
                                help='🔤 密码子翻译表 | Codon translation table')
    advanced_params.add_argument('--bssm-file',
                                help='📊 BSSM参数文件路径 | BSSM parameter file path')
    
    # 序列范围参数 | Sequence range parameters
    range_group = parser.add_argument_group('序列范围参数 | Sequence range parameters')
    range_group.add_argument('--from-pos', type=int, default=1,
                            help='🎯 起始位置 | Starting position')
    range_group.add_argument('--to-pos', type=int,
                            help='🏁 结束位置 | Ending position')
    range_group.add_argument('--width', type=int, default=1000000,
                            help='📏 序列处理宽度 | Sequence processing width')
    
    # 输出限制参数 | Output limitation parameters
    limit_group = parser.add_argument_group('输出限制参数 | Output limitation parameters')
    limit_group.add_argument('--first-alignments', type=int, default=0,
                            help='🔢 每个基因组序列最大比对数 (0=无限制) | Max alignments per genomic sequence (0=unlimited)')
    
    # 工具路径参数 | Tool path parameters
    tool_group = parser.add_argument_group('工具路径参数 | Tool path parameters')
    tool_group.add_argument('--gth-path', default='gth',
                           help='🔧 GenomeThreader程序路径 | GenomeThreader program path')
    tool_group.add_argument('--gthconsensus-path', default='gthconsensus',
                           help='🔧 gthconsensus程序路径 | gthconsensus program path')
    
    args = parser.parse_args()
    
    # 验证至少提供了一个序列文件
    if not any([args.cdna, args.protein, args.est]):
        parser.error("❌ 错误：必须至少提供一个序列文件 (--cdna, --protein, 或 --est) | Error: Must provide at least one sequence file (--cdna, --protein, or --est)")
    
    # 验证链方向参数
    if args.forward_only and args.reverse_only:
        parser.error("❌ 错误：不能同时指定 --forward-only 和 --reverse-only | Error: Cannot specify both --forward-only and --reverse-only")
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = GenomeThreaderAnalyzer(
        genomic_file=args.genomic,
        output_dir=args.output,
        cdna_file=args.cdna,
        protein_file=args.protein,
        est_file=args.est,
        species=args.species,
        threads=args.threads,
        min_alignment_score=args.min_score,
        gc_min_coverage=args.gc_min_coverage,
        paralogs=args.paralogs,
        intron_cutout=args.intron_cutout,
        auto_intron_cutout=args.auto_intron_cutout,
        fastdp=args.fastdp,
        forward_only=args.forward_only,
        reverse_only=args.reverse_only,
        cdna_forward=args.cdna_forward,
        xml_output=args.xml_output,
        intermediate=args.intermediate,
        skip_alignment_out=args.skip_alignment_out,
        score_matrix=args.score_matrix,
        translation_table=args.translation_table,
        bssm_file=args.bssm_file,
        from_pos=args.from_pos,
        to_pos=args.to_pos,
        width=args.width,
        first_alignments=args.first_alignments,
        gth_path=args.gth_path,
        gthconsensus_path=args.gthconsensus_path
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
