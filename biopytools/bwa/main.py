# ===== FILE: bwa_align/main.py =====
"""
BWA比对主程序模块|BWA Alignment Main Module
"""

import sys
from .config import AlignConfig
from .utils import AlignLogger, CommandRunner, check_dependencies, find_fastq_pairs
from .genome import GenomeIndexer
from .alignment import BWAAlignmentProcessor
from .coverage import CoverageAnalyzer
from .stats import AlignmentStatsGenerator

class BWAAligner:
    """BWA比对器主类|Main BWA Aligner Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = AlignConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志|Initialize logging
        self.logger_manager = AlignLogger(self.config.log_dir)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器|Initialize processors
        self.genome_indexer = GenomeIndexer(self.config, self.logger, self.cmd_runner)
        self.aligner = BWAAlignmentProcessor(self.config, self.logger, self.cmd_runner)
        self.coverage_analyzer = CoverageAnalyzer(self.config, self.logger, self.cmd_runner)
        self.stats_generator = AlignmentStatsGenerator(self.config, self.logger, self.cmd_runner)
    
    def check_dependencies(self):
        """检查依赖软件|Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的比对分析流程|Run complete alignment analysis pipeline"""
        try:
            self.logger.info("=" * 80)
            self.logger.info("开始BWA全基因组比对分析流程|Starting BWA alignment pipeline")
            self.logger.info("=" * 80)

            # Step 1: 检查依赖|Check dependencies
            self.logger.info("步骤1/6: 检查依赖软件|Step 1/6: Checking dependencies")
            self.check_dependencies()

            # Step 2: 检查并构建基因组索引|Check and build genome index
            self.logger.info("步骤2/6: 检查基因组索引|Step 2/6: Checking genome index")
            if not self.genome_indexer.check_and_build_index():
                raise RuntimeError("基因组索引构建失败|Genome index building failed")

            # Step 3: 查找FASTQ文件对|Find FASTQ pairs
            self.logger.info("步骤3/6: 查找FASTQ文件|Step 3/6: Finding FASTQ files")
            fastq_pairs = find_fastq_pairs(self.config.input_dir,
                                          self.config.pattern,
                                          self.logger)

            if not fastq_pairs:
                raise RuntimeError("未找到匹配的FASTQ文件对|No matching FASTQ pairs found")

            # Step 4: 批量比对|Batch alignment
            self.logger.info(f"步骤4/6: 批量比对 ({len(fastq_pairs)} 个样品)|"
                           f"Step 4/6: Batch alignment ({len(fastq_pairs)} samples)")

            processed_samples = []
            for i, (sample_name, read1, read2) in enumerate(fastq_pairs, 1):
                self.logger.info(f"{'='*80}")
                self.logger.info(f"处理样品 {i}/{len(fastq_pairs)}: {sample_name}|"
                               f"Processing sample {i}/{len(fastq_pairs)}: {sample_name}")
                self.logger.info(f"{'='*80}")

                # 比对|Alignment
                bam_file = self.aligner.align_sample(sample_name, read1, read2)
                if not bam_file:
                    self.logger.error(f"样品比对失败|Sample alignment failed: {sample_name}")
                    continue

                processed_samples.append(sample_name)

                # Step 5: 统计|Statistics
                self.stats_generator.generate_stats(sample_name, bam_file)

                # Step 6: 覆盖度分析|Coverage analysis
                self.coverage_analyzer.analyze_sample_coverage(sample_name, bam_file)

            # 生成汇总报告|Generate summary report
            self.logger.info("步骤6/6: 生成汇总报告|Step 6/6: Generating summary report")
            self.stats_generator.generate_summary_report(processed_samples)

            self.logger.info("=" * 80)
            self.logger.info("BWA全基因组比对分析完成|BWA alignment analysis completed")
            self.logger.info("=" * 80)
            self.logger.info(f"成功处理样品数|Successfully processed: {len(processed_samples)}/{len(fastq_pairs)}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

            return True

        except Exception as e:
            self.logger.error(f"分析流程出错|Analysis pipeline error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False

def main():
    """命令行入口函数|Command-line entry function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='BWA全基因组比对分析工具|BWA Whole Genome Alignment Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -g genome.fa -i fastq_dir -p _1.clean.fq.gz
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument('-g', '--genome', required=True,
                         help='参考基因组文件|Reference genome file')
    required.add_argument('-i', '--input', required=True,
                         help='输入FASTQ文件目录|Input FASTQ directory')
    required.add_argument('-p', '--pattern', required=True,
                         help='FASTQ文件匹配模式|FASTQ file pattern')

    # 输出参数|Output parameters
    output = parser.add_argument_group('输出参数|Output parameters')
    output.add_argument('-o', '--output-dir', default='./bwa_output',
                       help='输出目录|Output directory')

    # 性能参数|Performance parameters
    perf = parser.add_argument_group('性能参数|Performance parameters')
    perf.add_argument('-t', '--threads', type=int, default=88,
                     help='线程数|Number of threads')

    # BWA算法参数|BWA algorithm parameters
    bwa_algo = parser.add_argument_group('BWA算法参数|BWA algorithm parameters')
    bwa_algo.add_argument('--bwa-k', type=int, default=19,
                         help='最小种子长度|Minimum seed length')
    bwa_algo.add_argument('--bwa-w', type=int, default=100,
                         help='带宽|Band width')
    bwa_algo.add_argument('--bwa-d', type=int, default=100,
                         help='X-dropoff|Off-diagonal X-dropoff')
    bwa_algo.add_argument('--bwa-r', type=float, default=1.5,
                         help='内部种子因子|Internal seed factor')
    bwa_algo.add_argument('--bwa-c', type=int, default=500,
                         help='种子出现次数阈值|Seed occurrence threshold')
    bwa_algo.add_argument('--bwa-D', type=float, default=0.50, dest='bwa_drop_ratio',
                         help='短链丢弃比例|Short chain drop fraction')
    bwa_algo.add_argument('--bwa-W', type=int, default=0, dest='bwa_min_chain',
                         help='最小链长|Minimum chain length')
    bwa_algo.add_argument('--bwa-m', type=int, default=50,
                         help='配对拯救轮数|Mate rescue rounds')
    bwa_algo.add_argument('--bwa-S', action='store_true', dest='bwa_s',
                         help='跳过配对拯救|Skip mate rescue')
    bwa_algo.add_argument('--bwa-P', action='store_true', dest='bwa_p',
                         help='跳过配对|Skip pairing')

    # BWA打分参数|BWA scoring parameters
    bwa_score = parser.add_argument_group('BWA打分参数|BWA scoring parameters')
    bwa_score.add_argument('--bwa-A', type=int, default=1, dest='bwa_score_match',
                          help='匹配得分|Match score')
    bwa_score.add_argument('--bwa-B', type=int, default=4, dest='bwa_score_mismatch',
                          help='错配罚分|Mismatch penalty')
    bwa_score.add_argument('--bwa-O', default="6,6", dest='bwa_gap_open',
                          help='Gap开放罚分|Gap open penalty')
    bwa_score.add_argument('--bwa-E', default="1,1", dest='bwa_gap_ext',
                          help='Gap延伸罚分|Gap extension penalty')
    bwa_score.add_argument('--bwa-L', default="5,5", dest='bwa_clip',
                          help='末端剪切罚分|Clipping penalty')
    bwa_score.add_argument('--bwa-U', type=int, default=17, dest='bwa_unpaired',
                          help='未配对罚分|Unpaired penalty')

    # BWA输出参数|BWA output parameters
    bwa_io = parser.add_argument_group('BWA输出参数|BWA output parameters')
    bwa_io.add_argument('--bwa-M', action='store_true', dest='bwa_mark_secondary',
                       help='标记次要比对|Mark shorter split hits as secondary')
    bwa_io.add_argument('--bwa-T', type=int, default=30, dest='bwa_min_score',
                       help='最小输出得分|Minimum score to output')
    bwa_io.add_argument('--bwa-a', action='store_true', dest='bwa_all_align',
                       help='输出所有比对|Output all alignments')
    bwa_io.add_argument('--bwa-C', action='store_true', dest='bwa_append_comment',
                       help='附加FASTQ注释|Append FASTA/FASTQ comment')
    bwa_io.add_argument('--bwa-V', action='store_true', dest='bwa_ref_header',
                       help='输出参考序列头|Output reference FASTA header')
    bwa_io.add_argument('--bwa-Y', action='store_true', dest='bwa_soft_clip',
                       help='软剪切补充比对|Soft clipping for supplementary alignments')

    # 后处理参数|Post-processing parameters
    post = parser.add_argument_group('后处理参数|Post-processing parameters')
    post.add_argument('--markdup', action='store_true',
                     help='标记重复序列|Mark duplicate reads')
    post.add_argument('--remove-dup', action='store_true',
                     help='移除重复序列|Remove duplicate reads')

    # 覆盖度参数|Coverage parameters
    cov = parser.add_argument_group('覆盖度参数|Coverage parameters')
    cov.add_argument('--min-base-quality', type=int, default=0,
                    help='最小碱基质量|Minimum base quality')
    cov.add_argument('--min-mapping-quality', type=int, default=0,
                    help='最小比对质量|Minimum mapping quality')
    cov.add_argument('--max-depth', type=int, default=0,
                    help='最大深度限制|Max depth limit')

    # 滑窗参数|Window parameters
    window = parser.add_argument_group('滑窗参数|Window parameters')
    window.add_argument('--window-size', type=int, default=1000000,
                       help='窗口大小|Window size in bp')
    window.add_argument('--step-size', type=int, default=100000,
                       help='步长|Step size in bp')

    # 其他参数|Other parameters
    other = parser.add_argument_group('其他参数|Other parameters')
    other.add_argument('--resume', action='store_true',
                      help='断点续传|Resume skip completed samples')
    other.add_argument('--keep-sam', action='store_true',
                      help='保留SAM文件|Keep SAM files')

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    analyzer = BWAAligner(
        genome=args.genome,
        input_dir=args.input,
        pattern=args.pattern,
        output_dir=args.output_dir,
        threads=args.threads,
        # BWA算法参数
        bwa_k=args.bwa_k,
        bwa_w=args.bwa_w,
        bwa_d=args.bwa_d,
        bwa_r=args.bwa_r,
        bwa_c=args.bwa_c,
        bwa_drop_ratio=args.bwa_drop_ratio,
        bwa_min_chain=args.bwa_min_chain,
        bwa_m=args.bwa_m,
        bwa_s=args.bwa_s,
        bwa_p=args.bwa_p,
        # BWA打分参数
        bwa_score_match=args.bwa_score_match,
        bwa_score_mismatch=args.bwa_score_mismatch,
        bwa_gap_open=args.bwa_gap_open,
        bwa_gap_ext=args.bwa_gap_ext,
        bwa_clip=args.bwa_clip,
        bwa_unpaired=args.bwa_unpaired,
        # BWA输出参数
        bwa_mark_secondary=args.bwa_mark_secondary,
        bwa_min_score=args.bwa_min_score,
        bwa_all_align=args.bwa_all_align,
        bwa_append_comment=args.bwa_append_comment,
        bwa_ref_header=args.bwa_ref_header,
        bwa_soft_clip=args.bwa_soft_clip,
        # 后处理参数
        markdup=args.markdup,
        remove_dup=args.remove_dup,
        # 覆盖度参数
        min_base_quality=args.min_base_quality,
        min_mapping_quality=args.min_mapping_quality,
        max_depth=args.max_depth,
        # 滑窗参数
        window_size=args.window_size,
        step_size=args.step_size,
        # 其他参数
        resume=args.resume,
        keep_sam=args.keep_sam
    )
    
    success = analyzer.run_analysis()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()