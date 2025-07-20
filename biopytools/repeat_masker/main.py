"""
重复序列分析主程序模块 | Repeat Sequence Analysis Main Module
"""

import argparse
import sys
from .config import RepeatConfig
from .utils import RepeatLogger, CommandRunner, SequenceValidator
from .analysis import RepeatModelerRunner, RepeatMaskerRunner, TRFRunner, EDTARunner
from .results import SummaryGenerator

class RepeatAnalyzer:
    """重复序列分析主类 | Main Repeat Sequence Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = RepeatConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = RepeatLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化验证器 | Initialize validator
        self.validator = SequenceValidator(self.logger)
        
        # 初始化各个分析器 | Initialize analyzers
        self.modeler_runner = RepeatModelerRunner(self.config, self.logger, self.cmd_runner)
        self.masker_runner = RepeatMaskerRunner(self.config, self.logger, self.cmd_runner)
        self.trf_runner = TRFRunner(self.config, self.logger, self.cmd_runner)
        self.edta_runner = EDTARunner(self.config, self.logger, self.cmd_runner)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的重复序列分析流程 | Run complete repeat sequence analysis pipeline"""
        try:
            self.logger.info("="*20 + " 开始基因组重复序列分析 | Starting Genome Repeat Sequence Analysis " + "="*20)
            
            # 验证输入文件 | Validate input files
            self.logger.info("验证输入文件 | Validating input files")
            if not self.validator.validate_fasta(self.config.genome_file):
                sys.exit(1)
            
            # 输出基因组统计信息 | Output genome statistics
            genome_stats = self.validator.get_sequence_stats(self.config.genome_file)
            self.logger.info(f"基因组统计 | Genome statistics:")
            self.logger.info(f"  序列数 | Number of sequences: {genome_stats['num_sequences']}")
            self.logger.info(f"  总长度 | Total length: {genome_stats['total_length']:,} bp")
            self.logger.info(f"  GC含量 | GC content: {genome_stats['gc_content']:.2f}%")
            self.logger.info(f"  N含量 | N content: {genome_stats['n_content']:.2f}%")
            
            # 步骤1: RepeatModeler从头预测 | Step 1: RepeatModeler de novo prediction
            if self.config.run_de_novo and self.config.rm_run_modeler:
                self.logger.info("\n步骤1: RepeatModeler从头预测重复序列 | Step 1: RepeatModeler de novo repeat prediction")
                
                if not self.modeler_runner.build_database():
                    self.logger.warning("RepeatModeler数据库构建失败，跳过从头预测 | RepeatModeler database building failed, skipping de novo prediction")
                elif not self.modeler_runner.run_repeat_modeler():
                    self.logger.warning("RepeatModeler运行失败 | RepeatModeler execution failed")
            
            # 步骤2: RepeatMasker数据库搜索和屏蔽 | Step 2: RepeatMasker database search and masking
            if self.config.run_database_search:
                self.logger.info("\n步骤2: RepeatMasker重复序列鉴定和屏蔽 | Step 2: RepeatMasker repeat identification and masking")
                
                if not self.masker_runner.run_repeat_masker():
                    self.logger.error("RepeatMasker运行失败 | RepeatMasker execution failed")
                    sys.exit(1)
            
            # 步骤3: TRF串联重复分析 | Step 3: TRF tandem repeat analysis
            if self.config.run_tandem_repeats and self.config.trf_run_analysis:
                self.logger.info("\n步骤3: TRF串联重复序列分析 | Step 3: TRF tandem repeat analysis")
                
                if not self.trf_runner.run_trf():
                    self.logger.warning("TRF分析失败 | TRF analysis failed")
            
            # 步骤4: EDTA转座元件分析 | Step 4: EDTA transposable element analysis
            if self.config.run_edta and self.config.edta_run_analysis:
                self.logger.info("\n步骤4: EDTA转座元件分析 | Step 4: EDTA transposable element analysis")
                
                if not self.edta_runner.run_edta():
                    self.logger.warning("EDTA分析失败 | EDTA analysis failed")
            
            # 步骤5: 生成分析报告 | Step 5: Generate analysis report
            self.logger.info("\n步骤5: 生成分析报告 | Step 5: Generating analysis report")
            self.summary_generator.generate_summary_report()
            self.summary_generator.export_statistics_table()
            
            self.logger.info("\n" + "="*20 + " 重复序列分析完成 | Repeat Sequence Analysis Completed " + "="*20)
            self.logger.info(f"结果保存在 | Results saved in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="基因组重复序列分析工具 (模块化版本) | Genome Repeat Sequence Analysis Tool (Modular Version)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 必需参数 | Required arguments
    parser.add_argument("-g", "--genome", required=True, 
                       help="输入基因组FASTA文件 | Input genome FASTA file")
    
    # 基本参数 | Basic parameters
    parser.add_argument("-o", "--output", default="./repeat_output", 
                       help="输出目录 | Output directory")
    parser.add_argument("-s", "--species", default="human", 
                       help="物种名称 | Species name")
    parser.add_argument("-l", "--lib", 
                       help="自定义重复序列库文件 | Custom repeat library file")
    
    # RepeatMasker参数 | RepeatMasker parameters
    parser.add_argument("-t", "--threads", type=int, default=8, 
                       help="RepeatMasker线程数 | RepeatMasker thread count")
    parser.add_argument("--hard-mask", action="store_true", 
                       help="使用硬屏蔽（N字符）而不是软屏蔽 | Use hard masking (N) instead of soft masking")
    parser.add_argument("--no-gff", action="store_true", 
                       help="不生成GFF注释文件 | Do not generate GFF annotation file")
    parser.add_argument("--include-low-complexity", action="store_true", 
                       help="包含低复杂度序列 | Include low complexity sequences")
    
    # RepeatModeler参数 | RepeatModeler parameters
    parser.add_argument("-m", "--modeler-threads", type=int, default=8, 
                       help="RepeatModeler线程数 | RepeatModeler thread count")
    parser.add_argument("--no-ltr-struct", action="store_true", 
                       help="不使用LTRStruct | Do not use LTRStruct")
    parser.add_argument("--no-modeler", action="store_true", 
                       help="跳过RepeatModeler从头预测 | Skip RepeatModeler de novo prediction")
    
    # TRF参数 | TRF parameters
    parser.add_argument("--trf-match", type=int, default=2, 
                       help="TRF匹配权重 | TRF match weight")
    parser.add_argument("--trf-mismatch", type=int, default=7, 
                       help="TRF错配惩罚 | TRF mismatch penalty")
    parser.add_argument("--trf-indel", type=int, default=7, 
                       help="TRF插入缺失惩罚 | TRF indel penalty")
    parser.add_argument("--trf-min-score", type=int, default=80, 
                       help="TRF最小比对得分 | TRF minimum alignment score")
    parser.add_argument("--trf-max-period", type=int, default=10, 
                       help="TRF最大周期长度 | TRF maximum period size")
    parser.add_argument("--trf-min-copies", type=int, default=50, 
                       help="TRF最小重复次数 | TRF minimum copy number")
    parser.add_argument("--trf-max-length", type=int, default=500, 
                       help="TRF最大期望长度 | TRF maximum expected length")
    parser.add_argument("--no-trf", action="store_true", 
                       help="跳过TRF串联重复分析 | Skip TRF tandem repeat analysis")
    
    # EDTA参数 | EDTA parameters
    parser.add_argument("-e", "--edta-threads", type=int, default=8, 
                       help="EDTA线程数 | EDTA thread count")
    parser.add_argument("--edta-species", choices=["rice", "maize", "others"], default="others",
                       help="EDTA物种类型 | EDTA species type")
    parser.add_argument("--edta-step", choices=["all", "filter", "final", "anno"], default="all",
                       help="EDTA执行步骤 | EDTA execution step")
    parser.add_argument("--edta-sensitive", action="store_true", 
                       help="EDTA敏感模式 | EDTA sensitive mode")
    parser.add_argument("--no-edta-anno", action="store_true", 
                       help="不进行EDTA注释 | Do not perform EDTA annotation")
    parser.add_argument("--no-edta-eval", action="store_true", 
                       help="不进行EDTA评估 | Do not perform EDTA evaluation")
    parser.add_argument("--edta-overwrite", action="store_true", 
                       help="覆盖现有EDTA结果 | Overwrite existing EDTA results")
    parser.add_argument("--no-edta", action="store_true", 
                       help="跳过EDTA转座元件分析 | Skip EDTA transposable element analysis")
    
    # 分析选项 | Analysis options
    parser.add_argument("--no-de-novo", action="store_true", 
                       help="跳过从头预测分析 | Skip de novo prediction analysis")
    parser.add_argument("--no-database", action="store_true", 
                       help="跳过数据库搜索 | Skip database search")
    parser.add_argument("--min-length", type=int, default=50, 
                       help="最小重复序列长度 | Minimum repeat length")
    parser.add_argument("--max-divergence", type=float, default=0.25, 
                       help="最大分化度 (0-1) | Maximum divergence (0-1)")
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = RepeatAnalyzer(
        genome_file=args.genome,
        output_dir=args.output,
        species=args.species,
        custom_lib=args.lib,
        
        # RepeatMasker参数 | RepeatMasker parameters
        rm_threads=args.threads,
        rm_soft_mask=not args.hard_mask,
        rm_generate_gff=not args.no_gff,
        rm_exclude_low_complexity=not args.include_low_complexity,
        
        # RepeatModeler参数 | RepeatModeler parameters
        rm_model_threads=args.modeler_threads,
        rm_use_ltr_struct=not args.no_ltr_struct,
        rm_run_modeler=not args.no_modeler,
        
        # TRF参数 | TRF parameters
        trf_match_weight=args.trf_match,
        trf_mismatch_penalty=args.trf_mismatch,
        trf_indel_penalty=args.trf_indel,
        trf_min_score=args.trf_min_score,
        trf_max_period=args.trf_max_period,
        trf_min_copies=args.trf_min_copies,
        trf_max_length=args.trf_max_length,
        trf_run_analysis=not args.no_trf,
        
        # EDTA参数 | EDTA parameters
        edta_threads=args.edta_threads,
        edta_species=args.edta_species,
        edta_step=args.edta_step,
        edta_sensitive=1 if args.edta_sensitive else 0,
        edta_anno=0 if args.no_edta_anno else 1,
        edta_evaluate=0 if args.no_edta_eval else 1,
        edta_overwrite=1 if args.edta_overwrite else 0,
        edta_run_analysis=not args.no_edta,
        
        # 分析选项 | Analysis options
        run_de_novo=not args.no_de_novo,
        run_database_search=not args.no_database,
        run_tandem_repeats=not args.no_trf,
        run_edta=not args.no_edta,
        min_repeat_length=args.min_length,
        max_divergence=args.max_divergence
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
