"""
Augustus多转录组预测主程序模块 | Augustus Multiple RNA-seq Prediction Main Module
"""

import argparse
import sys
from .config import AugustusConfig
from .utils import AugustusLogger, CommandRunner, check_dependencies
from .data_processing import HISAT2Manager, BAMProcessor, HintsGenerator
from .annotation import AugustusPredictor
from .results import ResultsComparator, ReportGenerator

class AugustusMultiRNASeq:
    """Augustus多转录组预测主类 | Main Augustus Multiple RNA-seq Prediction Class"""
    
    def __init__(self, skip_deps_check=False, **kwargs):
        # 检查依赖工具 | Check dependencies
        if skip_deps_check:
            print("跳过依赖检查 | Skipping dependency check")
            self.has_augustus_tools = True  # 假设工具可用
        else:
            self.has_augustus_tools = check_dependencies()
        
        # 初始化配置 | Initialize configuration
        self.config = AugustusConfig(**kwargs)
        
        # 如果没有augustus工具，强制禁用BAM过滤 | Force disable BAM filtering if no augustus tools
        if not self.has_augustus_tools:
            self.config.filter_bam = False
            print("注意：由于缺少filterBam工具，BAM过滤已被禁用 | Note: BAM filtering disabled due to missing filterBam tool")
        
        self.config.validate()
        
        # 解析样本配置 | Parse sample configuration
        self.config.parse_sample_config()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = AugustusLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 记录工具可用性 | Log tool availability
        if not self.has_augustus_tools and not skip_deps_check:
            self.logger.warning("Augustus辅助工具不可用，将使用基础模式运行 | Augustus auxiliary tools unavailable, running in basic mode")
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.hisat2_manager = HISAT2Manager(self.config, self.logger, self.cmd_runner)
        self.bam_processor = BAMProcessor(self.config, self.logger, self.cmd_runner, self.has_augustus_tools)
        self.hints_generator = HintsGenerator(self.config, self.logger, self.cmd_runner, self.has_augustus_tools)
        self.augustus_predictor = AugustusPredictor(self.config, self.logger, self.cmd_runner)
        self.results_comparator = ResultsComparator(self.config, self.logger, self.cmd_runner)
        self.report_generator = ReportGenerator(self.config, self.logger)
    
    def run_prediction(self):
        """运行完整的Augustus多转录组预测流程 | Run complete Augustus multiple RNA-seq prediction pipeline"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始多转录组Augustus基因预测流程 | Starting multiple RNA-seq Augustus gene prediction pipeline")
            self.logger.info("=" * 60)
            
            # 基本信息
            self.logger.info(f"基因组文件 | Genome file: {self.config.genome_file}")
            self.logger.info(f"Augustus模型 | Augustus model: {self.config.species_model}")
            self.logger.info(f"转录组样本数 | RNA-seq samples: {len(self.config.samples)}")
            self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
            
            # 步骤1: 构建HISAT2索引
            self.logger.info("\n" + "=" * 60)
            self.logger.info("步骤1: 检查和构建HISAT2索引 | Step 1: Check and build HISAT2 index")
            self.logger.info("=" * 60)
            
            if not self.hisat2_manager.build_hisat2_index():
                self.logger.error("HISAT2索引构建失败 | HISAT2 index building failed")
                sys.exit(1)
            
            # 步骤2: 处理转录组样本
            self.logger.info("\n" + "=" * 60)
            self.logger.info("步骤2: 处理转录组样本 | Step 2: Process RNA-seq samples")
            self.logger.info("=" * 60)
            
            bam_files = []
            for sample in self.config.samples:
                # 比对样本
                success, bam_file = self.hisat2_manager.align_sample(sample)
                if not success:
                    self.logger.error(f"样本比对失败 | Sample alignment failed: {sample['name']}")
                    continue
                
                # 过滤BAM文件
                filtered_bam = self.bam_processor.filter_bam(bam_file)
                bam_files.append(filtered_bam)
                
                # 获取比对统计
                self.bam_processor.get_alignment_stats(filtered_bam)
            
            self.logger.info(f"所有样本处理完成，共 {len(bam_files)} 个BAM文件 | All samples processed, {len(bam_files)} BAM files total")
            
            # 步骤3: 生成hints文件
            self.logger.info("\n" + "=" * 60)
            self.logger.info("步骤3: 生成hints文件 | Step 3: Generate hints files")
            self.logger.info("=" * 60)
            
            hints_files = []
            for bam_file in bam_files:
                hints_file = self.hints_generator.generate_hints_for_sample(bam_file)
                if hints_file:
                    hints_files.append(hints_file)
            
            # 合并和过滤hints
            filtered_hints = self.hints_generator.merge_and_filter_hints(hints_files)
            
            # 步骤4: 运行Augustus预测
            self.logger.info("\n" + "=" * 60)
            self.logger.info("步骤4: 运行Augustus预测 | Step 4: Run Augustus prediction")
            self.logger.info("=" * 60)
            
            # 方案A: 仅使用Augustus模型
            model_only_result = self.augustus_predictor.run_augustus_model_only()
            
            # 方案B: 结合转录组数据
            rnaseq_result = self.augustus_predictor.run_augustus_with_rnaseq(filtered_hints)
            
            # 步骤5: 结果比较和报告生成
            self.logger.info("\n" + "=" * 60)
            self.logger.info("步骤5: 结果比较和报告生成 | Step 5: Results comparison and report generation")
            self.logger.info("=" * 60)
            
            # 比较预测结果
            comparison_results = self.results_comparator.compare_predictions(
                model_only_result, rnaseq_result, filtered_hints
            )
            
            # 生成最终报告
            self.report_generator.generate_final_report(comparison_results, bam_files)
            
            # 完成
            self.logger.info("\n" + "=" * 60)
            self.logger.info("多转录组Augustus基因预测流程完成 | Multiple RNA-seq Augustus gene prediction pipeline completed")
            self.logger.info("=" * 60)
            
            self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
            self.logger.info("主要结果文件 | Main result files:")
            self.logger.info("  - augustus_model_only.gff3: 仅模型预测结果 | Model-only prediction")
            self.logger.info("  - augustus_with_multi_rnaseq.gff3: 结合转录组预测结果 | Prediction with RNA-seq")
            self.logger.info("  - filtered_hints.gff: 过滤后的hints文件 | Filtered hints file")
            self.logger.info("  - final_report.txt: 最终分析报告 | Final analysis report")
            
            self.logger.info("\n建议后续分析 | Suggested further analysis:")
            self.logger.info("- 使用BUSCO评估基因集完整性 | Use BUSCO to assess gene set completeness")
            self.logger.info("- 使用gffcompare与参考注释比较 | Use gffcompare to compare with reference annotation")
            
        except Exception as e:
            self.logger.error(f"预测流程在执行过程中意外终止 | Prediction pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="多转录组Augustus基因预测脚本 | Multiple RNA-seq Augustus Gene Prediction Script",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 必需参数 | Required arguments
    parser.add_argument("-g", "--genome", required=True, 
                       help="基因组fasta文件路径 | Genome fasta file path")
    parser.add_argument("-s", "--species", required=True, 
                       help="Augustus训练的物种模型名称 | Augustus trained species model name")
    
    # 输入方式 (二选一) | Input method (choose one)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-i", "--input-dir", 
                            help="输入FASTQ文件目录 | Input FASTQ files directory")
    input_group.add_argument("-c", "--config", 
                            help="样本配置文件路径 | Sample configuration file path")
    
    # 文件模式参数 (仅在使用-i时有效) | File pattern (only valid when using -i)
    parser.add_argument("-p", "--pattern", default="*.R1.fastq.gz",
                       help="R1文件匹配模式 | R1 file matching pattern (only used with -i/--input-dir)")
    
    # 可选参数 | Optional arguments
    parser.add_argument("-o", "--output", default="augustus_multi_rnaseq", 
                       help="输出目录 | Output directory")
    parser.add_argument("-t", "--threads", type=int, default=8, 
                       help="线程数 | Number of threads")
    parser.add_argument("-x", "--hisat2-index", default=None, 
                       help="HISAT2索引路径前缀 | HISAT2 index path prefix")
    parser.add_argument("-m", "--min-intron-support", type=int, default=2, 
                       help="内含子hints的最小支持度 | Minimum support for intron hints")
    parser.add_argument("-f", "--no-filter-bam", action="store_true", 
                       help="跳过BAM文件过滤步骤 | Skip BAM filtering step")
    parser.add_argument("-a", "--no-alternatives", action="store_true", 
                       help="不使用转录组证据生成可变剪切 | Do not use RNA-seq evidence for alternative splicing")
    parser.add_argument("-z", "--splicesites", default="atac", choices=["atac", "gtag", "gcag"], 
                       help="允许的剪切位点类型 | Allowed splice site types")
    parser.add_argument("--skip-deps-check", action="store_true",
                       help="跳过依赖检查 | Skip dependency check")
    
    args = parser.parse_args()
    
    # 创建预测器并运行 | Create predictor and run
    predictor = AugustusMultiRNASeq(
        genome_file=args.genome,
        species_model=args.species,
        input_dir=args.input_dir,
        config_file=args.config,
        pattern=args.pattern,
        output_dir=args.output,
        threads=args.threads,
        hisat2_index=args.hisat2_index,
        min_intron_support=args.min_intron_support,
        filter_bam=not args.no_filter_bam,
        alternatives_from_evidence=not args.no_alternatives,
        allow_hinted_splicesites=args.splicesites,
        skip_deps_check=args.skip_deps_check
    )
    
    predictor.run_prediction()

if __name__ == "__main__":
    main()
