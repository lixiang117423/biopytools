"""
ANNOVAR注释主程序模块|ANNOVAR Annotation Main Module
"""

import argparse
import sys
import os
from .config import ANNOVARConfig
from .utils import ANNOVARLogger, CommandRunner
from .data_processing import GFF3Processor, SequenceExtractor, VCFProcessor
from .annotation import VariantAnnotator
from .results import SummaryGenerator
from .results_processor import ANNOVARResultsProcessor

class ANNOVARAnnotator:
    """ANNOVAR注释主类|Main ANNOVAR Annotator Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = ANNOVARConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = ANNOVARLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        # 初始化各个处理器|Initialize processors
        self.gff3_processor = GFF3Processor(self.config, self.logger, self.cmd_runner)
        self.sequence_extractor = SequenceExtractor(self.config, self.logger, self.cmd_runner)
        self.vcf_processor = VCFProcessor(self.config, self.logger, self.cmd_runner)
        self.variant_annotator = VariantAnnotator(self.config, self.logger, self.cmd_runner)
        self.summary_generator = SummaryGenerator(self.config, self.logger)

        # 初始化结果处理器|Initialize results processor
        self.results_processor = ANNOVARResultsProcessor(self.logger, self.config.output_dir)

    def step1_gff3_to_genepred(self):
        """步骤1: GFF3转GenPred|Step 1: GFF3 to GenPred"""
        return self.gff3_processor.gff3_to_genepred()

    def step2_extract_transcript_sequences(self):
        """步骤2: 提取转录本序列|Step 2: Extract transcript sequences"""
        return self.sequence_extractor.extract_transcript_sequences()

    def step3_filter_and_convert_vcf(self):
        """步骤3: 过滤并转换VCF|Step 3: Filter and convert VCF"""
        return self.vcf_processor.filter_and_convert_vcf()

    def step4_annotate_variants(self):
        """步骤4: 注释变异|Step 4: Annotate variants"""
        return self.variant_annotator.annotate_variants()

    def run_single_step(self, step_num: int):
        """运行单个步骤|Run single step"""
        step_functions = {
            1: (self.step1_gff3_to_genepred, "GFF3转GenPred|GFF3 to GenPred"),
            2: (self.step2_extract_transcript_sequences, "提取转录本序列|Extract transcript sequences"),
            3: (self.step3_filter_and_convert_vcf, "过滤并转换VCF|Filter and convert VCF"),
            4: (self.step4_annotate_variants, "注释变异|Annotate variants")
        }

        if step_num not in step_functions:
            self.logger.error(f"无效的步骤编号|Invalid step number: {step_num}")
            return False

        step_func, step_name = step_functions[step_num]
        self.logger.info(f"执行步骤{step_num}|Executing step {step_num}: {step_name}")

        success = step_func()
        if success:
            self.logger.info(f"步骤{step_num}完成|Step {step_num} completed: {step_name}")
        else:
            self.logger.error(f"步骤{step_num}失败|Step {step_num} failed: {step_name}")

        return success

    def run_full_pipeline(self):
        """运行完整的ANNOVAR注释流程|Run complete ANNOVAR annotation pipeline"""
        self.logger.info("开始ANNOVAR注释流程|Starting ANNOVAR annotation pipeline")

        steps = [
            (self.step1_gff3_to_genepred, "GFF3转GenPred|GFF3 to GenPred"),
            (self.step2_extract_transcript_sequences, "提取转录本序列|Extract transcript sequences"),
            (self.step3_filter_and_convert_vcf, "处理并转换VCF|Process and convert VCF"),
            (self.step4_annotate_variants, "变异注释|Variant annotation")
        ]

        for i, (step_func, step_name) in enumerate(steps, 1):
            self.logger.info(f"执行步骤{i}|Executing step {i}: {step_name}")

            if not step_func():
                self.logger.error(f"步骤{i}失败|Step {i} failed: {step_name}")
                return False

            self.logger.info(f"步骤{i}完成|Step {i} completed: {step_name}")

        self.logger.info("ANNOVAR注释流程全部完成|ANNOVAR annotation pipeline completed!")
        self.summary_generator.generate_summary_report()

        # 自动处理注释结果|Automatically process annotation results
        if hasattr(self.config, 'vcf_basename'):
            self.logger.info("开始处理注释结果|Starting to process annotation results")
            self.process_annotation_results()

        return True
    
    def run_analysis(self):
        """运行分析|Run analysis"""
        try:
            if self.config.step:
                # 运行单个步骤|Run single step
                success = self.run_single_step(self.config.step)
                if success:
                    self.logger.info(f"步骤 {self.config.step} 执行成功|Step {self.config.step} executed successfully")
                else:
                    self.logger.error(f"步骤 {self.config.step} 执行失败|Step {self.config.step} execution failed")
                    sys.exit(1)
            else:
                # 运行完整流程|Run full pipeline
                success = self.run_full_pipeline()
                if not success:
                    sys.exit(1)

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            sys.exit(1)

    def process_annotation_results(self, apply_filters: bool = False):
        """处理注释结果|Process annotation results"""
        if not hasattr(self.config, 'vcf_basename'):
            self.logger.warning("缺少VCF基础名称，无法处理结果|Missing VCF basename, cannot process results")
            return {}

        try:
            processed_files = self.results_processor.process_available_results(
                self.config.vcf_basename, apply_filters
            )

            # 更新配置中的输出文件列表|Update output files list in configuration
            if not hasattr(self.config, 'processed_output_files'):
                self.config.processed_output_files = []

            self.config.processed_output_files.extend(processed_files.values())

            self.logger.info(f"注释结果处理完成|Annotation results processing completed")
            for file_type, file_path in processed_files.items():
                self.logger.info(f"  {file_type.upper()}: {file_path}")

            return processed_files

        except Exception as e:
            self.logger.error(f"处理注释结果时出错|Error processing annotation results: {str(e)}")
            return {}

    def process_exonic_results_only(self, exonic_file: str = None, output_prefix: str = None):
        """仅处理外显子注释结果|Process exonic annotation results only"""
        if exonic_file is None:
            if not hasattr(self.config, 'vcf_basename'):
                self.logger.error("缺少VCF基础名称，无法自动查找外显子文件|Missing VCF basename, cannot automatically find exonic file")
                return None
            exonic_file = os.path.join(self.config.output_dir, f"{self.config.vcf_basename}.exonic_variant_function")

        if output_prefix is None and hasattr(self.config, 'vcf_basename'):
            output_prefix = self.config.vcf_basename

        return self.results_processor.process_exonic_results(exonic_file, output_prefix)

    def process_all_results_only(self, variant_function_file: str = None, output_prefix: str = None,
                                apply_filters: bool = False, filters: dict = None):
        """仅处理所有变异注释结果|Process all variant annotation results only"""
        if variant_function_file is None:
            if not hasattr(self.config, 'vcf_basename'):
                self.logger.error("缺少VCF基础名称，无法自动查找变异功能文件|Missing VCF basename, cannot automatically find variant function file")
                return None
            variant_function_file = os.path.join(self.config.output_dir, f"{self.config.vcf_basename}.variant_function")

        if output_prefix is None and hasattr(self.config, 'vcf_basename'):
            output_prefix = self.config.vcf_basename

        return self.results_processor.process_all_results(
            variant_function_file, output_prefix, apply_filters, filters
        )

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='ANNOVAR VCF注释自动化脚本(模块化版本)|ANNOVAR VCF Annotation Automation Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='VCF变异文件路径|VCF variant file path')
    parser.add_argument('-g', '--gff3', required=True,
                       help='GFF3注释文件路径|GFF3 annotation file path')
    parser.add_argument('-G', '--genome', required=True,
                       help='基因组序列文件路径|Genome sequence file path')
    parser.add_argument('-b', '--build-ver', required=True,
                       help='基因组构建版本标识符(如: OV, KY131) - 不应包含路径分隔符|'
                            'Genome build version identifier (e.g., OV, KY131) - should not contain path separators')

    # 可选参数|Optional arguments
    parser.add_argument('-a', '--annovar-path',
                       default='/share/org/YZWL/yzwl_lixg/software/annovar/annovar',
                       help='ANNOVAR软件安装路径|ANNOVAR software installation path')
    parser.add_argument('-d', '--database-path',
                       default='./database',
                       help='ANNOVAR数据库路径|ANNOVAR database path')
    parser.add_argument('-o', '--output-dir',
                       default='./annovar_output',
                       help='输出目录|Output directory')
    parser.add_argument('-q', '--qual-threshold',
                       type=int, default=20,
                       help='VCF质量过滤阈值(仅在启用VCF过滤时生效)|'
                            'VCF quality filtering threshold (only effective when VCF filtering is enabled)')

    # 步骤控制|Step control
    parser.add_argument('-s', '--step', type=int, choices=[1, 2, 3, 4],
                       help='只运行指定步骤|Run only specified step '
                            '(1: gff3转换|gff3 conversion, 2: 提取序列|extract sequences, '
                            '3: VCF处理|VCF processing, 4: 注释|annotation)')

    # 处理选项|Processing options
    parser.add_argument('--skip-gff-cleaning', action='store_true',
                       help='跳过GFF3文件的格式清理(attributes清理和坐标修复)|'
                            'Skip GFF3 file format cleaning (attributes cleaning and coordinate fixing)')
    parser.add_argument('--skip-gff-fix', action='store_true',
                       help='跳过GFF3文件的自动修复(CDS phase等问题)|'
                            'Skip automatic GFF3 file fixes (CDS phase and other issues)')

    parser.add_argument('--skip-vcf-filter', action='store_true', default=True,
                       help='跳过VCF过滤步骤，直接使用输入的VCF文件(默认启用)|'
                            'Skip VCF filtering step, use input VCF file directly (enabled by default)')
    parser.add_argument('--enable-vcf-filter', action='store_true',
                       help='启用VCF过滤步骤(使用bcftools)|'
                            'Enable VCF filtering step (using bcftools)')

    args = parser.parse_args()

    # 处理VCF过滤选项|Handle VCF filtering options
    skip_vcf_filter = args.skip_vcf_filter and not args.enable_vcf_filter

    # 创建注释器并运行|Create annotator and run
    annotator = ANNOVARAnnotator(
        vcf_file=args.input,
        gff3_file=args.gff3,
        genome_file=args.genome,
        build_ver=args.build_ver,
        annovar_path=args.annovar_path,
        database_path=args.database_path,
        output_dir=args.output_dir,
        qual_threshold=args.qual_threshold,
        skip_gff_cleaning=args.skip_gff_cleaning,
        skip_gff_fix=args.skip_gff_fix,
        skip_vcf_filter=skip_vcf_filter,
        step=args.step
    )

    annotator.run_analysis()

if __name__ == "__main__":
    main()