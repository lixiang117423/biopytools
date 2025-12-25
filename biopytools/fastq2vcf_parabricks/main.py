"""
Fastq到VCF (Parabricks) 主程序模块 | Fastq to VCF (Parabricks) Main Module
"""

import argparse
import sys
import os
import time
from pathlib import Path

from .config import Fastq2VcfParabricksConfig
from .utils import Fastq2VcfLogger, CommandRunner, FileManager, SystemChecker
from .data_processing import QualityController, GenomeIndexer, ParabricksMapper, JointCaller, VariantFilter


class Fastq2VcfParabricksProcessor:
    """Fastq到VCF (Parabricks) 主处理器 | Main Fastq to VCF (Parabricks) Processor"""

    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = Fastq2VcfParabricksConfig(**kwargs)
        self.config.validate()

        # 创建日志目录 | Create log directory
        log_dir = os.path.join(self.config.project_base, "99.logs")
        Path(log_dir).mkdir(parents=True, exist_ok=True)

        # 初始化日志 | Initialize logging
        self.logger_manager = Fastq2VcfLogger(log_dir, self.config.verbose)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_dir, self.config.dry_run)

        # 初始化各个处理器 | Initialize processors
        self.quality_controller = QualityController(self.config, self.logger, self.cmd_runner)
        self.genome_indexer = GenomeIndexer(self.config, self.logger, self.cmd_runner)
        self.parabricks_mapper = ParabricksMapper(self.config, self.logger, self.cmd_runner)
        self.joint_caller = JointCaller(self.config, self.logger, self.cmd_runner)
        self.variant_filter = VariantFilter(self.config, self.logger, self.cmd_runner)

        # 记录开始时间 | Record start time
        self.pipeline_start_time = time.time()

    def step1_quality_control(self):
        """步骤1: 质量控制 | Step 1: Quality Control"""
        self.logger_manager.step("🧹 Step 1: 质量控制 | Quality Control")

        # 检查是否跳过质控 | Check if skipping QC
        if self.config.skip_qc:
            self.logger.info("⏭️ 用户指定跳过质控步骤 | User specified to skip QC step")
            return True

        return self.quality_controller.run_quality_control()

    def step2_build_genome_index(self):
        """步骤2: 构建基因组索引 | Step 2: Build Genome Index"""
        self.logger_manager.step("📊 Step 2: 构建基因组索引 | Build Genome Index")
        return self.genome_indexer.build_genome_index()

    def step3_sequence_mapping(self):
        """步骤3: 序列比对 | Step 3: Sequence Mapping"""
        self.logger_manager.step("🗺️ Step 3: 序列比对 | Sequence Mapping")
        return self.parabricks_mapper.run_parabricks_mapping()

    def step4_joint_calling(self):
        """步骤4: 联合变异检测 | Step 4: Joint Variant Calling"""
        self.logger_manager.step("🧬 Step 4: 联合变异检测 | Joint Variant Calling")
        success, vcf_path = self.joint_caller.run_joint_calling()
        if success:
            self.joint_caller.final_vcf_path = vcf_path
        return success, vcf_path

    def step5_variant_filtering(self, input_vcf: str = None):
        """步骤5: 变异过滤 | Step 5: Variant Filtering"""
        self.logger_manager.step("🧹 Step 5: 变异过滤 | Variant Filtering")

        if input_vcf is None:
            input_vcf = self.joint_caller.final_vcf_path

        if not input_vcf:
            self.logger.error("❌ 未找到输入VCF文件 | No input VCF file found")
            return False

        return self.variant_filter.filter_variants(input_vcf)

    def run_single_step(self, step_num: int):
        """运行单个步骤 | Run single step"""
        step_functions = {
            1: (self.step1_quality_control, "🧹 质量控制 | Quality Control"),
            2: (self.step2_build_genome_index, "📊 构建基因组索引 | Build Genome Index"),
            3: (self.step3_sequence_mapping, "🗺️ 序列比对 | Sequence Mapping"),
            4: (self.step4_joint_calling, "🧬 联合变异检测 | Joint Variant Calling"),
            5: (self.step5_variant_filtering, "🧹 变异过滤 | Variant Filtering")
        }

        if step_num not in step_functions:
            self.logger.error(f"❌ 无效的步骤编号 | Invalid step number: {step_num}")
            return False

        step_func, step_name = step_functions[step_num]
        self.logger.info(f"🚀 执行步骤 {step_num} | Executing step {step_num}: {step_name}")

        if step_num == 4:
            success, vcf_path = step_func()
            if success:
                self.logger.info(f"✅ 步骤 {step_num} 完成 | Step {step_num} completed: {step_name}")
                if vcf_path:
                    self.logger.info(f"📄 输出VCF: {vcf_path} | Output VCF: {vcf_path}")
            else:
                self.logger.error(f"❌ 步骤 {step_num} 失败 | Step {step_num} failed: {step_name}")
            return success
        else:
            success = step_func()
            if success:
                self.logger.info(f"✅ 步骤 {step_num} 完成 | Step {step_num} completed: {step_name}")
            else:
                self.logger.error(f"❌ 步骤 {step_num} 失败 | Step {step_num} failed: {step_name}")
            return success

    def run_full_pipeline(self):
        """运行完整的fastq到vcf流程 | Run complete fastq to vcf pipeline"""
        self.logger_manager.step("🎯 开始Fastq到VCF (Parabricks) 流程 | Starting Fastq to VCF (Parabricks) Pipeline")

        # 预检查 | Pre-flight checks
        self._pre_flight_checks()

        # 运行流程步骤 | Run pipeline steps
        steps = [
            (self.step1_quality_control, "🧹 质量控制 | Quality Control"),
            (self.step2_build_genome_index, "📊 构建基因组索引 | Build Genome Index"),
            (self.step3_sequence_mapping, "🗺️ 序列比对 | Sequence Mapping"),
            (self.step4_joint_calling, "🧬 联合变异检测 | Joint Variant Calling"),
            (self.step5_variant_filtering, "🧹 变异过滤 | Variant Filtering")
        ]

        final_vcf_path = None

        for i, (step_func, step_name) in enumerate(steps, 1):
            self.logger.info(f"🚀 执行步骤 {i} | Executing step {i}: {step_name}")

            if i == 4:  # 联合变异检测 | Joint variant calling
                success, vcf_path = step_func()
                if not success:
                    self.logger.error(f"❌ 步骤 {i} 失败 | Step {i} failed: {step_name}")
                    if vcf_path == "CLUSTER_MODE":  # 特殊返回值，表示集群模式
                        self.logger.info("📋 请按照生成的指南进行后续操作 | Please follow the generated guide for subsequent operations")
                        return True  # 这不算失败，只是需要手动操作
                    return False
                final_vcf_path = vcf_path
            else:
                success = step_func()
                if not success:
                    self.logger.error(f"❌ 步骤 {i} 失败 | Step {i} failed: {step_name}")
                    return False

            self.logger.info(f"✅ 步骤 {i} 完成 | Step {i} completed: {step_name}")

        # 生成最终报告 | Generate final report
        self._generate_final_report(final_vcf_path)

        self.logger.info("🎉 Fastq到VCF (Parabricks) 流程全部完成 | Fastq to VCF (Parabricks) pipeline completed!")
        return True

    def _pre_flight_checks(self):
        """预检查 | Pre-flight checks"""
        self.logger.info("🔍 系统预检查 | System pre-flight checks")

        # 检查必需工具 | Check required tools
        required_tools = ['bwa', 'samtools', 'bcftools', 'gatk', 'biopytools']
        for tool in required_tools:
            if not SystemChecker.check_command_exists(tool, self.logger):
                self.logger.error(f"❌ 缺少必需工具: {tool} | Missing required tool: {tool}")
                sys.exit(1)

        # 检查系统资源 | Check system resources
        SystemChecker.check_disk_space(self.config.project_base, 200, self.logger)
        SystemChecker.check_memory(64, self.logger)

        # 配置摘要 | Configuration summary
        self.logger.info("配置摘要 | Configuration Summary:")
        self.logger.info(f"  项目路径 | Project path: {self.config.project_base}")
        self.logger.info(f"  原始数据目录 | Raw data directory: {self.config.raw_fastq_dir}")
        self.logger.info(f"  参考基因组 | Reference genome: {self.config.ref_genome_fa}")
        self.logger.info(f"  线程配置 | Thread config: Mapping={self.config.threads_mapping}, Filtering={self.config.threads_filter}")
        self.logger.info(f"  断点续传 | Checkpoint: {self.config.enable_checkpoint}")
        self.logger.info(f"  测试模式 | Dry run: {self.config.dry_run}")

        self.logger.info("✅ 预检查通过 | Pre-flight checks passed")

    def _generate_final_report(self, final_vcf_path: str = None):
        """生成最终报告 | Generate final report"""
        report_file = os.path.join(self.config.project_base, "ANALYSIS_REPORT.txt")
        total_time = time.strftime('%H:%M:%S', time.gmtime(time.time() - self.pipeline_start_time))

        with open(report_file, 'w') as f:
            f.write("========================================================================\n")
            f.write("         Fastq到VCF (Parabricks) 分析流程 - 最终报告\n")
            f.write("         Fastq to VCF (Parabricks) Analysis Pipeline - Final Report\n")
            f.write("========================================================================\n")
            f.write(f"分析日期 | Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"项目路径 | Project Path: {self.config.project_base}\n")
            f.write(f"总运行时间 | Total Runtime: {total_time}\n\n")

            f.write("------------------------------------------------------------------------\n")
            f.write("📁 输入数据 | Input Data\n")
            f.write("------------------------------------------------------------------------\n")
            f.write(f"原始FASTQ目录 | Raw FASTQ Directory: {self.config.raw_fastq_dir}\n")
            f.write(f"参考基因组 | Reference Genome: {self.config.ref_genome_fa}\n")
            f.write(f"样本数量 | Sample Count: {FileManager.count_files(self.config.gvcf_dir, '*.g.vcf.gz')}\n\n")

            f.write("------------------------------------------------------------------------\n")
            f.write("⚙️ 处理参数 | Processing Parameters\n")
            f.write("------------------------------------------------------------------------\n")
            f.write(f"比对线程 | Mapping Threads: {self.config.threads_mapping}\n")
            f.write(f"过滤线程 | Filtering Threads: {self.config.threads_filter}\n\n")

            f.write("过滤参数 | Filtering Parameters:\n")
            f.write(f"  - SNP最小深度 | SNP Min Depth: {self.config.snp_min_dp}\n")
            f.write(f"  - SNP最小质量 | SNP Min Quality: {self.config.snp_min_qual}\n")
            f.write(f"  - InDel最小深度 | InDel Min Depth: {self.config.indel_min_dp}\n")
            f.write(f"  - InDel最小质量 | InDel Min Quality: {self.config.indel_min_qual}\n\n")

            f.write("------------------------------------------------------------------------\n")
            f.write("📂 输出目录 | Output Directories\n")
            f.write("------------------------------------------------------------------------\n")
            f.write(f"清洁数据 | Clean Data: {self.config.clean_fastq_dir}\n")
            f.write(f"比对结果 | Mapping Results: {self.config.mapping_dir}\n")
            f.write(f"变异检测 | Variant Calling: {self.config.joint_dir}\n")
            f.write(f"过滤结果 | Filtered Results: {self.config.filter_dir}\n\n")

            f.write("------------------------------------------------------------------------\n")
            f.write("📊 结果文件 | Result Files\n")
            f.write("------------------------------------------------------------------------\n")

            if final_vcf_path and os.path.exists(final_vcf_path):
                f.write(f"最终VCF | Final VCF: {final_vcf_path}\n")

            if os.path.exists(self.config.filter_dir):
                f.write("\n过滤后VCF文件 | Filtered VCF Files:\n")
                for vcf_file in FileManager.find_files(self.config.filter_dir, "*.vcf.gz"):
                    if os.path.exists(vcf_file):
                        size = FileManager.get_file_size(vcf_file)
                        f.write(f"  - {os.path.basename(vcf_file)}: {size}\n")

            f.write("\n========================================================================\n")
            f.write("              分析完成 - 感谢使用本流程\n")
            f.write("              Analysis Completed - Thank You for Using This Pipeline\n")
            f.write("========================================================================\n")

        self.logger.info(f"报告已生成 | Report generated: {report_file}")

    def run_analysis(self, step=None):
        """运行分析 | Run analysis"""
        try:
            if step:
                # 运行单个步骤 | Run single step
                success = self.run_single_step(step)
                if success:
                    self.logger.info(f"🎉 步骤 {step} 执行成功 | Step {step} executed successfully")
                else:
                    self.logger.error(f"💥 步骤 {step} 执行失败 | Step {step} execution failed")
                    sys.exit(1)
            else:
                # 运行完整流程 | Run full pipeline
                success = self.run_full_pipeline()
                if not success:
                    sys.exit(1)

        except Exception as e:
            self.logger.error(f"💥 程序执行出错 | Program execution error: {str(e)}")
            sys.exit(1)


def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 Fastq到VCF (Parabricks) 自动化分析脚本 | Fastq to VCF (Parabricks) Automation Script',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='📂 原始FASTQ文件目录路径 | Raw FASTQ files directory path')
    parser.add_argument('-r', '--ref-genome', required=True,
                       help='🧬 参考基因组文件路径 | Reference genome file path')
    parser.add_argument('-p', '--project-base', required=True,
                       help='🏗️ 项目根目录路径 | Project base directory path')

    # 可选参数 | Optional arguments
    parser.add_argument('--clean-fastq-dir',
                       help='🧹 清洁FASTQ文件目录路径 | Clean FASTQ files directory path')
    parser.add_argument('--mapping-dir',
                       help='🗺️ 比对结果目录路径 | Mapping results directory path')
    parser.add_argument('--gvcf-dir',
                       help='🧬 gVCF文件目录路径 | gVCF files directory path')
    parser.add_argument('--bam-dir',
                       help='📦 BAM文件目录路径 | BAM files directory path')
    parser.add_argument('--joint-dir',
                       help='🧬 联合检测结果目录路径 | Joint calling results directory path')
    parser.add_argument('--filter-dir',
                       help='🧹 过滤结果目录路径 | Filtering results directory path')
    parser.add_argument('-o', '--output-dir',
                       help='📁 输出目录路径 | Output directory path')

    # 线程配置 | Thread configuration
    parser.add_argument('--threads-mapping', type=int, default=88,
                       help='🔧 比对线程数 | Number of mapping threads')
    parser.add_argument('--threads-gtx', type=int, default=88,
                       help='🔧 GTX线程数 | Number of GTX threads')
    parser.add_argument('--threads-filter', type=int, default=88,
                       help='🔧 过滤线程数 | Number of filtering threads')

    # 过滤参数 | Filtering parameters
    parser.add_argument('--snp-min-dp', type=int, default=5,
                       help='🎯 SNP最小深度 | SNP minimum depth')
    parser.add_argument('--snp-min-qual', type=int, default=30,
                       help='🎯 SNP最小质量 | SNP minimum quality')
    parser.add_argument('--indel-min-dp', type=int, default=5,
                       help='🎯 InDel最小深度 | InDel minimum depth')
    parser.add_argument('--indel-min-qual', type=int, default=30,
                       help='🎯 InDel最小质量 | InDel minimum quality')

    # 样本阈值 | Sample thresholds
    parser.add_argument('--gatk-threshold', type=int, default=4,
                       help='🎯 GATK模式样本数阈值 (< N 使用GATK) | GATK sample count threshold')
    parser.add_argument('--gtx-single-threshold', type=int, default=200,
                       help='🎯 GTX单机模式样本数阈值 (< N 使用GTX单机) | GTX single machine sample count threshold')
    parser.add_argument('--gtx-window-size', type=int, default=20000000,
                       help='🎯 GTX分块窗口大小 (bp) | GTX chunk window size in bp')

    # 工具路径 | Tool paths
    parser.add_argument('--gtx-bin',
                       default='/share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx',
                       help='🛠️ GTX可执行文件路径 | GTX executable path')

    # 高级选项 | Advanced options
    parser.add_argument('--step', type=int, choices=[1, 2, 3, 4, 5],
                       help='🎯 只运行指定步骤 | Run only specified step '
                            '(1:🧹 质控 | QC, 2:📊 索引 | Index, 3:🗺️ 比对 | Mapping, '
                            '4:🧬 联合检测 | Joint calling, 5:🧹 过滤 | Filtering)')
    parser.add_argument('--no-checkpoint', action='store_true',
                       help='⏭️ 禁用断点续传 | Disable checkpoint resume')
    parser.add_argument('--dry-run', action='store_true',
                       help='🧪 测试模式，不执行实际命令 | Test mode, do not execute actual commands')
    parser.add_argument('--verbose', action='store_true',
                       help='📝 详细输出模式 | Verbose output mode')
    parser.add_argument('--skip-qc', action='store_true',
                       help='⏭️ 跳过质控步骤 | Skip QC step')
    parser.add_argument('--skip-mapping', action='store_true',
                       help='⏭️ 跳过比对步骤 | Skip mapping step')

    # 文件模式参数 | File pattern parameters
    parser.add_argument('--read1-pattern-fastp', default="_1.fq.gz",
                       help='📄 质控R1文件匹配模式 (默认: "_1.fq.gz") | QC R1 file pattern (default: "_1.fq.gz")')
    parser.add_argument('--read2-pattern-fastp', default="_2.fq.gz",
                       help='📄 质控R2文件匹配模式 (默认: "_2.fq.gz") | QC R2 file pattern (default: "_2.fq.gz")')

    args = parser.parse_args()

    # 创建处理器并运行 | Create processor and run
    processor = Fastq2VcfParabricksProcessor(
        raw_fastq_dir=args.input,
        ref_genome_fa=args.ref_genome,
        project_base=args.project_base,
        clean_fastq_dir=args.clean_fastq_dir,
        mapping_dir=args.mapping_dir,
        gvcf_dir=args.gvcf_dir,
        bam_dir=args.bam_dir,
        joint_dir=args.joint_dir,
        filter_dir=args.filter_dir,
        output_dir=args.output_dir,
        threads_mapping=args.threads_mapping,
        threads_gtx=args.threads_gtx,
        threads_filter=args.threads_filter,
        snp_min_dp=args.snp_min_dp,
        snp_min_qual=args.snp_min_qual,
        indel_min_dp=args.indel_min_dp,
        indel_min_qual=args.indel_min_qual,
        gatk_threshold=args.gatk_threshold,
        gtx_single_threshold=args.gtx_single_threshold,
        gtx_window_size=args.gtx_window_size,
        gtx_bin=args.gtx_bin,
        enable_checkpoint=not args.no_checkpoint,
        dry_run=args.dry_run,
        verbose=args.verbose,
        skip_qc=args.skip_qc,
        skip_mapping=args.skip_mapping,
        read1_pattern_fastp=args.read1_pattern_fastp,
        read2_pattern_fastp=args.read2_pattern_fastp
    )

    # 传递step参数给run_analysis方法 | Pass step parameter to run_analysis method
    processor.run_analysis(step=args.step)


if __name__ == "__main__":
    main()