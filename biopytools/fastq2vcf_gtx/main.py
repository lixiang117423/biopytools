"""
Fastq到VCF (GTX) 主程序模块|Fastq to VCF (GTX) Main Module
"""

import argparse
import sys
import os
import time
import subprocess
import shutil
import yaml
from pathlib import Path

from .config import Fastq2VcfGTXConfig
from .utils import Fastq2VcfGTXLogger, CommandRunner, FileManager, SystemChecker
from .data_processing import QualityController, GenomeIndexer, GTXMapper, JointCaller, VariantFilter


class Fastq2VcfGTXProcessor:
    """Fastq到VCF (GTX) 主处理器|Main Fastq to VCF (GTX) Processor"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = Fastq2VcfGTXConfig(**kwargs)
        self.config.validate()

        # 创建日志目录|Create log directory
        log_dir = os.path.join(self.config.output_dir, "99_logs")
        Path(log_dir).mkdir(parents=True, exist_ok=True)

        # 初始化日志|Initialize logging
        self.logger_manager = Fastq2VcfGTXLogger(log_dir, self.config.verbose, self.config.quiet)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_dir, self.config.dry_run)

        # 初始化各个处理器|Initialize processors
        self.quality_controller = QualityController(self.config, self.logger, self.cmd_runner)
        self.genome_indexer = GenomeIndexer(self.config, self.logger, self.cmd_runner)
        self.gtx_mapper = GTXMapper(self.config, self.logger, self.cmd_runner)
        self.joint_caller = JointCaller(self.config, self.logger, self.cmd_runner)
        self.variant_filter = VariantFilter(self.config, self.logger, self.cmd_runner)

        # 记录开始时间|Record start time
        self.pipeline_start_time = time.time()

    def step1_quality_control(self):
        """步骤1: 质量控制|Step 1: Quality Control"""
        self.logger_manager.step("Step 1: 质量控制|Quality Control")

        # 调试信息|Debug info
        self.logger.info(f"调试信息|Debug info: skip_qc = {self.config.skip_qc}")
        self.logger.info(f"调试信息|Debug info: read1_pattern_fastp = {self.config.read1_pattern_fastp}")
        self.logger.info(f"调试信息|Debug info: read2_pattern_fastp = {self.config.read2_pattern_fastp}")
        self.logger.info(f"调试信息|Debug info: raw_fastq_dir = {self.config.raw_fastq_dir}")

        # 检查是否跳过质控|Check if skipping QC
        if self.config.skip_qc:
            self.logger.info("用户指定跳过质控步骤|User specified to skip QC step")
            # 即使跳过QC也要执行GTX索引构建|Even if skipping QC, execute GTX index building
            self._force_build_gtx_index()
            return True

        # 执行质量控制|Execute quality control
        qc_success = self.quality_controller.run_quality_control()

        # fastp完成后立即执行GTX索引构建|Immediately execute GTX index building after fastp completion
        if qc_success:
            self.logger.info("Fastp质量控制完成，立即执行GTX索引构建|Fastp QC completed, immediately executing GTX index building")
            self._force_build_gtx_index()

        return qc_success

    def _force_build_gtx_index(self):
        """强制构建GTX索引|Force build GTX index"""
        self.logger_manager.step("强制构建GTX索引|Force Build GTX Index")

        # 确保基因组目录存在|Ensure genome directory exists
        genome_dir = self.config.genome_index_dir
        FileManager.ensure_directory(genome_dir)

        # 获取参考基因组文件名|Get reference genome filename
        genome_filename = os.path.basename(self.config.ref_genome_fa)
        genome_file_in_project = os.path.join(genome_dir, genome_filename)

        # 使用项目目录中的基因组文件|Use genome file in project directory
        if os.path.exists(genome_file_in_project):
            self.logger.info(f"使用项目目录中的基因组文件|Using genome file in project directory: {genome_file_in_project}")
            target_genome_file = genome_file_in_project
        else:
            self.logger.info(f"使用原始基因组文件|Using original genome file: {self.config.ref_genome_fa}")
            target_genome_file = self.config.ref_genome_fa

        # 构建GTX索引命令|Build GTX index command
        tmp_dir = os.path.join(self.config.output_dir, ".tmp")
        FileManager.ensure_directory(tmp_dir)

        gtx_index_cmd = f"faketime '2020-10-20 00:00:00' {self.config.gtx_bin} index {target_genome_file} --tmp-dir {tmp_dir}"

        self.logger.info(f"强制执行GTX索引构建|Force executing GTX index build")
        self.logger.info(f"命令|Command: {gtx_index_cmd}")
        self.logger.info(f"基因组文件|Genome file: {target_genome_file}")
        self.logger.info(f"临时目录|Temp directory: {tmp_dir}")

        # 执行索引构建|Execute index building
        if self.cmd_runner.run(gtx_index_cmd, "强制构建GTX索引|Force build GTX index"):
            self.logger.info("GTX索引构建成功|GTX index building successful")

            # 验证索引文件|Verify index files
            gtx_index_files = [
                f"{target_genome_file}.gtx",
                f"{target_genome_file}.gtx.bwt",
                f"{target_genome_file}.gtx.sa",
                f"{target_genome_file}.gtx.ann",
                f"{target_genome_file}.gtx.amb"
            ]

            self.logger.info("验证GTX索引文件|Verifying GTX index files:")
            for idx_file in gtx_index_files:
                exists = os.path.exists(idx_file)
                size = FileManager.get_file_size(idx_file) if exists else "0 B"
                self.logger.info(f"   {idx_file}: {'OK' if exists else 'MISSING'} ({size})")
        else:
            self.logger.error("GTX索引构建失败|GTX index building failed")

    def step2_build_genome_index(self):
        """步骤2: 构建基因组索引|Step 2: Build Genome Index"""
        self.logger.info("=" * 60)
        self.logger.info("构建基因组索引|Building Genome Index")
        self.logger.info("=" * 60)

        # GTX索引已在step1中构建，这里跳过|GTX index already built in step1, skip here
        self.logger.info("GTX索引已在step1中构建，跳过|GTX index already built in step1, skipping")
        return True

    def step3_sequence_mapping(self):
        """步骤3: 序列比对|Step 3: Sequence Mapping"""
        self.logger_manager.step("Step 3: 序列比对|Sequence Mapping")
        return self.gtx_mapper.run_gtx_mapping()

    def step4_joint_calling(self):
        """步骤4: 联合变异检测|Step 4: Joint Variant Calling"""
        self.logger_manager.step("Step 4: 联合变异检测|Joint Variant Calling")
        success, vcf_path = self.joint_caller.run_joint_calling()
        if success:
            self.joint_caller.final_vcf_path = vcf_path
        return success, vcf_path

    def step5_variant_filtering(self, input_vcf: str = None):
        """步骤5: 变异过滤|Step 5: Variant Filtering"""
        self.logger_manager.step("Step 5: 变异过滤|Variant Filtering")

        if input_vcf is None:
            input_vcf = self.joint_caller.final_vcf_path

        if not input_vcf:
            self.logger.error("未找到输入VCF文件|No input VCF file found")
            return False

        return self.variant_filter.filter_variants(input_vcf)

    def run_single_step(self, step_num: int):
        """运行单个步骤|Run single step"""
        step_functions = {
            1: (self.step1_quality_control, "质量控制|Quality Control"),
            2: (self.step2_build_genome_index, "构建基因组索引|Build Genome Index"),
            3: (self.step3_sequence_mapping, "序列比对|Sequence Mapping"),
            4: (self.step4_joint_calling, "联合变异检测|Joint Variant Calling"),
            5: (self.step5_variant_filtering, "变异过滤|Variant Filtering")
        }

        if step_num not in step_functions:
            self.logger.error(f"无效的步骤编号|Invalid step number: {step_num}")
            return False

        step_func, step_name = step_functions[step_num]
        self.logger.info(f"执行步骤{step_num}|Executing step {step_num}: {step_name}")

        if step_num == 4:
            success, vcf_path = step_func()
            if success:
                self.logger.info(f"步骤{step_num} 完成|Step {step_num} completed: {step_name}")
                if vcf_path:
                    self.logger.info(f"输出VCF: {vcf_path}|Output VCF: {vcf_path}")
            else:
                self.logger.error(f"步骤{step_num} 失败|Step {step_num} failed: {step_name}")
            return success
        else:
            success = step_func()
            if success:
                self.logger.info(f"步骤{step_num} 完成|Step {step_num} completed: {step_name}")
            else:
                self.logger.error(f"步骤{step_num} 失败|Step {step_num} failed: {step_name}")
            return success

    def run_full_pipeline(self):
        """运行完整的fastq到vcf流程|Run complete fastq to vcf pipeline"""
        self.logger_manager.step("开始Fastq到VCF (GTX) 流程|Starting Fastq to VCF (GTX) Pipeline")

        # 预检查|Pre-flight checks
        self._pre_flight_checks()

        # 运行流程步骤|Run pipeline steps
        steps = [
            (self.step1_quality_control, "质量控制|Quality Control"),
            (self.step2_build_genome_index, "构建基因组索引|Build Genome Index"),
            (self.step3_sequence_mapping, "序列比对|Sequence Mapping"),
            (self.step4_joint_calling, "联合变异检测|Joint Variant Calling"),
            (self.step5_variant_filtering, "变异过滤|Variant Filtering")
        ]

        final_vcf_path = None

        for i, (step_func, step_name) in enumerate(steps, 1):
            self.logger.info(f"执行步骤{i}|Executing step {i}: {step_name}")

            if i == 4:  # 联合变异检测|Joint variant calling
                success, vcf_path = step_func()
                if not success:
                    self.logger.error(f"步骤{i} 失败|Step {i} failed: {step_name}")
                    if vcf_path == "CLUSTER_MODE":  # 特殊返回值，表示集群模式
                        self.logger.info("请按照生成的指南进行后续操作|Please follow the generated guide for subsequent operations")
                        return True  # 这不算失败，只是需要手动操作
                    return False
                final_vcf_path = vcf_path
            else:
                success = step_func()
                if not success:
                    self.logger.error(f"步骤{i} 失败|Step {i} failed: {step_name}")
                    return False

            self.logger.info(f"步骤{i} 完成|Step {i} completed: {step_name}")

        # 生成最终报告|Generate final report
        self._generate_final_report(final_vcf_path)

        self.logger.info("Fastq到VCF (GTX) 流程全部完成|Fastq to VCF (GTX) pipeline completed!")
        return True

    def _update_all_genome_paths(self):
        """更新所有组件的参考基因组路径|Update reference genome paths for all components"""
        # 确保基因组文件在项目目录下|Ensure genome file is in project directory
        genome_dir = self.config.genome_index_dir
        target_genome_path = os.path.join(genome_dir, os.path.basename(self.config.ref_genome_fa))

        # 如果目标文件不存在，拷贝过去|Copy if target doesn't exist
        if not os.path.exists(target_genome_path):
            FileManager.ensure_directory(genome_dir)
            self.logger.info(f"拷贝基因组文件到 {target_genome_path}|Copying genome file to {target_genome_path}")
            import shutil
            shutil.copy2(self.config.ref_genome_fa, target_genome_path)

        # 强制更新配置中的路径|Force update path in config
        old_path = self.config.ref_genome_fa
        self.config.ref_genome_fa = os.path.normpath(os.path.abspath(target_genome_path))
        self.logger.info(f"强制更新基因组路径: {old_path} -> {self.config.ref_genome_fa}|Force updated genome path")

        # 更新所有组件的配置|Update all components' configs
        self.quality_controller.config.ref_genome_fa = self.config.ref_genome_fa
        self.genome_indexer.config.ref_genome_fa = self.config.ref_genome_fa
        self.gtx_mapper.config.ref_genome_fa = self.config.ref_genome_fa
        self.joint_caller.config.ref_genome_fa = self.config.ref_genome_fa
        self.variant_filter.config.ref_genome_fa = self.config.ref_genome_fa

        self.logger.info(f"已更新所有组件的基因组路径为: {self.config.ref_genome_fa}|Updated all components genome path to: {self.config.ref_genome_fa}")

    def _pre_flight_checks(self):
        """预检查|Pre-flight checks"""
        self.logger.info("系统预检查|System pre-flight checks")

        # 首先更新所有组件的基因组路径|First update genome paths for all components
        self.logger.info("开始更新所有组件的基因组路径...|Starting genome path update for all components...")
        self._update_all_genome_paths()
        self.logger.info("基因组路径更新完成|Genome path update completed")

        # 检查必需工具|Check required tools
        required_tools = ['samtools', 'bcftools', 'biopytools']
        for tool in required_tools:
            if not SystemChecker.check_command_exists(tool, self.logger):
                self.logger.error(f"缺少必需工具: {tool}|Missing required tool: {tool}")
                sys.exit(1)

        # 检查系统资源|Check system resources
        SystemChecker.check_disk_space(self.config.output_dir, 200, self.logger)
        SystemChecker.check_memory(64, self.logger)

        # 配置摘要|Configuration summary
        self.logger.info("配置摘要|Configuration Summary:")
        self.logger.info(f"  项目路径|Project path: {self.config.output_dir}")
        self.logger.info(f"  原始数据目录|Raw data directory: {self.config.raw_fastq_dir}")
        self.logger.info(f"  参考基因组|Reference genome: {self.config.ref_genome_fa}")
        self.logger.info(f"  线程配置|Thread config: {self.config.threads} threads")
        self.logger.info(f"  断点续传|Checkpoint: {self.config.enable_checkpoint}")
        self.logger.info(f"  测试模式|Dry run: {self.config.dry_run}")

        self.logger.info("预检查通过|Pre-flight checks passed")

        # 生成流程元数据|Generate pipeline metadata
        self._generate_pipeline_info()

    def _generate_pipeline_info(self):
        """生成流程元数据|Generate pipeline metadata"""
        info_dir = os.path.join(self.config.output_dir, "00_pipeline_info")
        FileManager.ensure_directory(info_dir)

        # 采集工具版本信息|Collect tool version information
        versions = {}
        tools_to_check = {
            'gtx': self.config.gtx_bin,
            'samtools': shutil.which('samtools') or 'samtools',
            'bcftools': shutil.which('bcftools') or 'bcftools',
        }
        for name, path in tools_to_check.items():
            try:
                r = subprocess.run(
                    [path, '--version'],
                    capture_output=True, text=True, timeout=10
                )
                versions[name] = {
                    'version': r.stdout.strip().split('\n')[0],
                    'path': path
                }
            except Exception:
                versions[name] = {'version': 'unknown', 'path': path}

        from . import __version__ as module_version

        info = {
            'pipeline': {
                'name': 'biopytools fastq2vcf_gtx',
                'version': module_version
            },
            'tools': versions,
            'parameters': {
                'threads': self.config.threads,
                'snp_min_dp': self.config.snp_min_dp,
                'snp_min_qual': self.config.snp_min_qual,
                'indel_min_dp': self.config.indel_min_dp,
                'indel_min_qual': self.config.indel_min_qual,
                'gtx_window_size': self.config.gtx_window_size,
                'gtx_single_threshold': self.config.gtx_single_threshold,
            },
        }
        with open(os.path.join(info_dir, 'software_versions.yml'), 'w') as f:
            yaml.dump(info, f, default_flow_style=False)

        self.logger.info(f"流程元数据已生成|Pipeline metadata generated: {info_dir}")

    def _generate_final_report(self, final_vcf_path: str = None):
        """生成最终报告|Generate final report"""
        report_file = os.path.join(self.config.output_dir, "ANALYSIS_REPORT.txt")
        total_time = time.strftime('%H:%M:%S', time.gmtime(time.time() - self.pipeline_start_time))

        with open(report_file, 'w') as f:
            f.write("========================================================================\n")
            f.write("         Fastq到VCF (GTX) 分析流程 - 最终报告\n")
            f.write("         Fastq to VCF (GTX) Analysis Pipeline - Final Report\n")
            f.write("========================================================================\n")
            f.write(f"分析日期|Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"项目路径|Project Path: {self.config.output_dir}\n")
            f.write(f"总运行时间|Total Runtime: {total_time}\n\n")

            f.write("------------------------------------------------------------------------\n")
            f.write(" 输入数据|Input Data\n")
            f.write("------------------------------------------------------------------------\n")
            f.write(f"原始FASTQ目录|Raw FASTQ Directory: {self.config.raw_fastq_dir}\n")
            f.write(f"参考基因组|Reference Genome: {self.config.ref_genome_fa}\n")
            f.write(f"样本数量|Sample Count: {FileManager.count_files(self.config.gvcf_dir, '*.g.vcf.gz')}\n\n")

            f.write("------------------------------------------------------------------------\n")
            f.write("处理参数|Processing Parameters\n")
            f.write("------------------------------------------------------------------------\n")
            f.write(f"线程数|Threads: {self.config.threads}\n\n")

            f.write("过滤参数|Filtering Parameters:\n")
            f.write(f"  - SNP最小深度|SNP Min Depth: {self.config.snp_min_dp}\n")
            f.write(f"  - SNP最小质量|SNP Min Quality: {self.config.snp_min_qual}\n")
            f.write(f"  - InDel最小深度|InDel Min Depth: {self.config.indel_min_dp}\n")
            f.write(f"  - InDel最小质量|InDel Min Quality: {self.config.indel_min_qual}\n\n")

            f.write("------------------------------------------------------------------------\n")
            f.write("输出目录|Output Directories\n")
            f.write("------------------------------------------------------------------------\n")
            f.write(f"清洁数据|Clean Data: {self.config.clean_fastq_dir}\n")
            f.write(f"比对结果|Mapping Results: {self.config.mapping_dir}\n")
            f.write(f"变异检测|Variant Calling: {self.config.joint_dir}\n")
            f.write(f"过滤结果|Filtered Results: {self.config.filter_dir}\n\n")

            f.write("------------------------------------------------------------------------\n")
            f.write("结果文件|Result Files\n")
            f.write("------------------------------------------------------------------------\n")

            if final_vcf_path and os.path.exists(final_vcf_path):
                f.write(f"最终VCF|Final VCF: {final_vcf_path}\n")

            if os.path.exists(self.config.filter_dir):
                f.write("\n过滤后VCF文件|Filtered VCF Files:\n")
                for vcf_file in FileManager.find_files(self.config.filter_dir, "*.vcf.gz"):
                    if os.path.exists(vcf_file):
                        size = FileManager.get_file_size(vcf_file)
                        f.write(f"  - {os.path.basename(vcf_file)}: {size}\n")

            f.write("\n========================================================================\n")
            f.write("              分析完成 - 感谢使用本流程\n")
            f.write("              Analysis Completed - Thank You for Using This Pipeline\n")
            f.write("========================================================================\n")

        self.logger.info(f"报告已生成|Report generated: {report_file}")

    def run_analysis(self, step=None):
        """运行分析|Run analysis"""
        try:
            if step:
                # 运行单个步骤|Run single step
                success = self.run_single_step(step)
                if success:
                    self.logger.info(f"步骤{step} 执行成功|Step {step} executed successfully")
                else:
                    self.logger.error(f"步骤{step} 执行失败|Step {step} execution failed")
                    sys.exit(1)
            else:
                # 运行完整流程|Run full pipeline
                success = self.run_full_pipeline()
                if not success:
                    sys.exit(1)

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='Fastq到VCF (GTX) 自动化分析脚本|Fastq to VCF (GTX) Automation Script',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数|Required arguments
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument('-i', '--input',
                         help='输入reads文件目录路径|Input reads directory path')
    required.add_argument('-g', '--genome', required=True,
                         help='参考基因组文件路径|Reference genome file path')
    required.add_argument('-o', '--output-dir', default='.',
                         help='输出目录路径|Output directory path')

    # 目录配置|Directory configuration
    directories = parser.add_argument_group('目录配置|Directory configuration')
    directories.add_argument('--clean-fastq-dir',
                            help='清洁FASTQ文件目录路径|Clean FASTQ files directory path')
    directories.add_argument('--mapping-dir',
                            help='比对结果目录路径|Mapping results directory path')
    directories.add_argument('--gvcf-dir',
                            help='gVCF文件目录路径|gVCF files directory path')
    directories.add_argument('--bam-dir',
                            help='BAM文件目录路径|BAM files directory path')
    directories.add_argument('--joint-dir',
                            help='联合检测结果目录路径|Joint calling results directory path')
    directories.add_argument('--filter-dir',
                            help='过滤结果目录路径|Filtering results directory path')

    # 线程配置|Thread configuration
    threads = parser.add_argument_group('线程配置|Thread configuration')
    threads.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Number of threads')

    # 过滤参数|Filtering parameters
    filtering = parser.add_argument_group('过滤参数|Filtering parameters')
    filtering.add_argument('--min-depth', '--snp-min-dp', type=int, default=5,
                         dest='snp_min_dp',
                         help='SNP/InDel最小深度|Minimum depth for SNP/InDel')
    filtering.add_argument('-q', '--quality', '--min-qual', '--snp-min-qual',
                         type=int, default=30,
                         dest='snp_min_qual',
                         help='SNP/InDel最小质量|Minimum quality for SNP/InDel')
    filtering.add_argument('--indel-min-dp', type=int, default=5,
                         help='InDel最小深度|InDel minimum depth')
    filtering.add_argument('--indel-min-qual', type=int, default=30,
                         help='InDel最小质量|InDel minimum quality')

    # 样本阈值|Sample thresholds
    thresholds = parser.add_argument_group('样本阈值|Sample thresholds')
    thresholds.add_argument('--gtx-single-threshold', type=int, default=200,
                          help='GTX单机模式样本数阈值|GTX single machine sample count threshold')
    thresholds.add_argument('--gtx-window-size', type=int, default=20000000,
                          help='GTX分块窗口大小|GTX chunk window size in bp')

    # 工具路径|Tool paths
    paths = parser.add_argument_group('工具路径|Tool paths')
    paths.add_argument('--gtx-bin',
                     default='~/software/gtx/bin/gtx',
                     help='GTX可执行文件路径|GTX executable path')
    paths.add_argument('--gtx-cmd-gen-script',
                     default="${HOME}/software/scripts/51.生成GTX按染色体合并gVCF的脚本.sh",
                     help='GTX命令生成脚本路径|GTX command generation script path')

    # GTX WGS参数|GTX WGS parameters
    gtx_params = parser.add_argument_group('GTX WGS参数|GTX WGS parameters')
    gtx_params.add_argument('--use-gtx-wgs', action='store_true', default=True,
                          help='使用GTX WGS|Use GTX WGS')
    gtx_params.add_argument('--no-gtx-wgs', dest='use_gtx_wgs', action='store_false',
                          help='禁用GTX WGS|Disable GTX WGS')
    gtx_params.add_argument('--gtx-pcr-indel-model',
                          default='CONSERVATIVE',
                          help='GTX PCR InDel模型|GTX PCR InDel model')
    gtx_params.add_argument('--gtx-min-confidence', type=int, default=30,
                          help='GTX最小置信度|GTX minimum confidence')
    gtx_params.add_argument('--gtx-min-base-qual', type=int, default=20,
                          help='GTX最小碱基质量|GTX minimum base quality')
    gtx_params.add_argument('--gtx-ploidy', type=int, default=2,
                          help='GTX倍性|GTX ploidy')

    # 执行控制|Execution control
    execution = parser.add_argument_group('执行控制|Execution control')
    execution.add_argument('-f', '--force', action='store_true',
                         help='强制覆盖已存在文件|Force overwrite existing files')
    execution.add_argument('--dry-run', action='store_true',
                         help='模拟运行|Dry run')
    execution.add_argument('--keep-intermediate', action='store_true',
                         help='保留中间文件|Keep intermediate files')
    execution.add_argument('--step', type=int, choices=[1, 2, 3, 4, 5],
                         help='只运行指定步骤|Run only specified step')
    execution.add_argument('--no-checkpoint', action='store_true',
                         help='禁用断点续传|Disable checkpoint resume')

    # 日志选项|Logging options
    logging_group = parser.add_argument_group('日志选项|Logging options')
    logging_group.add_argument('-v', '--verbose', action='count', default=0,
                            help='详细输出模式|Verbose mode')
    logging_group.add_argument('--quiet', action='store_true',
                            help='静默模式|Quiet mode')
    logging_group.add_argument('--log-file',
                            help='日志文件路径|Log file path')
    logging_group.add_argument('--log-level',
                            choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                            default='INFO',
                            help='日志级别|Log level')

    # 其他选项|Other options
    other = parser.add_argument_group('其他选项|Other options')
    other.add_argument('--skip-qc', action='store_true',
                      help='跳过质控步骤|Skip QC step')
    other.add_argument('--skip-mapping', action='store_true',
                      help='跳过比对步骤|Skip mapping step')
    other.add_argument('--read1-pattern-fastp', default="_1.fq.gz",
                      help='质控R1文件匹配模式|QC R1 file pattern')
    other.add_argument('--read2-pattern-fastp', default="_2.fq.gz",
                      help='质控R2文件匹配模式|QC R2 file pattern')

    args = parser.parse_args()

    # 处理verbose参数 (count -> bool)|Handle verbose parameter (count -> bool)
    verbose_bool = args.verbose > 0

    # 互斥验证: -i 和 --clean-fastq-dir 必须至少提供一个
    # Mutual exclusion: either -i or --clean-fastq-dir must be provided
    if not args.input and not args.clean_fastq_dir:
        parser.error("需要提供 --input/-i 或 --clean-fastq-dir|Must provide --input/-i or --clean-fastq-dir")

    # 当用户显式提供 --clean-fastq-dir 时，自动跳过质控
    # When user explicitly provides --clean-fastq-dir, auto-skip QC
    if args.clean_fastq_dir:
        raw_fastq_dir = None
        skip_qc = True
    else:
        raw_fastq_dir = args.input
        skip_qc = args.skip_qc

    # 创建处理器并运行|Create processor and run
    processor = Fastq2VcfGTXProcessor(
        raw_fastq_dir=raw_fastq_dir,
        ref_genome_fa=args.genome,
        output_dir=args.output_dir,
        clean_fastq_dir=args.clean_fastq_dir,
        mapping_dir=args.mapping_dir,
        gvcf_dir=args.gvcf_dir,
        bam_dir=args.bam_dir,
        joint_dir=args.joint_dir,
        filter_dir=args.filter_dir,
        threads=args.threads,
        snp_min_dp=args.snp_min_dp,
        snp_min_qual=args.snp_min_qual,
        indel_min_dp=args.indel_min_dp,
        indel_min_qual=args.indel_min_qual,
        gtx_single_threshold=args.gtx_single_threshold,
        gtx_window_size=args.gtx_window_size,
        use_gtx_wgs=args.use_gtx_wgs,
        gtx_pcr_indel_model=args.gtx_pcr_indel_model,
        gtx_min_confidence=args.gtx_min_confidence,
        gtx_min_base_qual=args.gtx_min_base_qual,
        gtx_ploidy=args.gtx_ploidy,
        gtx_bin=args.gtx_bin,
        gtx_cmd_gen_script=args.gtx_cmd_gen_script,
        enable_checkpoint=not args.no_checkpoint,
        dry_run=args.dry_run,
        verbose=verbose_bool,
        quiet=args.quiet,
        force=args.force,
        log_file=args.log_file,
        log_level=args.log_level,
        skip_qc=skip_qc,
        skip_mapping=args.skip_mapping,
        read1_pattern_fastp=args.read1_pattern_fastp,
        read2_pattern_fastp=args.read2_pattern_fastp
    )

    # 传递step参数给run_analysis方法|Pass step parameter to run_analysis method
    processor.run_analysis(step=args.step)


if __name__ == "__main__":
    main()