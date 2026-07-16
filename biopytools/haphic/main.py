"""
HapHiC主程序模块|HapHiC Main Module
"""

import subprocess
import sys
import time
import os
import copy
from pathlib import Path
from typing import Optional

from .config import HapHiCConfig
from .utils import HapHiCLogger, QualityController, ResultValidator, BWAAligner, build_conda_command
from .pipeline import HapHiCPipeline
from ..common.paths import resolve_legacy_path


class HapHiCProcessor:
    """HapHiC处理器主类|Main HapHiC Processor Class"""

    def __init__(self, **kwargs):
        """初始化HapHiC处理器|Initialize HapHiC processor"""
        # 初始化配置|Initialize configuration
        self.config = HapHiCConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = HapHiCLogger(
            log_file=self.config.log_file,
            verbose=self.config.verbose
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化组件|Initialize components
        self.quality_controller = QualityController(self.logger, self.config)
        self.result_validator = ResultValidator(self.logger, self.config)

        # 初始化BWA比对器|Initialize BWA aligner
        if self.config.hic_file_type == "fastq":
            self.bwa_aligner = BWAAligner(self.config, self.logger)
        else:
            self.bwa_aligner = None

        self.logger.info("HapHiC处理器已初始化|HapHiC processor initialized")
        self._log_configuration()

    def _log_configuration(self):
        """记录配置信息|Log configuration information"""
        self.logger.info("配置参数|Configuration Parameters:")

        summary = self.config.get_summary()
        for key, value in summary.items():
            self.logger.info(f" {key}: {value}")

    def run_pipeline(self) -> bool:
        """运行完整流程|Run complete pipeline"""
        try:
            self.logger.info("开始HapHiC流程|Starting HapHiC pipeline")
            start_time = time.time()

            # 预检查|Pre-flight checks
            if not self._pre_flight_checks():
                return False

            # BWA比对 (如果需要)|BWA alignment (if needed)
            if self.config.hic_file_type == "fastq":
                # 检查是否已存在比对结果
                alignment_dir = resolve_legacy_path(self.config.output_dir, "00_mapping")
                potential_bam = os.path.join(alignment_dir, "HiC.bam")
                potential_filtered_bam = os.path.join(alignment_dir, "HiC.filtered.bam")

                self.logger.info(f"检查已存在的比对结果|Checking for existing alignment results")
                self.logger.debug(f"输出目录|Output directory: {self.config.output_dir}")
                self.logger.debug(f"比对目录|Alignment directory: {alignment_dir}")
                self.logger.debug(f"潜在BAM文件|Potential BAM file: {potential_bam}")
                self.logger.debug(f"潜在过滤BAM文件|Potential filtered BAM file: {potential_filtered_bam}")

                # 检查目录是否存在
                if os.path.exists(alignment_dir):
                    self.logger.debug(f"比对目录存在|Alignment directory exists")
                    try:
                        dir_contents = os.listdir(alignment_dir)
                        self.logger.debug(f"比对目录内容|Alignment directory contents: {dir_contents}")
                    except Exception as e:
                        self.logger.debug(f"无法列出比对目录内容|Cannot list alignment directory contents: {e}")
                else:
                    self.logger.debug(f"比对目录不存在|Alignment directory does not exist")

                # 更宽松的检查：即使没有索引文件，如果BAM文件存在且大小合理，也跳过比对
                if os.path.exists(potential_filtered_bam):
                    try:
                        size = os.path.getsize(potential_filtered_bam)
                        if size > 1000000:  # 至少1MB
                            self.logger.info(f"发现已存在的过滤BAM文件，跳过BWA比对步骤|Found existing filtered BAM file, skipping BWA alignment: {potential_filtered_bam} ({size:,} bytes)")
                            # HapHiC需要按read name排序的BAM文件，使用原始的过滤文件
                            # HapHiC requires read name-sorted BAM, use the original filtered file
                            self.logger.info(f"使用按read name排序的过滤BAM文件以满足HapHiC要求|Using read name-sorted filtered BAM file to meet HapHiC requirements")
                            self.config.bam_file = potential_filtered_bam
                            # HapHiC不需要BAM索引，跳过索引创建以避免错误
                            # HapHiC does not need BAM index, skip index creation to avoid errors
                        else:
                            self.logger.info(f"过滤BAM文件太小，重新执行比对|Filtered BAM file too small, re-running alignment: {size:,} bytes")
                    except Exception as e:
                        self.logger.warning(f"无法检查过滤BAM文件大小|Cannot check filtered BAM file size: {e}")
                elif os.path.exists(potential_bam):
                    try:
                        size = os.path.getsize(potential_bam)
                        if size > 1000000:  # 至少1MB
                            self.logger.info(f"发现已存在的BAM文件，跳过BWA比对步骤|Found existing BAM file, skipping BWA alignment: {potential_bam} ({size:,} bytes)")
                            # HapHiC需要按read name排序的BAM文件，不是coordinate排序的
                            # coordinate排序的文件通常是HiC.sorted.bam，我们需要使用原始的HiC.bam
                            # HapHiC requires read name-sorted BAM, not coordinate-sorted
                            # coordinate-sorted files are usually HiC.sorted.bam, we need to use the original HiC.bam
                            self.logger.info(f"使用按read name排序的BAM文件以满足HapHiC要求|Using read name-sorted BAM file to meet HapHiC requirements")
                            self.config.bam_file = potential_bam
                            # HapHiC不需要BAM索引，跳过索引创建以避免错误
                            # HapHiC does not need BAM index, skip index creation to avoid errors
                        else:
                            self.logger.info(f"BAM文件太小，重新执行比对|BAM file too small, re-running alignment: {size:,} bytes")
                    except Exception as e:
                        self.logger.warning(f"无法检查BAM文件大小|Cannot check BAM file size: {e}")
                else:
                    self.logger.info("执行BWA比对|Executing BWA alignment")
                    bam_file = self.bwa_aligner.run_alignment()
                    # 更新配置中的bam_file用于后续步骤|Update bam_file in config for subsequent steps
                    self.config.bam_file = bam_file
            else:
                # 对于BAM输入，检查索引文件是否存在
                if os.path.exists(self.config.hic_file + ".bai"):
                    self.logger.info(f"BAM索引已存在，跳过索引创建|BAM index already exists, skipping index creation")
                self.config.bam_file = self.config.hic_file

            # 创建流程管理器|Create pipeline manager
            pipeline = HapHiCPipeline(self.config, self.logger)

            # 运行流程|Run pipeline
            success = pipeline.run_complete_pipeline()

            # 检查是否有纠错后的contig，需要重新比对|Check if corrected contigs exist and need realignment
            if success and self.config.correct_nrounds > 0:
                corrected_asm = os.path.join(resolve_legacy_path(self.config.output_dir, "01_cluster"), "corrected_asm.fa")
                if os.path.exists(corrected_asm):
                    self.logger.info("检测到纠错后的contig，需要重新比对Hi-C数据|Detected corrected contigs, need to realign Hi-C data")

                    # 为纠错后的流程创建新的目录结构|Create new directory structure for corrected workflow
                    corrected_output_dir = resolve_legacy_path(self.config.output_dir, "07_corrected_contig_haphic")
                    os.makedirs(corrected_output_dir, exist_ok=True)
                    self.logger.info(f"创建纠错流程目录|Creating corrected workflow directory: {corrected_output_dir}")

                    # 创建纠错流程的配置|Create corrected workflow configuration
                    corrected_config = copy.deepcopy(self.config)
                    corrected_config.output_dir = corrected_output_dir
                    corrected_config.asm_file = corrected_asm
                    corrected_config.prefix = f"{self.config.prefix}_corrected"

                    # 重新比对Hi-C数据到纠错后的组装|Realign Hi-C data to corrected assembly
                    self.logger.info("开始重新比对Hi-C数据到纠错后的组装|Starting Hi-C realignment to corrected assembly")

                    # 检查是否已有纠错后的BAM文件|Check if corrected BAM already exists
                    corrected_alignment_dir = resolve_legacy_path(corrected_output_dir, "00_mapping_corrected")
                    potential_corrected_bam = os.path.join(corrected_alignment_dir, "HiC.filtered.bam")

                    if os.path.exists(potential_corrected_bam):
                        try:
                            size = os.path.getsize(potential_corrected_bam)
                            if size > 1000000:  # 至少1MB
                                self.logger.info(f"发现已存在的纠错BAM文件，跳过重新比对|Found existing corrected BAM file, skipping realignment: {potential_corrected_bam} ({size:,} bytes)")
                                corrected_config.bam_file = potential_corrected_bam
                            else:
                                self.logger.info(f"纠错BAM文件太小，重新执行比对|Corrected BAM file too small, re-running alignment: {size:,} bytes")
                                # 重新运行比对|Re-run alignment
                                self.bwa_aligner = BWAAligner(corrected_config, self.logger)
                                corrected_bam_file = self.bwa_aligner.run_alignment(is_realignment=True)
                                corrected_config.bam_file = corrected_bam_file
                        except Exception as e:
                            self.logger.warning(f"无法检查纠错BAM文件大小|Cannot check corrected BAM file size: {e}")
                            # 重新运行比对|Re-run alignment
                            self.bwa_aligner = BWAAligner(corrected_config, self.logger)
                            corrected_bam_file = self.bwa_aligner.run_alignment(is_realignment=True)
                            corrected_config.bam_file = corrected_bam_file
                    else:
                        # 重新运行比对|Re-run alignment
                        self.bwa_aligner = BWAAligner(corrected_config, self.logger)
                        corrected_bam_file = self.bwa_aligner.run_alignment(is_realignment=True)
                        corrected_config.bam_file = corrected_bam_file

                    # 创建纠错流程的pipeline|Create corrected workflow pipeline
                    corrected_pipeline = HapHiCPipeline(corrected_config, self.logger)

                    # 运行完整的HapHiC pipeline|Run complete HapHiC pipeline
                    self.logger.info("开始针对纠错后contig的完整HapHiC流程|Starting complete HapHiC workflow for corrected contigs")
                    corrected_success = corrected_pipeline.run_complete_pipeline()

                    if corrected_success:
                        self.logger.info("纠错后contig的HapHiC流程完成|Corrected contigs HapHiC workflow completed")
                    else:
                        self.logger.warning("纠错后contig的HapHiC流程失败|Corrected contigs HapHiC workflow failed")

                    # 更新总体成功状态|Update overall success status
                    success = corrected_success

            # 结果验证|Result validation
            if success:
                self._post_processing()

            # 统计信息|Statistics
            elapsed_time = time.time() - start_time
            self._log_statistics(elapsed_time)

            return success

        except Exception as e:
            self.logger.error(f"HapHiC流程执行失败|HapHiC pipeline execution failed: {e}")
            return False

    
    def _pre_flight_checks(self) -> bool:
        """预检查|Pre-flight checks"""
        self.logger.info("系统预检查|System pre-flight checks")

        # 验证输入文件|Validate input files
        if not self.quality_controller.validate_inputs():
            return False

        # 检查系统资源|Check system resources
        if not self.quality_controller.check_system_resources():
            self.logger.warning("系统资源检查有警告，但将继续执行|System resource check has warnings, but will continue")

        # 检查HapHiC工具|Check HapHiC tool
        # 如果是命令名（不是路径），跳过文件存在性检查|Skip file check if it's a command name
        if os.path.sep not in self.config.haphic_bin:
            # haphic_bin是命令名，尝试使用conda run验证|haphic_bin is command name, try conda run validation
            pass  # 跳过文件检查，直接进行conda run验证|Skip file check, go to conda run validation
        else:
            # haphic_bin是路径，检查文件是否存在|haphic_bin is path, check if file exists
            if not os.path.exists(self.config.haphic_bin):
                self.logger.error(f"HapHiC工具不存在|HapHiC tool not found: {self.config.haphic_bin}")
                return False

            # 检查工具可执行性|Check tool executable permission
            if not os.access(self.config.haphic_bin, os.X_OK):
                self.logger.error(f"HapHiC工具不可执行|HapHiC tool not executable: {self.config.haphic_bin}")
                return False

        # 验证工具可用性（使用conda run包装）|Verify tool availability (with conda run wrapper)
        try:
            # 传递完整路径，让build_conda_command自动检测conda环境
            # Pass full path, let build_conda_command auto-detect conda environment
            wrapped_cmd = build_conda_command(self.config.haphic_bin, ["--help"])

            result = subprocess.run(
                wrapped_cmd,
                capture_output=True,
                text=True,
                timeout=60,  # conda run可能需要更长时间启动，60秒足够|conda run may need longer, 60s is sufficient
                shell=False  # 列表形式必须使用shell=False|Must use shell=False with list
            )

            if result.returncode != 0:
                self.logger.error("HapHiC工具无法正常运行|HapHiC tool cannot run properly")
                if result.stderr:
                    self.logger.debug(f"错误信息|Error message: {result.stderr}")
                return False

            self.logger.info(f"HapHiC工具验证通过|HapHiC tool validated: {self.config.haphic_bin}")

        except subprocess.TimeoutExpired:
            self.logger.warning("HapHiC工具验证超时（可能是conda run启动慢），但将继续执行|HapHiC tool validation timeout (likely slow conda run), will continue")
            # 工具存在且可执行，只是验证超时，仍然继续
        except Exception as e:
            self.logger.error(f"HapHiC工具验证失败|HapHiC tool validation failed: {e}")
            return False

        self.logger.info("预检查通过|Pre-flight checks passed")
        return True

    def _post_processing(self):
        """后处理|Post processing"""
        self.logger.info("后处理和验证|Post-processing and validation")

        # 验证输出文件|Validate output files
        self.result_validator.validate_output_files()

        # 获取统计信息|Get statistics
        stats = self.result_validator.get_statistics()
        if stats:
            self.logger.info("输出统计|Output Statistics:")
            for key, value in stats.items():
                self.logger.info(f" {key}: {value}")

        # 生成报告|Generate report
        self._generate_report()

    def _log_statistics(self, elapsed_time: float):
        """记录统计信息|Log statistics"""
        self.logger.info("处理统计|Processing Statistics:")
        self.logger.info(f" 总耗时|Total time: {elapsed_time:.2f}秒")

        output_files = self.config.get_output_files()
        for file_type, file_path in output_files.items():
            if os.path.exists(file_path):
                size = os.path.getsize(file_path)
                self.logger.info(f" {file_type}: {os.path.basename(file_path)} ({size:,} bytes)")

    def _try_create_bam_index(self, bam_file: str):
        """已弃用：HapHiC不需要BAM索引|Deprecated: HapHiC does not need BAM index"""
        self.logger.info("HapHiC不需要BAM索引，跳过索引创建|HapHiC does not need BAM index, skipping index creation")
        pass

    def _generate_report(self):
        """生成报告|Generate report"""
        report_file = os.path.join(self.config.output_dir, "99_logs", f"{self.config.prefix}_haphic_report.txt")

        try:
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("HapHiC基因组Scaffolding分析报告|HapHiC Genome Scaffolding Report\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"分析时间|Analysis Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

                f.write("配置参数|Configuration Parameters:\n")
                f.write("-" * 30 + "\n")
                summary = self.config.get_summary()
                for key, value in summary.items():
                    f.write(f"{key}: {value}\n")

                f.write("\n输出文件|Output Files:\n")
                f.write("-" * 30 + "\n")
                output_files = self.config.get_output_files()
                for file_type, file_path in output_files.items():
                    if os.path.exists(file_path):
                        size = os.path.getsize(file_path)
                        f.write(f"{file_type}: {file_path} ({size:,} bytes)\n")

                f.write("\n" + "=" * 60 + "\n")
                f.write("报告生成完成|Report generation completed\n")

            self.logger.info(f"报告已生成|Report generated: {report_file}")

        except Exception as e:
            self.logger.error(f"报告生成失败|Report generation failed: {e}")


