"""
Fastq到VCF (GTX) 数据处理模块|Fastq to VCF (GTX) Data Processing Module
"""

import os
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple
import glob

from .config import Fastq2VcfGTXConfig
from .utils import CommandRunner, FileManager, CheckpointManager, Fastq2VcfGTXLogger


class QualityController:
    """质控处理器|Quality Control Processor"""

    def __init__(self, config: Fastq2VcfGTXConfig, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_quality_control(self) -> bool:
        """运行质量控制|Run quality control"""
        step_name = "quality_control"

        if self.config.enable_checkpoint and CheckpointManager(
            os.path.join(self.config.output_dir, "00_pipeline_info", "checkpoints"),
            self.logger
        ).exists(step_name):
            self.logger.info("检查点已存在，跳过质控|Checkpoint exists, skipping QC")
            return True

        if self.config.skip_qc:
            self.logger.info("用户指定跳过质控步骤|User specified to skip QC step")
            return True

        if self.config.raw_fastq_dir is None:
            self.logger.info("未提供原始FASTQ目录，跳过质控|No raw FASTQ directory provided, skipping QC")
            return True

        self.logger.info("开始质控处理|Starting quality control")

        # 统计原始文件数量|Count original files
        raw_count = FileManager.count_files(self.config.raw_fastq_dir, "*.fq.gz")
        if raw_count == 0:
            raw_count = FileManager.count_files(self.config.raw_fastq_dir, "*.fastq.gz")

        if raw_count == 0:
            self.logger.error("未找到原始 FASTQ 文件 (*.fq.gz or *.fastq.gz)|No raw FASTQ files found (*.fq.gz or *.fastq.gz)")
            return False

        self.logger.info(f"检测到 {raw_count} 个原始 FASTQ 文件|Found {raw_count} raw FASTQ files")

        # 确保输出目录存在|Ensure output directory exists
        FileManager.ensure_directory(self.config.clean_fastq_dir)

        # 运行质控命令|Run QC command
        command = [
            "biopytools", "fastp",
            "-i", self.config.raw_fastq_dir,
            "-o", self.config.clean_fastq_dir,
            "-t", str(self.config.threads),
            "--read1-suffix", self.config.read1_pattern_fastp,
            "--read2-suffix", self.config.read2_pattern_fastp
        ]

        success = self.cmd_runner.run(command, "质量控制|Quality Control")

        if success:
            clean_count = FileManager.count_files(self.config.clean_fastq_dir, "*.fq.gz")
            self.logger.info(f"质控完成: {clean_count} 个清洁文件|QC completed: {clean_count} clean files")

            # 创建检查点|Create checkpoint
            if self.config.enable_checkpoint:
                checkpoint_mgr = CheckpointManager(
                    os.path.join(self.config.output_dir, "00_pipeline_info", "checkpoints"),
                    self.logger
                )
                checkpoint_mgr.create(step_name)
        else:
            self.logger.error("质控处理失败|QC processing failed")

        return success


class GenomeIndexer:
    """基因组索引构建器|Genome Index Builder"""

    def __init__(self, config: Fastq2VcfGTXConfig, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def build_genome_index(self) -> bool:
        """构建基因组索引|Build genome index"""
        step_name = "genome_index"

        self.logger.info("准备执行索引构建步骤|Preparing to execute index building step")
        self.logger.info(f"当前配置中的基因组路径: {self.config.ref_genome_fa}|Genome path in current config")
        self.logger.info(f"GTX可执行文件路径: {self.config.gtx_bin}|GTX executable path")

        if self.config.enable_checkpoint and CheckpointManager(
            os.path.join(self.config.output_dir, "00_pipeline_info", "checkpoints"),
            self.logger
        ).exists(step_name):
            self.logger.info("检查点已存在，跳过索引构建|Checkpoint exists, skipping index building")
            return True

        self.logger.info("构建基因组索引|Building genome index")
        need_index = False

        # # 检查BWA索引|Check BWA index
        # if not os.path.exists(f"{self.config.ref_genome_fa}.bwt"):
        #     self.logger.info("构建 BWA 索引|Building BWA index")
        #     need_index = True
        #     if not self.cmd_runner.run(f"bwa index {self.config.ref_genome_fa}", "构建BWA索引|Build BWA index"):
        #         return False
        # else:
        #     self.logger.info("BWA 索引已存在|BWA index already exists")

        # # 检查SAMtools索引|Check SAMtools index
        # if not os.path.exists(f"{self.config.ref_genome_fa}.fai"):
        #     self.logger.info("构建 SAMtools 索引|Building SAMtools index")
        #     need_index = True
        #     if not self.cmd_runner.run(f"samtools faidx {self.config.ref_genome_fa}", "构建SAMtools索引|Build SAMtools index"):
        #         return False
        # else:
        #     self.logger.info("SAMtools 索引已存在|SAMtools index already exists")

        # 检查GTX索引|Check GTX index
        # 确认当前基因组路径|Confirm current genome path
        self.logger.info(f"当前参考基因组路径: {self.config.ref_genome_fa}|Current reference genome path")
        self.logger.info(f"文件是否存在: {'OK' if os.path.exists(self.config.ref_genome_fa) else 'MISSING'}|File exists")

        gtx_index_files = [
            # f"{self.config.ref_genome_fa}.gtx",
            # f"{self.config.ref_genome_fa}.bwt.2bit.64",
            f"{self.config.ref_genome_fa}.amb",
            f"{self.config.ref_genome_fa}.ann",
            f"{self.config.ref_genome_fa}.fai",
            f"{self.config.ref_genome_fa}.pac"
        ]

        gtx_index_exists = all(os.path.exists(f) for f in gtx_index_files)
        self.logger.info(f"检查GTX索引状态: {'已存在' if gtx_index_exists else '不存在'}|GTX index status: {'exists' if gtx_index_exists else 'missing'}")

        if not gtx_index_exists:
            self.logger.info("构建 GTX 索引|Building GTX index")
            need_index = True

            # 显示构建命令|Show build command
            gtx_index_cmd = f"faketime '2020-10-20 00:00:00' {self.config.gtx_bin} index {self.config.ref_genome_fa}"
            self.logger.info(f"将执行索引构建命令: {gtx_index_cmd}|Will execute index build command")

            # 使用faketime构建GTX索引|Build GTX index with faketime
            self.logger.info(f"正在执行GTX索引构建命令|Executing GTX index build command: {gtx_index_cmd}")
            if not self.cmd_runner.run(gtx_index_cmd, "构建GTX索引|Build GTX index"):
                self.logger.error("GTX 索引构建失败|GTX index building failed")

                # 构建失败后检查索引文件状态|Check index file status after build failure
                self.logger.info("检查构建后的索引文件状态:|Check index file status after build:")
                for idx_file in gtx_index_files:
                    exists = os.path.exists(idx_file)
                    self.logger.info(f"   {idx_file}: {'OK' if exists else 'MISSING'}")

                return False
            else:
                # 构建成功后验证索引文件|Verify index files after successful build
                self.logger.info("验证构建后的索引文件:|Verify index files after build:")
                for idx_file in gtx_index_files:
                    exists = os.path.exists(idx_file)
                    self.logger.info(f"   {idx_file}: {'OK' if exists else 'MISSING'}")
        else:
            self.logger.info("GTX 索引已存在|GTX index already exists")

        if not need_index:
            self.logger.info("所有索引均已存在|All indexes already exist")
        else:
            self.logger.info("索引构建完成|Index building completed")

        # 创建检查点|Create checkpoint
        if self.config.enable_checkpoint:
            checkpoint_mgr = CheckpointManager(
                os.path.join(self.config.output_dir, "00_pipeline_info", "checkpoints"),
                self.logger
            )
            checkpoint_mgr.create(step_name)

        return True


class GTXMapper:
    """GTX比对器|GTX Mapper"""

    def __init__(self, config: Fastq2VcfGTXConfig, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def refresh_config(self):
        """刷新配置，确保使用最新的路径|Refresh config to ensure using latest paths"""
        # 在使用前重新规范化路径|Re-normalize paths before use
        self.config.ref_genome_fa = os.path.normpath(os.path.abspath(self.config.ref_genome_fa))
        self.logger.info(f"刷新配置后的基因组路径: {self.config.ref_genome_fa}|Genome path after config refresh")

    def run_gtx_mapping(self) -> bool:
        """运行GTX比对|Run GTX mapping"""
        step_name = "mapping"

        # 刷新配置以确保使用最新路径|Refresh config to ensure latest paths
        self.refresh_config()

        if self.config.enable_checkpoint and CheckpointManager(
            os.path.join(self.config.output_dir, "00_pipeline_info", "checkpoints"),
            self.logger
        ).exists(step_name):
            self.logger.info("检查点已存在，跳过比对|Checkpoint exists, skipping mapping")
            return True

        if self.config.skip_mapping:
            self.logger.info("用户指定跳过比对步骤|User specified to skip mapping step")
            return True

        self.logger.info("使用GTX WGS进行序列比对 (CPU优化，比对+变异检测一体化)|Using GTX WGS for sequence mapping (CPU optimized, alignment+variant calling integrated)")
        self.logger.info(f"使用 {self.config.threads} 线程进行比对|Using {self.config.threads} threads for mapping")

        # 确保输出目录存在|Ensure output directory exists
        FileManager.ensure_directory(self.config.mapping_dir)
        FileManager.ensure_directory(self.config.gvcf_dir)
        FileManager.ensure_directory(self.config.bam_dir)

        if self.config.use_gtx_wgs:
            # 使用GTX WGS完整流程|Use GTX WGS complete pipeline
            success = self._run_gtx_wgs_pipeline()
        else:
            # 使用标准比对流程（用于兼容性）| Use standard mapping pipeline (for compatibility)
            success = self._run_standard_mapping()

        if success:
            gvcf_count = FileManager.count_files(self.config.gvcf_dir, "*.g.vcf.gz")
            bam_count = FileManager.count_files(self.config.bam_dir, "*.bam")
            self.logger.info(f"比对完成: {gvcf_count} 个 gVCF 文件, {bam_count} 个 BAM 文件|Mapping completed: {gvcf_count} gVCF files, {bam_count} BAM files")

            # 创建检查点|Create checkpoint
            if self.config.enable_checkpoint:
                checkpoint_mgr = CheckpointManager(
                    os.path.join(self.config.output_dir, "00_pipeline_info", "checkpoints"),
                    self.logger
                )
                checkpoint_mgr.create(step_name)
        else:
            self.logger.error("GTX比对失败|GTX mapping failed")

        return success

    def _run_gtx_wgs_pipeline(self) -> bool:
        """运行GTX WGS完整流程 (CPU优化，比对+变异检测一体化)|Run GTX WGS complete pipeline"""
        self.logger.info("使用GTX WGS完整流程 (CPU优化，比对+变异检测一体化)|Using GTX WGS complete pipeline (CPU optimized, alignment+variant calling integrated)")
        self.logger.info(f"使用 {self.config.threads} 线程处理|Using {self.config.threads} threads")

        # 根据输入来源决定查找目录和文件模式|Determine search directory and file patterns based on input source
        if self.config.raw_fastq_dir is None:
            # 用户显式提供 --clean-fastq-dir，从清洁数据目录查找
            # User explicitly provided --clean-fastq-dir, search in clean directory
            search_dir = self.config.clean_fastq_dir
            r1_patterns = ['*_1.clean.fq.gz', '*_1.fq.gz', '*_1.fastq.gz']
            r2_patterns = ['*_2.clean.fq.gz', '*_2.fq.gz', '*_2.fastq.gz']
            self.logger.info(f"使用清洁数据目录|Using clean data directory: {search_dir}")
        elif self.config.skip_qc:
            # 跳过质控，从原始目录查找（支持多种文件名格式）|Skip QC, search in raw directory
            search_dir = self.config.raw_fastq_dir
            r1_patterns = ['*_1.clean.fq.gz', '*_1.fq.gz', '*_1.fastq.gz']
            r2_patterns = ['*_2.clean.fq.gz', '*_2.fq.gz', '*_2.fastq.gz']
            self.logger.info(f"跳过质控模式，从原始目录查找FASTQ文件|Skip QC mode, searching for FASTQ files in raw directory: {search_dir}")
        else:
            # 使用质控后的数据|Use cleaned data
            search_dir = self.config.clean_fastq_dir
            r1_patterns = ['*_1.clean.fq.gz', '*_1.fq.gz']
            r2_patterns = ['*_2.clean.fq.gz', '*_2.fq.gz']
            self.logger.info(f"质控模式，从清洁数据目录查找文件|QC mode, searching for files in clean directory: {search_dir}")

        # 查找所有R1文件|Find all R1 files
        r1_files = []
        for pattern in r1_patterns:
            r1_files = FileManager.find_files(search_dir, pattern)
            if r1_files:
                break

        if not r1_files:
            self.logger.error(f"未找到任何R1文件|No R1 files found in {search_dir} with patterns: {r1_patterns}")
            return False

        total_samples = len(r1_files)
        self.logger.info(f"找到 {total_samples} 个样品需要处理|Found {total_samples} samples to process")

        current = 0
        failed_samples = []
        success_count = 0

        # 处理每个样品|Process each sample
        for r1_file in r1_files:
            current += 1

            # 提取样品名|Extract sample name
            sample_name = os.path.basename(r1_file)
            # 移除所有可能的R1后缀|Remove all possible R1 suffixes
            for suffix in ['_1.clean.fq.gz', '_1.fq.gz', '_1.fastq.gz']:
                if sample_name.endswith(suffix):
                    sample_name = sample_name[:-len(suffix)]
                    break

            # 构建R2文件路径|Build R2 file path
            r2_file = None
            # 直接构造R2文件名（不使用通配符）|Construct R2 filename directly (without wildcard)
            # 修复: 避免使用glob.glob()的模糊匹配，防止ER_1匹配到ER_10等文件
            for r2_suffix in ['_2.clean.fq.gz', '_2.fq.gz', '_2.fastq.gz']:
                potential_r2_file = os.path.join(search_dir, f"{sample_name}{r2_suffix}")
                if os.path.exists(potential_r2_file):
                    r2_file = potential_r2_file
                    break

            if not r2_file:
                self.logger.error(f"未找到样品 {sample_name} 的R2文件|R2 file not found for sample {sample_name}")
                self.logger.debug(f"尝试的模式|Tried patterns: {r2_patterns}")
                self.logger.debug(f"搜索目录|Search directory: {search_dir}")
                failed_samples.append(sample_name)
                continue

            # 解析软链接获取真实路径|Resolve symlinks to get real paths
            r1_file_real = os.path.realpath(r1_file)
            r2_file_real = os.path.realpath(r2_file)

            # 检查真实路径是否存在|Check if real paths exist
            if not os.path.exists(r1_file_real):
                self.logger.error(f"R1文件真实路径不存在|R1 file real path does not exist: {r1_file_real}")
                self.logger.debug(f"原始路径|Original path: {r1_file}")
                failed_samples.append(sample_name)
                continue

            if not os.path.exists(r2_file_real):
                self.logger.error(f"R2文件真实路径不存在|R2 file real path does not exist: {r2_file_real}")
                self.logger.debug(f"原始路径|Original path: {r2_file}")
                failed_samples.append(sample_name)
                continue

            # 定义输出文件|Define output files
            output_vcf = os.path.join(self.config.gvcf_dir, f"{sample_name}.g.vcf.gz")
            output_bam = os.path.join(self.config.bam_dir, f"{sample_name}.sorted.bam")

            # 检查是否已完成|Check if already completed
            if os.path.exists(output_vcf) and os.path.exists(output_bam):
                self.logger.info(f"[{current}/{total_samples}] 样品 {sample_name} 已处理，跳过|Sample {sample_name} already processed, skipping")
                success_count += 1
                continue

            self.logger.info(f"[{current}/{total_samples}] 处理样品: {sample_name}|Processing sample: {sample_name}")

            # 构建Read Group|Build Read Group
            read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tLB:{sample_name}"

            # 确保临时目录存在|Ensure temp directory exists
            tmp_dir = os.path.join(self.config.output_dir, '.tmp')
            FileManager.ensure_directory(tmp_dir)

            if not self.config.dry_run:
                # 调试信息|Debug info
                self.logger.info(f"GTX使用的基因组路径: {self.config.ref_genome_fa}|Genome path used by GTX")

                # 运行GTX WGS|Run GTX WGS
                command = (
                    f"faketime '2020-10-20 00:00:00' {self.config.gtx_bin} wgs "
                    f"-R \"{read_group}\" "
                    f"-o {output_vcf} "
                    f"-b {output_bam} "
                    f"-t {self.config.threads} "
                    f"-g "
                    f"--tmp-dir {tmp_dir} "
                    f"--pcr-indel-model {self.config.gtx_pcr_indel_model} "
                    f"--standard-min-confidence-threshold-for-calling {self.config.gtx_min_confidence} "
                    f"--min-base-quality-score {self.config.gtx_min_base_qual} "
                    f"--ploidy {self.config.gtx_ploidy} "
                    f"{self.config.ref_genome_fa} "
                    f"{r1_file_real} "  # 使用真实路径|Use real path
                    f"{r2_file_real}"  # 使用真实路径|Use real path
                )

                success = self.cmd_runner.run_with_progress(command, f"GTX WGS处理 {sample_name}|GTX WGS processing {sample_name}")

                if success:
                    self.logger.info(f"  样品 {sample_name} 完成|Sample {sample_name} completed")
                    success_count += 1

                    # 显示文件大小|Show file sizes
                    if os.path.exists(output_vcf):
                        vcf_size = FileManager.get_file_size(output_vcf)
                        self.logger.info(f"    VCF: {vcf_size}|VCF: {vcf_size}")

                    if os.path.exists(output_bam):
                        bam_size = FileManager.get_file_size(output_bam)
                        self.logger.info(f"    BAM: {bam_size}|BAM: {bam_size}")
                else:
                    self.logger.error(f"  样品 {sample_name} 处理失败|Sample {sample_name} processing failed")
                    failed_samples.append(sample_name)
            else:
                self.logger.info(f"  [DRY RUN] 跳过样品 {sample_name}|[DRY RUN] Skip sample {sample_name}")
                success_count += 1

        # 处理结果统计|Processing result statistics
        self.logger.info(f"GTX WGS处理完成|GTX WGS processing completed:")
        self.logger.info(f"  成功: {success_count}/{total_samples}|Success: {success_count}/{total_samples}")
        self.logger.info(f"  失败: {len(failed_samples)}/{total_samples}|Failed: {len(failed_samples)}/{total_samples}")

        if failed_samples:
            self.logger.warning("失败的样品|Failed samples:")
            for sample in failed_samples:
                self.logger.warning(f"  - {sample}")

        return len(failed_samples) == 0

    def _run_standard_mapping(self) -> bool:
        """标准比对流程（兼容性）| Standard mapping pipeline (compatibility)"""
        self.logger.info("使用标准比对模式（兼容性）| Using standard mapping mode (compatibility)")
        # 这里可以实现标准比对流程，但GTX版本主要使用GTX WGS
        # This could implement standard mapping, but GTX version mainly uses GTX WGS
        self.logger.warning("标准比对模式在GTX版本中不建议使用，建议设置use_gtx_wgs=True|Standard mapping mode not recommended in GTX version,建议设置use_gtx_wgs=True")
        return True


class JointCaller:
    """联合变异检测器|Joint Variant Caller"""

    def __init__(self, config: Fastq2VcfGTXConfig, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.final_vcf_path = ""

    def run_joint_calling(self) -> Tuple[bool, str]:
        """运行联合变异检测|Run joint variant calling"""
        step_name = "joint_calling"

        # 恢复检查点状态|Restore checkpoint status
        if self.config.enable_checkpoint:
            checkpoint_mgr = CheckpointManager(
                os.path.join(self.config.output_dir, "00_pipeline_info", "checkpoints"),
                self.logger
            )
            if checkpoint_mgr.exists(step_name):
                self.logger.info("检查点已存在，跳过联合检测|Checkpoint exists, skipping joint calling")

                # 尝试恢复VCF路径|Try to restore VCF path
                for vcf_name in ["gtx_joint_raw.vcf.gz", "joint_genotyping_raw.vcf.gz", "joint_genotyping_merged_filtered.vcf.gz"]:
                    vcf_path = os.path.join(self.config.joint_dir, vcf_name)
                    if os.path.exists(vcf_path):
                        self.final_vcf_path = vcf_path
                        return True, vcf_path

                return True, ""

        # 统计gVCF文件数量|Count gVCF files
        sample_count = FileManager.count_files(self.config.gvcf_dir, "*.g.vcf.gz")
        if sample_count == 0:
            self.logger.error("未找到任何 gVCF 文件|No gVCF files found")
            return False, ""

        self.logger.info(f"检测到 {sample_count} 个 gVCF 样本|Found {sample_count} gVCF samples")

        # 策略选择|Strategy selection
        self.logger.info("样本数分析策略选择|Sample count analysis for strategy selection:")
        self.logger.info(f"  < {self.config.gtx_single_threshold} → GTX单机模式|GTX single machine mode")
        self.logger.info(f"  >= {self.config.gtx_single_threshold} → GTX集群模式|GTX cluster mode")

        success = False
        if sample_count >= self.config.gtx_single_threshold:
            self.logger.warning("大规模样本模式，需要手动处理GTX集群任务|Large-scale sample mode, manual GTX cluster processing required")
            success, self.final_vcf_path = self._generate_gtx_cluster_scripts(sample_count)
        else:
            success, self.final_vcf_path = self._run_gtx_single_machine(sample_count)

        if success and self.final_vcf_path:
            # 创建检查点|Create checkpoint
            if self.config.enable_checkpoint:
                checkpoint_mgr = CheckpointManager(
                    os.path.join(self.config.output_dir, "00_pipeline_info", "checkpoints"),
                    self.logger
                )
                checkpoint_mgr.create(step_name)

        return success, self.final_vcf_path

    def _run_gtx_single_machine(self, sample_count: int) -> Tuple[bool, str]:
        """运行GTX单机模式|Run GTX single machine mode"""
        self.logger.info("使用 GTX 单机模式|Using GTX single machine mode")

        output_vcf = os.path.join(self.config.joint_dir, "gtx_joint_raw.vcf.gz")
        tmp_dir = os.path.join(self.config.output_dir, ".tmp", "gtx")
        FileManager.ensure_directory(tmp_dir)

        # 构建GTX命令|Build GTX command
        gtx_args = [
            f"-r {self.config.ref_genome_fa}",
            f"-o {output_vcf}",
            f"-t {self.config.threads}",
            f"--tmp-dir {tmp_dir}"
        ]

        # 添加gVCF文件|Add gVCF files
        gvcf_files = FileManager.find_files(self.config.gvcf_dir, "*.g.vcf.gz")
        for gvcf_file in gvcf_files:
            gtx_args.append(f"-v {gvcf_file}")

        command = f"faketime '2020-10-20 00:00:00' {self.config.gtx_bin} joint {' '.join(gtx_args)}"

        self.logger.info(f"准备处理 {len(gvcf_files)} 个样本|Preparing to process {len(gvcf_files)} samples")

        success = self.cmd_runner.run_with_progress(command, "GTX联合变异检测|GTX Joint Variant Calling", timeout=3600*24)  # 24小时超时

        if success and os.path.exists(output_vcf):
            self.final_vcf_path = output_vcf
            self.logger.info(f"GTX 输出: {output_vcf}|GTX output: {output_vcf}")
            return True, output_vcf
        else:
            self.logger.error("GTX未生成预期的VCF文件|GTX did not generate expected VCF file")
            return False, ""

    def _generate_gtx_cluster_scripts(self, sample_count: int) -> Tuple[bool, str]:
        """生成GTX集群脚本|Generate GTX cluster scripts"""
        self.logger.warning(f"大规模样本模式 (>= {self.config.gtx_single_threshold})|Large-scale sample mode")

        chunks_dir = os.path.join(self.config.joint_dir, "chunks")
        gtx_job_script = os.path.join(self.config.joint_dir, "01.run_gtx_jobs.sh")
        merge_py_script = os.path.join(self.config.output_dir, "00_pipeline_info", "scripts", "02.merge_vcf.py")
        final_merged_vcf = os.path.join(self.config.joint_dir, "merged_all.vcf.gz")

        FileManager.ensure_directory(chunks_dir)
        FileManager.ensure_directory(os.path.dirname(merge_py_script))

        # 生成合并脚本|Generate merge script
        self._generate_merge_script(merge_py_script)

        # 生成操作指南|Generate operation guide
        self._generate_manual_guide(sample_count, gtx_job_script, merge_py_script, final_merged_vcf)

        # 返回特殊状态码，表示需要手动处理|Return special status code indicating manual processing
        return False, ""

    def _generate_merge_script(self, script_path: str):
        """生成VCF合并脚本|Generate VCF merge script"""
        merge_script_content = f'''#!/usr/bin/env python3
"""
VCF合并脚本 - 支持自然排序和并行处理|VCF merge script - supports natural sorting and parallel processing
"""
import os
import sys
import glob
import re
import subprocess
import tempfile
from pathlib import Path

def natural_sort_key(filename):
    """自然排序关键字函数|Natural sort key function"""
    basename = os.path.basename(filename)
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'([0-9]+)', basename)]

def validate_vcf(vcf_file):
    """验证VCF文件完整性|Validate VCF file integrity"""
    try:
        result = subprocess.run(
            ['bcftools', 'index', '--nrecords', vcf_file],
            capture_output=True,
            text=True,
            check=True
        )
        return True
    except subprocess.CalledProcessError:
        print(f" 警告: {{vcf_file}} 验证失败|Warning: {{vcf_file}} validation failed", file=sys.stderr)
        return False

def main():
    if len(sys.argv) < 3:
        print("用法: python3 merge_vcf.py <input_dir> <output_vcf> [threads]", file=sys.stderr)
        print("Usage: python3 merge_vcf.py <input_dir> <output_vcf> [threads]", file=sys.stderr)
        sys.exit(1)

    input_dir = Path(sys.argv[1])
    output_file = Path(sys.argv[2])
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else {self.config.threads}

    print(f"使用 {{threads}} 个线程进行合并|Using {{threads}} threads for merging")

    # 查找VCF文件|Find VCF files
    vcf_pattern = input_dir / "*.joint.vcf.gz"
    vcf_files = sorted(glob.glob(str(vcf_pattern)), key=natural_sort_key)

    if not vcf_files:
        print(f"错误: 未找到 *.joint.vcf.gz 文件在 {{input_dir}}|Error: No *.joint.vcf.gz files found in {{input_dir}}", file=sys.stderr)
        sys.exit(1)

    print(f"发现 {{len(vcf_files)}} 个VCF文件|Found {{len(vcf_files)}} VCF files")

    # 验证VCF文件|Validate VCF files
    print("验证VCF文件完整性|Validating VCF file integrity...")
    valid_files = [f for f in vcf_files if validate_vcf(f)]

    if len(valid_files) != len(vcf_files):
        print(f" 警告: {{len(vcf_files) - len(valid_files)}} 个文件验证失败|Warning: {{len(vcf_files) - len(valid_files)}} files failed validation", file=sys.stderr)
        response = input("是否继续使用有效文件? Continue with valid files? (y/N): ")
        if response.lower() != 'y':
            sys.exit(1)
        vcf_files = valid_files

    print(f"{{len(vcf_files)}} 个文件验证通过|{{len(vcf_files)}} files validated")

    # 创建文件列表|Create file list
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
        for vcf in vcf_files:
            tmp.write(f"{{vcf}}\\n")
        list_path = tmp.name

    try:
        # 合并VCF|Merge VCF
        print(f"合并VCF文件到: {{output_file}}|Merging VCF files to: {{output_file}}")
        subprocess.check_call(
            f"bcftools concat -f {{list_path}} -a -O z -o {{output_file}} --threads {{threads}}",
            shell=True
        )

        # 创建索引|Create index
        print("创建索引|Creating index...")
        subprocess.check_call(f"tabix -p vcf {{output_file}}", shell=True)

        print(f"合并完成: {{output_file}}|Merge completed: {{output_file}}")

        # 显示统计信息|Show statistics
        result = subprocess.run(
            f"bcftools stats {{output_file}}|grep 'number of records:'",
            shell=True,
            capture_output=True,
            text=True
        )
        if result.stdout:
            print(f"{{result.stdout.strip()}}")

    except subprocess.CalledProcessError as e:
        print(f"合并失败: {{e}}|Merge failed: {{e}}", file=sys.stderr)
        sys.exit(1)
    finally:
        os.remove(list_path)

if __name__ == "__main__":
    main()
'''

        with open(script_path, 'w') as f:
            f.write(merge_script_content)

        os.chmod(script_path, 0o755)
        self.logger.info(f"合并脚本已生成: {script_path}|Merge script generated: {script_path}")

    def _generate_manual_guide(self, sample_count: int, gtx_job_script: str, merge_py_script: str, final_merged_vcf: str):
        """生成操作指南|Generate operation manual"""
        manual_guide = f'''

============================================================================
自动化流程已暂停 - 进入手动投递模式|Automated flow paused - Manual submission mode
============================================================================
样本数: {sample_count}|Sample count: {sample_count}
配置参数|Configuration parameters:
  - GTX单机阈值|GTX single threshold: {self.config.gtx_single_threshold}
  - GTX窗口大小|GTX window size: {self.config.gtx_window_size:,} bp

生成的脚本路径: Generated script paths:
  - VCF合并脚本: {merge_py_script}

操作步骤|Operation Steps:
----------------------------------------------------------------------------
 手动生成GTX分块任务|Manually generate GTX chunk jobs:
   请使用GTX命令生成脚本，窗口大小设置为 {self.config.gtx_window_size:,} bp
   Please use GTX command generation script with window size set to {self.config.gtx_window_size:,} bp

   参考命令|Reference command:
   bash /path/to/GTX_CMD_GEN_SCRIPT \\
       -g {self.config.gtx_bin} \\
       -r {self.config.ref_genome_fa} \\
       -i {self.config.gvcf_dir} \\
       -o {self.config.joint_dir}/chunks \\
       -w {self.config.gtx_window_size} \\
       -s {os.path.join(self.config.joint_dir, "01.run_gtx_jobs.sh")} \\
       -t {self.config.threads}

 投递GTX任务到集群|Submit GTX jobs to cluster:
   batch_sub -i {os.path.join(self.config.joint_dir, "01.run_gtx_jobs.sh")} \\
             -j gtx_joint \\
             -s 5 \\
             -m 800

3. 监控任务状态|Monitor job status:
   batch_stat -j gtx_joint

4. 任务完成后合并VCF|Merge VCF after jobs complete:
   python3 {merge_py_script} \\
           {os.path.join(self.config.joint_dir, "chunks")} \\
           {final_merged_vcf} \\
           {self.config.threads}

5. 验证合并结果|Validate merge result:
   bcftools stats {final_merged_vcf}|head -n 50

6. 运行变异过滤|Run variant filtering:
   biopytools filter-snp-indel \\
       -i {final_merged_vcf} \\
       -o {self.config.filter_dir} \\
       -t {self.config.threads} \\
       --snp-dp {self.config.snp_min_dp} \\
       --indel-dp {self.config.indel_min_dp}

============================================================================
提示|Tips:
  - 由于样本数 >= {self.config.gtx_single_threshold}，需要使用GTX集群模式
  - Since sample count >= {self.config.gtx_single_threshold}, GTX cluster mode is required
  - 窗口大小已设置为 {self.config.gtx_window_size:,} bp (约 {self.config.gtx_window_size/1000000:.1f} Mb)
  - Window size is set to {self.config.gtx_window_size:,} bp (~{self.config.gtx_window_size/1000000:.1f} Mb)
  - 建议先投递1-2个任务测试|Recommend submitting 1-2 jobs for testing first
  - 可用 tail -f 日志文件 查看日志|Use tail -f log_file to view logs
============================================================================
'''

        self.logger.info(manual_guide)


class VariantFilter:
    """变异过滤器|Variant Filter"""

    def __init__(self, config: Fastq2VcfGTXConfig, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def filter_variants(self, input_vcf: str) -> bool:
        """过滤变异|Filter variants"""
        step_name = "variant_filtering"

        if self.config.enable_checkpoint and CheckpointManager(
            os.path.join(self.config.output_dir, "00_pipeline_info", "checkpoints"),
            self.logger
        ).exists(step_name):
            self.logger.info("检查点已存在，跳过过滤|Checkpoint exists, skipping filtering")
            return True

        if not os.path.exists(input_vcf):
            self.logger.error(f"过滤输入文件不存在|Filtering input file does not exist: {input_vcf}")
            return False

        self.logger.info(f"输入 VCF: {input_vcf}|Input VCF: {input_vcf}")
        self.logger.info("过滤参数|Filtering parameters:")
        self.logger.info(f"  SNP  - 最小深度: {self.config.snp_min_dp}, 最小质量: {self.config.snp_min_qual}|SNP - min depth: {self.config.snp_min_dp}, min quality: {self.config.snp_min_qual}")
        self.logger.info(f"  InDel - 最小深度: {self.config.indel_min_dp}, 最小质量: {self.config.indel_min_qual}|InDel - min depth: {self.config.indel_min_dp}, min quality: {self.config.indel_min_qual}")

        # 确保输出目录存在|Ensure output directory exists
        FileManager.ensure_directory(self.config.filter_dir)

        command = (
            f"biopytools filter-snp-indel "
            f"-i {input_vcf} "
            f"-o {self.config.filter_dir} "
            f"-t {self.config.threads} "
            f"--snp-dp {self.config.snp_min_dp} "
            f"--indel-dp {self.config.indel_min_dp}"
        )

        success = self.cmd_runner.run_with_progress(command, "变异过滤|Variant Filtering")

        if success:
            self.logger.info("过滤完成|Filtering completed")

            # 显示结果统计|Show result statistics
            for vcf_file in FileManager.find_files(self.config.filter_dir, "*.vcf.gz"):
                if os.path.exists(vcf_file):
                    try:
                        count = subprocess.check_output(
                            f"bcftools view -H {vcf_file}|wc -l",
                            shell=True,
                            text=True
                        ).strip()
                        self.logger.info(f"  {os.path.basename(vcf_file)}: {count} 个变异|variants")
                    except Exception as e:
                        self.logger.warning(f"无法统计 {vcf_file}: {str(e)}|Failed to count {vcf_file}: {str(e)}")

            # 创建检查点|Create checkpoint
            if self.config.enable_checkpoint:
                checkpoint_mgr = CheckpointManager(
                    os.path.join(self.config.output_dir, "00_pipeline_info", "checkpoints"),
                    self.logger
                )
                checkpoint_mgr.create(step_name)
        else:
            self.logger.error("变异过滤失败|Variant filtering failed")

        return success