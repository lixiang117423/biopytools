"""
二代转录组比对模块|Short-read (2nd gen) Alignment Module

HISAT2 索引构建 + 剪接感知比对
"""

import os
import logging
import tempfile
from typing import List, Dict, Optional, TYPE_CHECKING
from .config import RnaseqValConfig
from .utils import CommandRunner, FileValidator, build_conda_command

if TYPE_CHECKING:
    pass


class HISAT2Indexer:
    """HISAT2 索引构建器|HISAT2 Index Builder"""

    def __init__(self, config: RnaseqValConfig, logger: logging.Logger, cmd_runner: CommandRunner):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: logger 实例|Logger instance
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)

    def build_index(self) -> str:
        """构建 HISAT2 基因组索引（含已知剪接位点）|Build HISAT2 genome index with splice sites

        Returns:
            str: 索引前缀路径|Index prefix path
        """
        genome_fa = self.config.genome_fa
        annotation_gtf = self.config.annotation_gtf
        threads = self.config.threads
        output_dir = self.config.output_dir
        index_dir = os.path.join(output_dir, "01_index")
        os.makedirs(index_dir, exist_ok=True)

        genome_name = os.path.splitext(os.path.basename(genome_fa))[0]
        index_prefix = os.path.join(index_dir, f"{genome_name}.hisat2")

        # 断点续传：检查索引是否已存在|Checkpoint: check if index exists
        if os.path.exists(f"{index_prefix}.1.ht2"):
            self.logger.info(f"HISAT2 索引已存在，跳过|HISAT2 index exists, skipping: {index_prefix}")
            return index_prefix

        self.logger.step("步骤: 构建 HISAT2 索引|Step: Building HISAT2 index")

        # 提取剪接位点|Extract splice sites
        ss_file = os.path.join(index_dir, f"{genome_name}.ss")
        self._extract_splice_sites(annotation_gtf, ss_file)

        # 提取外显子|Extract exons
        exon_file = os.path.join(index_dir, f"{genome_name}.exon")
        self._extract_exons(annotation_gtf, exon_file)

        # 构建索引|Build index
        cmd_parts = ["hisat2-build", "-p", str(threads)]
        if os.path.exists(ss_file):
            cmd_parts.extend(["--ss", ss_file])
        if os.path.exists(exon_file):
            cmd_parts.extend(["--exon", exon_file])
        cmd_parts.extend([genome_fa, index_prefix])

        cmd = " ".join(cmd_parts)
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(cmd, "构建 HISAT2 索引|Building HISAT2 index")

        if not success:
            raise RuntimeError("HISAT2 索引构建失败|HISAT2 index build failed")

        return index_prefix

    def _extract_splice_sites(self, gtf_file: str, output_file: str) -> bool:
        """从 GTF 提取剪接位点|Extract splice sites from GTF

        Args:
            gtf_file: GTF 注释文件路径|GTF annotation file path
            output_file: 输出文件路径|Output file path

        Returns:
            bool: 成功返回 True|True if succeeded
        """
        if os.path.exists(output_file):
            return True

        cmd = f"hisat2_extract_splice_sites.py {gtf_file} > {output_file}"
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(cmd, "提取剪接位点|Extracting splice sites")
        if not success:
            self.logger.warning("剪接位点提取失败，将不使用剪接位点|Splice sites extraction failed, building without splice sites")
            if os.path.exists(output_file):
                os.remove(output_file)
        return success

    def _extract_exons(self, gtf_file: str, output_file: str) -> bool:
        """从 GTF 提取外显子|Extract exons from GTF

        Args:
            gtf_file: GTF 注释文件路径|GTF annotation file path
            output_file: 输出文件路径|Output file path

        Returns:
            bool: 成功返回 True|True if succeeded
        """
        if os.path.exists(output_file):
            return True

        cmd = f"hisat2_extract_exons.py {gtf_file} > {output_file}"
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(cmd, "提取外显子|Extracting exons")
        if not success:
            self.logger.warning("外显子提取失败，将不使用外显子|Exons extraction failed, building without exons")
            if os.path.exists(output_file):
                os.remove(output_file)
        return success


class HISAT2Aligner:
    """HISAT2 比对器|HISAT2 Aligner"""

    def __init__(self, config: RnaseqValConfig, logger: logging.Logger, cmd_runner: CommandRunner):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: logger 实例|Logger instance
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.file_validator = FileValidator(logger)

    def align_sample(self, sample: Dict, index_prefix: str) -> Optional[str]:
        """比对单个二代样本|Align a single short-read sample

        分两步执行避免管道中多个 conda run 串联：
        Step 1: hisat2 -> tmp SAM
        Step 2: samtools sort -> BAM

        Args:
            sample: 样本信息字典|Sample info dict (name, fastq1, fastq2)
            index_prefix: HISAT2 索引前缀|HISAT2 index prefix

        Returns:
            Optional[str]: 排序后 BAM 路径，失败返回 None|Sorted BAM path, None on failure
        """
        sample_name = sample["name"]
        fastq1 = sample["fastq1"]
        fastq2 = sample["fastq2"]
        threads = self.config.threads
        timeout = self.config.sample_timeout

        out_dir = os.path.join(self.config.output_dir, "02_align_2nd")
        os.makedirs(out_dir, exist_ok=True)
        bam_file = os.path.join(out_dir, f"{sample_name}.sorted.bam")
        flagstat_file = os.path.join(out_dir, f"{sample_name}.flagstat")

        # 断点续传|Checkpoint
        if os.path.exists(bam_file) and os.path.exists(bam_file + ".bai") and os.path.exists(flagstat_file):
            if not self.config.force:
                self.logger.info(f"样本比对已完成，跳过|Sample alignment done, skipping: {sample_name}")
                return bam_file

        self.logger.info(f"比对样本|Aligning sample: {sample_name}")

        # 构建 HISAT2 命令|Build HISAT2 command
        strand_arg = self.config.get_hisat2_strand_arg()
        hisat2_parts = [
            "hisat2",
            f"-x {index_prefix}",
            f"-1 {fastq1}",
            f"-2 {fastq2}",
            f"-p {threads}",
            "--dta",
        ]
        if strand_arg:
            hisat2_parts.append(strand_arg)

        # Step 1: hisat2 输出到临时 SAM 文件|Step 1: hisat2 -> tmp SAM
        tmp_sam = os.path.join(out_dir, f".{sample_name}.tmp.sam")
        hisat2_parts.extend(["-S", tmp_sam])
        hisat2_cmd = " ".join(hisat2_parts)
        hisat2_cmd = build_conda_command(hisat2_cmd, self.config.conda_env)

        success = self.cmd_runner.run(
            hisat2_cmd,
            f"HISAT2 比对|HISAT2 alignment: {sample_name}",
            timeout=timeout
        )

        if not success:
            self.logger.error(f"样本比对失败|Sample alignment failed: {sample_name}")
            return None

        # Step 2: samtools sort|Step 2: samtools sort
        sort_cmd = f"samtools sort -@ {threads} -O BAM -o {bam_file} {tmp_sam}"
        sort_cmd = build_conda_command(sort_cmd, self.config.conda_env)

        if not self.cmd_runner.run(sort_cmd, f"samtools sort: {sample_name}"):
            return None

        # 清理临时文件|Clean up temp file
        if os.path.exists(tmp_sam):
            os.remove(tmp_sam)

        # samtools index
        index_cmd = f"samtools index {bam_file}"
        index_cmd = build_conda_command(index_cmd, self.config.conda_env)
        if not self.cmd_runner.run(index_cmd, f"BAM 索引|BAM indexing: {sample_name}"):
            return None

        # samtools flagstat
        flagstat_cmd = f"samtools flagstat {bam_file} > {flagstat_file}"
        flagstat_cmd = build_conda_command(flagstat_cmd, self.config.conda_env)
        self.cmd_runner.run(flagstat_cmd, f"BAM 统计|BAM flagstat: {sample_name}")

        self.logger.info(f"样本比对完成|Sample alignment done: {sample_name} -> {bam_file}")
        return bam_file

    def align_all_samples(self, samples: List[Dict], index_prefix: str) -> List[str]:
        """比对所有二代样本|Align all short-read samples

        Args:
            samples: 样本列表|Sample list
            index_prefix: HISAT2 索引前缀|HISAT2 index prefix

        Returns:
            List[str]: 成功比对的 BAM 文件路径列表|List of successful BAM paths
        """
        bam_files = []
        for sample in samples:
            bam = self.align_sample(sample, index_prefix)
            if bam:
                bam_files.append(bam)
            else:
                self.logger.warning(f"跳过失败样本|Skipping failed sample: {sample['name']}")

        self.logger.info(f"二代比对完成|2nd gen alignment done: {len(bam_files)}/{len(samples)} 成功|succeeded")
        return bam_files
