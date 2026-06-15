"""
三代转录组比对模块|Long-read (3rd gen) Alignment Module

minimap2 剪接感知比对，支持 PacBio HiFi 和 ONT
"""

import os
import logging
from typing import List, Dict, Optional
from .config import RnaseqValConfig
from .utils import CommandRunner, FileValidator, build_conda_command


class Minimap2Aligner:
    """minimap2 三代比对器|minimap2 Long-read Aligner"""

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

    def align_sample(self, sample: Dict) -> Optional[str]:
        """比对单个三代样本|Align a single long-read sample

        根据平台自动选择 minimap2 参数：
        - pacbio: -ax splice:hq
        - ont: -ax splice

        分两步执行避免管道中多个 conda run 串联：
        Step 1: minimap2 -> tmp SAM
        Step 2: samtools sort -> BAM

        Args:
            sample: 样本信息字典|Sample info dict (name, reads)

        Returns:
            Optional[str]: 排序后 BAM 路径，失败返回 None|Sorted BAM path, None on failure
        """
        sample_name = sample["name"]
        reads = sample["reads"]
        threads = self.config.threads
        timeout = self.config.sample_timeout
        platform = self.config.lr_platform

        out_dir = os.path.join(self.config.output_dir, "03_align_3rd")
        os.makedirs(out_dir, exist_ok=True)
        bam_file = os.path.join(out_dir, f"{sample_name}.sorted.bam")
        flagstat_file = os.path.join(out_dir, f"{sample_name}.flagstat")

        # 断点续传|Checkpoint
        if os.path.exists(bam_file) and os.path.exists(bam_file + ".bai") and os.path.exists(flagstat_file):
            if not self.config.force:
                self.logger.info(f"三代比对已完成，跳过|3rd gen alignment done, skipping: {sample_name}")
                return bam_file

        self.logger.info(f"比对三代样本|Aligning LR sample: {sample_name} ({platform})")

        # 根据平台选择预设参数|Choose preset based on platform
        if platform == "pacbio":
            preset = "splice:hq"
        else:
            preset = "splice"

        # Step 1: minimap2 输出到临时 SAM|Step 1: minimap2 -> tmp SAM
        tmp_sam = os.path.join(out_dir, f".{sample_name}.tmp.sam")
        mm2_cmd = (
            f"minimap2 -ax {preset} -G {self.config.max_intron} --secondary=no "
            f"-t {threads} {self.config.genome_fa} {reads} -o {tmp_sam}"
        )
        mm2_cmd = build_conda_command(mm2_cmd, self.config.conda_env)

        success = self.cmd_runner.run(
            mm2_cmd,
            f"minimap2 比对|minimap2 alignment: {sample_name} ({platform})",
            timeout=timeout
        )

        if not success:
            self.logger.error(f"三代比对失败|3rd gen alignment failed: {sample_name}")
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

        self.logger.info(f"三代比对完成|3rd gen alignment done: {sample_name} -> {bam_file}")
        return bam_file

    def align_all_samples(self, samples: List[Dict]) -> List[str]:
        """比对所有三代样本|Align all long-read samples

        Args:
            samples: 样本列表|Sample list

        Returns:
            List[str]: 成功比对的 BAM 文件路径列表|List of successful BAM paths
        """
        bam_files = []
        for sample in samples:
            bam = self.align_sample(sample)
            if bam:
                bam_files.append(bam)
            else:
                self.logger.warning(f"跳过失败样本|Skipping failed sample: {sample['name']}")

        self.logger.info(f"三代比对完成|3rd gen alignment done: {len(bam_files)}/{len(samples)} 成功|succeeded")
        return bam_files
