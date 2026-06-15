"""
二代转录本组装模块|Short-read Transcript Assembly Module

StringTie 参考引导的转录本组装 + 多样本合并
"""

import os
import logging
from typing import List, Dict, Optional, TYPE_CHECKING

from .utils import CommandRunner, FileValidator, build_conda_command

if TYPE_CHECKING:
    from .config import RnaseqValConfig


class StringTieAssembler:
    """StringTie 转录本组装器|StringTie Transcript Assembler"""

    def __init__(self, config: 'RnaseqValConfig', logger: logging.Logger, cmd_runner: CommandRunner):
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

    def assemble_sample(self, sample_name: str, bam_file: str) -> Optional[str]:
        """StringTie 单样本组装|StringTie assemble a single sample

        Args:
            sample_name: 样本名|Sample name
            bam_file: 排序后的 BAM 文件路径|Sorted BAM file path

        Returns:
            Optional[str]: GTF 文件路径，失败返回 None|GTF file path, None on failure
        """
        out_dir = os.path.join(self.config.output_dir, "04_assemble_2nd")
        os.makedirs(out_dir, exist_ok=True)
        gtf_file = os.path.join(out_dir, f"{sample_name}.gtf")

        # 断点续传|Checkpoint
        if os.path.exists(gtf_file) and os.path.getsize(gtf_file) > 0:
            if not self.config.force:
                self.logger.info(f"样本组装已完成，跳过|Sample assembly done, skipping: {sample_name}")
                return gtf_file

        # 构建 StringTie 命令|Build StringTie command
        strand_flag = self.config.get_stringtie_strand_flag()

        cmd_parts = [
            "stringtie",
            bam_file,
            "-G", self.config.annotation_gtf,
            "-o", gtf_file,
            "-c", str(self.config.stringtie_min_cov),
            "-j", str(self.config.stringtie_min_junction_reads),
            "-f", str(self.config.stringtie_min_isoform_fraction),
            "-p", str(self.config.threads),
        ]
        if strand_flag:
            cmd_parts.append(strand_flag)

        cmd = " ".join(cmd_parts)
        cmd = build_conda_command(cmd, self.config.conda_env)

        success = self.cmd_runner.run(
            cmd,
            f"StringTie 组装|StringTie assembly: {sample_name}",
            timeout=self.config.sample_timeout
        )

        if not success:
            self.logger.error(f"样本组装失败|Sample assembly failed: {sample_name}")
            return None

        return gtf_file

    def merge_assemblies(self, gtf_list: List[str]) -> Optional[str]:
        """StringTie --merge 合并多样本组装结果|StringTie merge multi-sample assemblies

        Args:
            gtf_list: 单样本 GTF 文件路径列表|List of per-sample GTF file paths

        Returns:
            Optional[str]: 合并后的 GTF 路径，失败返回 None|Merged GTF path, None on failure
        """
        out_dir = os.path.join(self.config.output_dir, "04_assemble_2nd")
        merged_gtf = os.path.join(out_dir, "merged.gtf")
        gtf_list_file = os.path.join(out_dir, "gtf_list.txt")

        # 断点续传|Checkpoint
        if os.path.exists(merged_gtf) and os.path.getsize(merged_gtf) > 0:
            if not self.config.force:
                self.logger.info(f"合并 GTF 已存在，跳过|Merged GTF exists, skipping")
                return merged_gtf

        if not gtf_list:
            self.logger.warning("无有效 GTF 文件可用于合并|No valid GTF files to merge")
            return None

        self.logger.step("步骤: 合并二代组装结果|Step: Merging 2nd gen assemblies")

        # 写入 GTF 列表文件|Write GTF list file
        with open(gtf_list_file, "w") as f:
            for gtf in gtf_list:
                f.write(gtf + "\n")

        # 构建合并命令|Build merge command
        cmd_parts = [
            "stringtie",
            "--merge",
            "-G", self.config.annotation_gtf,
            "-o", merged_gtf,
            "-p", str(self.config.threads),
        ]
        # 添加每个 GTF|Add each GTF
        for gtf in gtf_list:
            cmd_parts.append(gtf)

        cmd = " ".join(cmd_parts)
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(cmd, "StringTie 合并|StringTie merge")

        if not success:
            self.logger.error("StringTie 合并失败|StringTie merge failed")
            return None

        self.logger.info(f"合并 GTF 完成|Merged GTF done: {merged_gtf}")
        return merged_gtf

    def run_all(self, bam_files: List[str], sample_names: Optional[List[str]] = None) -> Optional[str]:
        """执行全部二代样本组装和合并|Run all 2nd gen assembly and merge

        Args:
            bam_files: BAM 文件路径列表|BAM file path list
            sample_names: 对应的样本名列表（可选，默认从 BAM 文件名提取）|Sample name list

        Returns:
            Optional[str]: 合并后的 GTF 路径|Merged GTF path
        """
        if not bam_files:
            self.logger.warning("无二代 BAM 文件，跳过组装|No 2nd gen BAM files, skipping assembly")
            return None

        if sample_names is None:
            sample_names = []
            for bam in bam_files:
                base = os.path.splitext(os.path.basename(bam))[0]
                # 移除 .sorted 后缀|Remove .sorted suffix
                if base.endswith(".sorted"):
                    base = base[: -len(".sorted")]
                sample_names.append(base)

        self.logger.step("步骤: 二代转录本组装|Step: 2nd gen transcript assembly")
        self.logger.info(f"共 {len(bam_files)} 个样本|Total {len(bam_files)} samples")

        # 逐样本组装|Per-sample assembly
        gtf_files = []
        for bam, name in zip(bam_files, sample_names):
            gtf = self.assemble_sample(name, bam)
            if gtf:
                gtf_files.append(gtf)

        if not gtf_files:
            self.logger.error("所有样本组装均失败|All sample assemblies failed")
            return None

        # 合并|Merge
        merged = self.merge_assemblies(gtf_files)
        return merged
