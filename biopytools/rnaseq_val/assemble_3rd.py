"""
三代转录本组装模块|Long-read Transcript Assembly Module

FLAIR (correct + collapse) 和 IsoQuant
输入已是 HiFi FASTQ，跳过 IsoSeq3 前处理
"""

import os
import logging
from typing import List, Dict, Optional, TYPE_CHECKING

from .utils import CommandRunner, FileValidator, build_conda_command

if TYPE_CHECKING:
    from .config import RnaseqValConfig


class FlairAssembler:
    """FLAIR 三代转录本组装器|FLAIR Long-read Transcript Assembler

    跳过 flair align（已有 minimap2 BAM），直接：
    1. BAM -> BED12
    2. flair correct（可用二代 junction 辅助）
    3. flair collapse
    """

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

    def run_sample(
        self,
        sample_name: str,
        bam_file: str,
        reads_file: str,
        sr_junctions: Optional[str] = None,
    ) -> Optional[str]:
        """运行 FLAIR 流程（correct + collapse）|Run FLAIR pipeline for a single sample

        Args:
            sample_name: 样本名|Sample name
            bam_file: 比对后的 BAM 文件|Aligned BAM file
            reads_file: 原始 reads FASTQ 文件|Original reads FASTQ file
            sr_junctions: 二代 junction BED 文件（可选）|SR junction BED file (optional)

        Returns:
            Optional[str]: 最终 GTF 路径，失败返回 None|Final GTF path, None on failure
        """
        out_dir = os.path.join(self.config.output_dir, "05_assemble_3rd", "flair", sample_name)
        os.makedirs(out_dir, exist_ok=True)
        final_gtf = os.path.join(out_dir, "collapse.isoforms.gtf")

        # 断点续传|Checkpoint
        if os.path.exists(final_gtf) and os.path.getsize(final_gtf) > 0:
            if not self.config.force:
                self.logger.info(f"FLAIR 组装已完成，跳过|FLAIR assembly done, skipping: {sample_name}")
                return final_gtf

        self.logger.info(f"运行 FLAIR 流程|Running FLAIR pipeline: {sample_name}")

        # Step 1: BAM -> BED12（使用 bedtools bamtobed 或 flair 内置转换）
        bed12_file = os.path.join(out_dir, "aligned.bed")
        if not self._bam_to_bed12(bam_file, bed12_file, sample_name):
            return None

        # Step 2: flair correct
        corrected_bed = os.path.join(out_dir, "corrected.bed")
        if not self._correct(bed12_file, corrected_bed, sample_name, sr_junctions):
            return None

        # Step 3: flair collapse
        if not self._collapse(corrected_bed, reads_file, out_dir, sample_name):
            return None

        self.logger.info(f"FLAIR 组装完成|FLAIR assembly done: {sample_name} -> {final_gtf}")
        return final_gtf

    def _bam_to_bed12(self, bam_file: str, bed12_file: str, sample_name: str) -> bool:
        """将 BAM 转换为 BED12 格式|Convert BAM to BED12 format

        Args:
            bam_file: BAM 文件路径|BAM file path
            bed12_file: 输出 BED12 文件路径|Output BED12 file path
            sample_name: 样本名|Sample name

        Returns:
            bool: 成功返回 True|True if succeeded
        """
        if os.path.exists(bed12_file):
            return True

        # 使用 bedtools bamtobed 转换
        cmd = f"bedtools bamtobed -i {bam_file} -bed12 > {bed12_file}"
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(cmd, f"BAM -> BED12: {sample_name}")

        if not success:
            self.logger.error(f"BAM -> BED12 转换失败|BAM -> BED12 conversion failed: {sample_name}")
            return False

        return True

    def _correct(
        self,
        bed12_file: str,
        corrected_bed: str,
        sample_name: str,
        sr_junctions: Optional[str] = None,
    ) -> bool:
        """flair correct：校正剪接位点|flair correct: splice site correction

        Args:
            bed12_file: 输入 BED12 文件|Input BED12 file
            corrected_bed: 输出校正后的 BED 文件|Output corrected BED file
            sample_name: 样本名|Sample name
            sr_junctions: 二代 junction 文件（可选）|SR junction file (optional)

        Returns:
            bool: 成功返回 True|True if succeeded
        """
        if os.path.exists(corrected_bed):
            return True

        cmd_parts = [
            "flair",
            "correct",
            "-q", bed12_file,
            "-g", self.config.genome_fa,
            "-f", self.config.annotation_gtf,
            "-o", corrected_bed,
            "--threads", str(self.config.threads),
        ]

        # 使用二代 junction 辅助校正（ONT 数据推荐）|Use SR junctions for correction (recommended for ONT)
        if sr_junctions and os.path.exists(sr_junctions):
            cmd_parts.extend(["--nvrna", sr_junctions])
            self.logger.info(f"使用二代 junction 辅助校正|Using SR junctions for correction")

        cmd = " ".join(cmd_parts)
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(cmd, f"FLAIR correct: {sample_name}")

        if not success:
            self.logger.error(f"FLAIR correct 失败|FLAIR correct failed: {sample_name}")
            return False

        return True

    def _collapse(
        self,
        corrected_bed: str,
        reads_file: str,
        out_dir: str,
        sample_name: str,
    ) -> bool:
        """flair collapse：聚类组装转录本|flair collapse: cluster and collapse transcripts

        Args:
            corrected_bed: 校正后的 BED 文件|Corrected BED file
            reads_file: 原始 reads FASTQ 文件|Original reads FASTQ file
            out_dir: 输出目录|Output directory
            sample_name: 样本名|Sample name

        Returns:
            bool: 成功返回 True|True if succeeded
        """
        collapse_prefix = os.path.join(out_dir, "collapse")

        cmd_parts = [
            "flair",
            "collapse",
            "-g", self.config.genome_fa,
            "-f", self.config.annotation_gtf,
            "-q", corrected_bed,
            "-r", reads_file,
            "-o", collapse_prefix,
            "--threads", str(self.config.threads),
        ]

        cmd = " ".join(cmd_parts)
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(cmd, f"FLAIR collapse: {sample_name}")

        if not success:
            self.logger.error(f"FLAIR collapse 失败|FLAIR collapse failed: {sample_name}")
            return False

        return True

    def _build_transcript_source_index(self, gtf_files: List[str]) -> Dict[str, str]:
        """建立 转录本坐标 -> 样本名 的索引|Build transcript coordinate -> sample name index

        Args:
            gtf_files: 各样本 GTF 文件路径列表|Per-sample GTF file path list

        Returns:
            Dict[str, str]: key 为 "chrom:start-end:strand"，value 为样本名
                Key is "chrom:start-end:strand", value is sample name
        """
        index = {}
        for gtf_path in gtf_files:
            sample_name = os.path.splitext(os.path.basename(gtf_path))[0]
            with open(gtf_path, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 9 or parts[2] != "transcript":
                        continue
                    key = f"{parts[0]}:{parts[3]}-{parts[4]}:{parts[6]}"
                    if key in index:
                        index[key] = f"{index[key]},{sample_name}"
                    else:
                        index[key] = sample_name
        return index

    def _annotate_support_samples(
        self, merged_gtf: str, source_index: Dict[str, str]
    ) -> None:
        """为合并后的 GTF 添加 support_sample 属性|Annotate merged GTF with support_sample

        根据转录本坐标反查来源样本，多样本支持的用逗号分隔

        Args:
            merged_gtf: 合并后的 GTF 文件路径|Merged GTF file path (in-place modification)
            source_index: 坐标 -> 样本名索引|Coordinate -> sample name index
        """
        annotated_file = merged_gtf + ".annotated_tmp"
        with open(merged_gtf, "r") as fin, open(annotated_file, "w") as fout:
            for line in fin:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    fout.write(line)
                    continue
                parts = stripped.split("\t")
                if len(parts) < 9:
                    fout.write(line)
                    continue
                if parts[2] == "transcript":
                    key = f"{parts[0]}:{parts[3]}-{parts[4]}:{parts[6]}"
                    sample = source_index.get(key, "")
                    if sample:
                        parts[8] = f'{parts[8].rstrip(";")}; support_sample "{sample}";'
                    fout.write("\t".join(parts) + "\n")
                else:
                    fout.write(line)
        # 原地替换|Replace in-place
        os.replace(annotated_file, merged_gtf)

    def _merge_assemblies(self, gtf_files: List[str]) -> Optional[str]:
        """合并多样本 FLAIR 组装结果并标注样本来源|Merge multi-sample FLAIR assemblies with source annotation

        使用 StringTie --merge 去重合并，合并后根据坐标反查为每条转录本标注 support_sample

        Args:
            gtf_files: 各样本 GTF 文件路径列表|Per-sample GTF file path list

        Returns:
            Optional[str]: 合并后的 GTF 路径，失败返回 None|Merged GTF path, None on failure
        """
        out_dir = os.path.join(self.config.output_dir, "05_assemble_3rd")
        merged_gtf = os.path.join(out_dir, "merged.gtf")

        # 断点续传|Checkpoint
        if os.path.exists(merged_gtf) and os.path.getsize(merged_gtf) > 0:
            if not self.config.force:
                self.logger.info(f"三代合并 GTF 已存在，跳过|LR merged GTF exists, skipping")
                return merged_gtf

        if not gtf_files:
            return None

        self.logger.step("步骤: 合并三代组装结果|Step: Merging 3rd gen assemblies")

        # Step 1: 建立坐标 -> 样本索引（在 merge 前构建，避免 merge 改写 transcript_id 后丢失信息）
        # Build coordinate -> sample index before merge (merge rewrites transcript_ids)
        source_index = self._build_transcript_source_index(gtf_files)
        self.logger.info(f"建立样本来源索引|Built source index: {len(source_index)} 条转录本坐标|transcript coordinates")

        # Step 2: StringTie --merge 去重合并|StringTie --merge for deduplication
        cmd_parts = [
            "stringtie",
            "--merge",
            "-G", self.config.annotation_gtf,
            "-o", merged_gtf,
            "-p", str(self.config.threads),
        ]
        for gtf in gtf_files:
            cmd_parts.append(gtf)

        cmd = " ".join(cmd_parts)
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(cmd, "StringTie 合并三代|StringTie merge 3rd gen")

        if not success:
            self.logger.error("三代 GTF 合并失败|LR GTF merge failed")
            return None

        # Step 3: 回写 support_sample 属性|Write back support_sample attribute
        if source_index:
            self._annotate_support_samples(merged_gtf, source_index)
            self.logger.info("已标注样本来源|Annotated support_sample attribute")

        self.logger.info(f"三代合并 GTF 完成|LR merged GTF done: {merged_gtf}")
        return merged_gtf

    def run_all(
        self,
        samples: List[Dict],
        bam_files: List[str],
        sr_junctions: Optional[str] = None,
    ) -> Optional[str]:
        """运行所有三代样本的 FLAIR 组装并合并|Run FLAIR assembly for all LR samples and merge

        多样本时使用 StringTie --merge 去重合并，并在 GTF attribute 中保留样本来源标注
        单样本时直接返回该样本的 GTF

        Args:
            samples: 样本信息列表|Sample info list (name, reads)
            bam_files: 对应的 BAM 文件列表|Corresponding BAM file list
            sr_junctions: 二代 junction BED 文件|SR junction BED file

        Returns:
            Optional[str]: 合并后的 GTF 路径（多样本）或单样本 GTF 路径
                Merged GTF path (multi-sample) or single sample GTF path, None on failure
        """
        if not samples:
            self.logger.warning("无三代样本，跳过 FLAIR|No LR samples, skipping FLAIR")
            return None

        self.logger.step("步骤: 三代转录本组装 (FLAIR)|Step: 3rd gen assembly (FLAIR)")
        self.logger.info(f"共 {len(samples)} 个样本|Total {len(samples)} samples")

        gtf_files = []
        for sample, bam in zip(samples, bam_files):
            gtf = self.run_sample(
                sample_name=sample["name"],
                bam_file=bam,
                reads_file=sample["reads"],
                sr_junctions=sr_junctions,
            )
            if gtf:
                gtf_files.append(gtf)

        self.logger.info(f"FLAIR 完成|FLAIR done: {len(gtf_files)}/{len(samples)} 成功|succeeded")

        if not gtf_files:
            return None

        # 多样本时合并去重|Merge and deduplicate for multiple samples
        if len(gtf_files) > 1:
            merged = self._merge_assemblies(gtf_files)
            if merged:
                return merged
            # 合并失败时回退：拼接不去重|Fallback: concatenate without dedup
            self.logger.warning("合并失败，回退为直接拼接|Merge failed, fallback to concatenation")
            return gtf_files[-1]

        return gtf_files[0]


class IsoQuantRunner:
    """IsoQuant 三代转录本定量与发现|IsoQuant LR transcript quantification and discovery"""

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

    def run(self, sample_name: str, bam_file: str) -> Optional[str]:
        """运行 IsoQuant|Run IsoQuant

        Args:
            sample_name: 样本名|Sample name
            bam_file: 比对后的 BAM 文件|Aligned BAM file

        Returns:
            Optional[str]: 输出目录路径，失败返回 None|Output directory path, None on failure
        """
        out_dir = os.path.join(self.config.output_dir, "05_assemble_3rd", "isoquant", sample_name)
        isoquant_gtf = os.path.join(out_dir, "isoquant.gtf")

        # 断点续传|Checkpoint
        if os.path.exists(isoquant_gtf):
            if not self.config.force:
                self.logger.info(f"IsoQuant 已完成，跳过|IsoQuant done, skipping: {sample_name}")
                return out_dir

        self.logger.info(f"运行 IsoQuant|Running IsoQuant: {sample_name}")

        # 根据 platform 确定参数|Determine data_type based on platform
        if self.config.lr_platform == "pacbio":
            data_type = "pacbio_ccs"
        else:
            data_type = "nanopore"

        cmd_parts = [
            "isoquant.py",
            "--reference", self.config.genome_fa,
            "--genedb", self.config.annotation_gtf,
            "--bam", bam_file,
            "--data_type", data_type,
            "--output", out_dir,
            "--threads", str(self.config.threads),
        ]

        cmd = " ".join(cmd_parts)
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(
            cmd,
            f"IsoQuant: {sample_name}",
            timeout=self.config.sample_timeout
        )

        if not success:
            self.logger.error(f"IsoQuant 失败|IsoQuant failed: {sample_name}")
            return None

        self.logger.info(f"IsoQuant 完成|IsoQuant done: {sample_name} -> {out_dir}")
        return out_dir
