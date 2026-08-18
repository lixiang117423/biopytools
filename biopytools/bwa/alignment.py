"""
序列比对模块|Sequence Alignment Module
"""

import os
from pathlib import Path

from ..common.conda_runner import build_conda_command

from ..common.conda_runner import build_conda_command

from ..common.conda_runner import build_conda_command

from ..common.conda_runner import build_conda_command

from ..common.conda_runner import build_conda_command

from ..common.conda_runner import build_conda_command

class BWAAlignmentProcessor:
    """BWA比对处理器|BWA Alignment Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def align_sample(self, sample_name: str, read1: str, read2: str):
        """
        比对单个样品|Align single sample
        
        返回最终的BAM文件路径|Returns final BAM file path
        """
        self.logger.info(f"开始比对样品|Starting alignment for sample: {sample_name}")
        
        # 定义输出文件|Define output files
        sam_file = self.config.bam_dir / f"{sample_name}.sam"
        raw_bam = self.config.bam_dir / f"{sample_name}.raw.bam"
        sorted_bam = self.config.bam_dir / f"{sample_name}.sorted.bam"
        final_bam = self.config.bam_dir / f"{sample_name}.bam"
        
        # 检查断点续传|Check resume
        if self.config.resume and final_bam.exists():
            self.logger.info(f"跳过已完成的样品|Skipping completed sample: {sample_name}")
            return str(final_bam)
        
        # Step 1: BWA mem比对|BWA mem alignment
        if not self._run_bwa_mem(sample_name, read1, read2, sam_file):
            return None
        
        # Step 2: SAM转BAM|SAM to BAM
        if not self._sam_to_bam(sam_file, raw_bam):
            return None
        
        # Step 3: 排序BAM|Sort BAM
        if not self._sort_bam(raw_bam, sorted_bam):
            return None
        
        # Step 4: 标记重复序列（可选）| Mark duplicates (optional)
        if self.config.markdup:
            markdup_bam = self.config.bam_dir / f"{sample_name}.markdup.bam"
            if not self._mark_duplicates(sorted_bam, markdup_bam):
                return None
            final_bam = markdup_bam
        else:
            # 重命名为最终BAM|Rename to final BAM
            sorted_bam.rename(final_bam)
        
        # Step 5: 构建索引|Build index
        if not self._build_bam_index(final_bam):
            return None
        
        # Step 6: 清理中间文件|Cleanup intermediate files
        self._cleanup_intermediate_files(sam_file, raw_bam, sorted_bam)
        
        self.logger.info(f"样品比对完成|Sample alignment completed: {sample_name}")
        
        return str(final_bam)
    
    def _run_bwa_mem(self, sample_name: str, read1: str, read2: str, sam_file: Path) -> bool:
        """运行BWA mem|Run BWA mem"""
        bwa_options = self.config.get_bwa_cmd_options()

        # 添加read group信息|Add read group info
        rg = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"

        # 构建命令(conda环境自动包装, SAM输出重定向到文件)
        # |Build command (auto conda wrap, redirect SAM to file)
        args = ["mem"] + bwa_options.split() + [
            "-R", rg,
            self.config.genome,
            read1, read2,
        ]
        cmd = build_conda_command(self.config.bwa_path, args)

        success, _, _ = self.cmd_runner.run(
            cmd,
            f"BWA mem比对|BWA mem alignment: {sample_name}",
            output_file=str(sam_file.absolute()),
        )
        return success
    
    def _sam_to_bam(self, sam_file: Path, bam_file: Path) -> bool:
        """SAM转BAM|Convert SAM to BAM"""
        args = [
            "view",
            "-@", str(self.config.threads),
            "-b",
            "-o", str(bam_file.absolute()),
            str(sam_file.absolute()),
        ]
        cmd = build_conda_command(self.config.samtools_path, args)

        success, _, _ = self.cmd_runner.run(cmd, "SAM转BAM|SAM to BAM conversion")
        return success

    def _sort_bam(self, input_bam: Path, output_bam: Path) -> bool:
        """排序BAM|Sort BAM"""
        args = [
            "sort",
            "-@", str(self.config.threads),
            "-o", str(output_bam.absolute()),
            str(input_bam.absolute()),
        ]
        cmd = build_conda_command(self.config.samtools_path, args)

        success, _, _ = self.cmd_runner.run(cmd, "排序BAM|Sorting BAM")
        return success

    def _mark_duplicates(self, input_bam: Path, output_bam: Path) -> bool:
        """标记重复序列|Mark duplicates"""
        stats_file = self.config.stats_dir / f"{input_bam.stem}.markdup_stats.txt"

        args = [
            "markdup",
            "-@", str(self.config.threads),
        ]
        if self.config.remove_dup:
            args.append("-r")
        args.extend([
            "-f", str(stats_file.absolute()),
            str(input_bam.absolute()), str(output_bam.absolute()),
        ])
        cmd = build_conda_command(self.config.samtools_path, args)

        success, _, _ = self.cmd_runner.run(cmd, "标记重复序列|Marking duplicates")
        return success

    def _build_bam_index(self, bam_file: Path) -> bool:
        """构建BAM索引|Build BAM index"""
        args = ["index", "-@", str(self.config.threads), str(bam_file.absolute())]
        cmd = build_conda_command(self.config.samtools_path, args)

        success, _, _ = self.cmd_runner.run(cmd, "构建BAM索引|Building BAM index")
        return success
    
    def _cleanup_intermediate_files(self, sam_file: Path, raw_bam: Path, sorted_bam: Path):
        """清理中间文件|Cleanup intermediate files"""
        self.logger.info("清理中间文件|Cleaning up intermediate files")

        # 删除SAM文件|Remove SAM file
        if not self.config.keep_sam and sam_file.exists():
            sam_file.unlink()
            self.logger.info(f"已删除SAM文件|Removed SAM file: {sam_file.name}")

        # 删除原始BAM|Remove raw BAM
        if raw_bam.exists():
            raw_bam.unlink()
            self.logger.info(f"已删除原始BAM|Removed raw BAM: {raw_bam.name}")

        # 如果做了markdup，删除排序BAM|Remove sorted BAM if markdup was done
        if self.config.markdup and sorted_bam.exists():
            sorted_bam.unlink()
            self.logger.info(f"已删除排序BAM|Removed sorted BAM: {sorted_bam.name}")
