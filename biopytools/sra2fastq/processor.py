"""
SRA文件处理模块 |SRA File Processing Module
"""

import os
import glob as glob_module
from pathlib import Path
from typing import List
from .utils import CommandRunner


class SRAProcessor:
    """SRA文件处理器|SRA File Processor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def _build_parallel_fastq_dump_cmd(self, sra_file: str) -> str:
        """构建parallel-fastq-dump命令|Build parallel-fastq-dump command"""
        cmd_parts = [self.config.tool_path]

        # SRA文件 (使用-s或--sra-id)|SRA file
        cmd_parts.append(f"--sra-id {sra_file}")

        # 线程数|Threads
        cmd_parts.append(f"--threads {self.config.threads}")

        # 输出目录|Output directory
        cmd_parts.append(f"--outdir {self.config.output_dir}")

        # 临时目录|Temporary directory
        if self.config.tmpdir:
            cmd_parts.append(f"--tmpdir {self.config.tmpdir}")

        # 以下参数会传递给底层的fastq-dump|Following params pass to underlying fastq-dump

        # 拆分双端测序|Split paired-end reads
        if self.config.split_files:
            cmd_parts.append("--split-files")

        # 压缩输出|Compress output
        if self.config.compress:
            cmd_parts.append("--gzip")

        # 跳过技术序列|Skip technical reads
        if self.config.skip_technical:
            cmd_parts.append("--skip-technical")

        # 剪切adapters|Clip adapters
        if self.config.clip:
            cmd_parts.append("--clip")

        # 最小读长过滤|Minimum read length filter
        if self.config.min_read_len > 0:
            cmd_parts.append(f"--minReadLen {self.config.min_read_len}")

        return " ".join(cmd_parts)

    def _build_fastq_dump_cmd(self, sra_file: str) -> str:
        """构建fastq-dump命令 (备选方案)|Build fastq-dump command (fallback)"""
        cmd_parts = [self.config.tool_path]

        # 拆分双端测序|Split paired-end reads
        if self.config.split_files:
            cmd_parts.append("--split-3")

        # 压缩输出|Compress output
        if self.config.compress:
            cmd_parts.append("--gzip")

        # 跳过技术序列|Skip technical reads
        if self.config.skip_technical:
            cmd_parts.append("--skip-technical")

        # 剪切adapters|Clip adapters
        if self.config.clip:
            cmd_parts.append("--clip")

        # 最小读长过滤|Minimum read length filter
        if self.config.min_read_len > 0:
            cmd_parts.append(f"--minReadLen {self.config.min_read_len}")

        # 输出目录|Output directory
        cmd_parts.append(f"--outdir {self.config.output_dir}")

        # 输入文件|Input file
        cmd_parts.append(sra_file)

        return " ".join(cmd_parts)

    def _is_already_converted(self, sra_file: str) -> bool:
        """
        检查SRA文件是否已转换（断点续传）|Check if SRA file already converted (checkpoint resume)

        通过检查输出目录中是否存在匹配的fastq文件来判断|Judge by checking for matching fastq files in output dir
        """
        base_name = Path(sra_file).stem
        output_dir = self.config.output_dir

        # 检查 .fq.gz 和 .fastq.gz 两种后缀|Check both .fq.gz and .fastq.gz suffixes
        for suffix in ['.fq.gz', '.fastq.gz']:
            # 单端: SRRxxx.fq.gz / 双端: SRRxxx_1.fq.gz, SRRxxx_2.fq.gz
            patterns = [
                os.path.join(output_dir, f"{base_name}{suffix}"),
                os.path.join(output_dir, f"{base_name}_1{suffix}"),
                os.path.join(output_dir, f"{base_name}_2{suffix}"),
                os.path.join(output_dir, f"{base_name}_pass{suffix}"),
            ]
            for pattern in patterns:
                if os.path.exists(pattern):
                    return True
                # 也检查 glob 匹配（处理 fastq-dump 可能添加的后缀）|Also check glob (handle suffixes fastq-dump may add)
                glob_pattern = os.path.join(output_dir, f"{base_name}*{suffix}")
                matches = glob_module.glob(glob_pattern)
                if matches:
                    return True

        return False

    def _rename_fastq_to_fq(self, sra_file: str) -> int:
        """
        将 .fastq.gz 重命名为 .fq.gz|Rename .fastq.gz to .fq.gz

        Returns:
            int: 重命名的文件数|Number of renamed files
        """
        base_name = Path(sra_file).stem
        output_dir = self.config.output_dir
        renamed_count = 0

        # glob 匹配该样本的所有 .fastq.gz 文件（单端/双端/带_pass等后缀）|Glob all .fastq.gz files for this sample
        glob_pattern = os.path.join(output_dir, f"{base_name}*.fastq.gz")
        for src in glob_module.glob(glob_pattern):
            dst = src.replace('.fastq.gz', '.fq.gz')
            os.rename(src, dst)
            self.logger.info(f"重命名|Renamed: {os.path.basename(src)} -> {os.path.basename(dst)}")
            renamed_count += 1

        return renamed_count

    def convert_single_file(self, sra_file: str) -> bool:
        """转换单个SRA文件|Convert single SRA file"""
        base_name = Path(sra_file).stem
        self.logger.info(f"{'='*60}")
        self.logger.info(f"处理文件|Processing file: {base_name}")
        self.logger.info(f"{'='*60}")

        # 断点续传：检查是否已完成|Checkpoint: skip if already done
        if self._is_already_converted(sra_file):
            self.logger.info(f"跳过已完成文件|Skipping already converted file: {base_name}")
            return True

        # 根据工具类型构建命令|Build command based on tool type
        if self.config.use_parallel:
            cmd = self._build_parallel_fastq_dump_cmd(sra_file)
            tool_name = "parallel-fastq-dump (多线程加速)|(multi-threaded)"
        else:
            cmd = self._build_fastq_dump_cmd(sra_file)
            tool_name = "fastq-dump (单线程)|(single-threaded)"

        self.logger.info(f"使用工具|Using tool: {tool_name}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")

        # 执行转换|Execute conversion
        success = self.cmd_runner.run(
            cmd,
            f"转换SRA文件|Converting SRA file: {base_name}"
        )

        if success:
            # 重命名 .fastq.gz -> .fq.gz|Rename .fastq.gz -> .fq.gz
            renamed = self._rename_fastq_to_fq(sra_file)
            if renamed > 0:
                self.logger.info(f"完成|Completed: {base_name} (重命名 {renamed} 个文件|renamed {renamed} files)")
            else:
                self.logger.info(f"完成|Completed: {base_name}")
        else:
            self.logger.error(f"失败|Failed: {base_name}")

        return success

    def convert_all_files(self) -> dict:
        """转换所有SRA文件|Convert all SRA files"""
        total = len(self.config.input_files)
        self.logger.info(f" 共找到 {total} 个SRA文件|Found {total} SRA files")

        results = {
            'success': [],
            'failed': [],
            'skipped': [],
            'total': total
        }

        for idx, sra_file in enumerate(self.config.input_files, 1):
            base_name = Path(sra_file).stem
            self.logger.info(f" 总进度|Overall Progress: [{idx}/{total}]")

            # 断点续传检查|Checkpoint check
            if self._is_already_converted(sra_file):
                self.logger.info(f"跳过已完成|Skipping already converted: {base_name}")
                results['success'].append(sra_file)
                results['skipped'].append(sra_file)
                # 即使是跳过的文件，也检查是否需要重命名（兼容旧文件）|Rename legacy .fastq.gz even for skipped files
                self._rename_fastq_to_fq(sra_file)
                continue

            if self.convert_single_file(sra_file):
                results['success'].append(sra_file)
            else:
                results['failed'].append(sra_file)

        return results
