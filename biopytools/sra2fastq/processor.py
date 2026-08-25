"""
SRA文件处理模块 |SRA File Processing Module
"""

import os
import glob as glob_module
from pathlib import Path
from typing import List
from ..common.conda_runner import build_conda_command
from .utils import CommandRunner


class SRAProcessor:
    """SRA文件处理器|SRA File Processor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def _build_parallel_fastq_dump_args(self, sra_file: str) -> list:
        """构建parallel-fastq-dump参数|Build parallel-fastq-dump arguments"""
        args = []

        # SRA文件 (使用-s或--sra-id)|SRA file
        args.extend(['--sra-id', sra_file])

        # 线程数|Threads
        args.extend(['--threads', str(self.config.threads)])

        # 输出目录|Output directory
        args.extend(['--outdir', self.config.output_dir])

        # 临时目录|Temporary directory
        if self.config.tmpdir:
            args.extend(['--tmpdir', self.config.tmpdir])

        # 以下参数会传递给底层的fastq-dump|Following params pass to underlying fastq-dump

        # 拆分双端测序|Split paired-end reads
        if self.config.split_files:
            args.append("--split-files")

        # 压缩输出|Compress output
        if self.config.compress:
            args.append("--gzip")

        # 跳过技术序列|Skip technical reads
        if self.config.skip_technical:
            args.append("--skip-technical")

        # 剪切adapters|Clip adapters
        if self.config.clip:
            args.append("--clip")

        # 最小读长过滤|Minimum read length filter
        if self.config.min_read_len > 0:
            args.extend(['--minReadLen', str(self.config.min_read_len)])

        return args

    def _build_fastq_dump_args(self, sra_file: str) -> list:
        """构建fastq-dump参数 (备选方案)|Build fastq-dump arguments (fallback)"""
        args = []

        # 拆分双端测序|Split paired-end reads
        if self.config.split_files:
            args.append("--split-3")

        # 压缩输出|Compress output
        if self.config.compress:
            args.append("--gzip")

        # 跳过技术序列|Skip technical reads
        if self.config.skip_technical:
            args.append("--skip-technical")

        # 剪切adapters|Clip adapters
        if self.config.clip:
            args.append("--clip")

        # 最小读长过滤|Minimum read length filter
        if self.config.min_read_len > 0:
            args.extend(['--minReadLen', str(self.config.min_read_len)])

        # 输出目录|Output directory
        args.extend(['--outdir', self.config.output_dir])

        # 输入文件|Input file
        args.append(sra_file)

        return args

    # 输出后缀: 压缩与未压缩都覆盖(断点续传对 --no-compress 同样生效)
    # |Output suffixes: cover both compressed and uncompressed (checkpoint works with --no-compress too)
    OUTPUT_SUFFIXES = ('.fq.gz', '.fastq.gz', '.fq', '.fastq')

    def _candidate_outputs(self, base_name: str, suffix: str) -> List[str]:
        """生成候选输出名(精确名 + 锚定 _pass 的 glob)|Generate candidate output names (exact + anchored _pass glob)

        只用精确文件名和紧跟样本名的 _pass 通配, 避免 SRR12 误配 SRR123 的前缀碰撞
        |Exact names plus _pass wildcard anchored right after the sample name,
        so SRR12 never collides with SRR123
        """
        candidates = [
            os.path.join(self.config.output_dir, f"{base_name}{suffix}"),
            os.path.join(self.config.output_dir, f"{base_name}_1{suffix}"),
            os.path.join(self.config.output_dir, f"{base_name}_2{suffix}"),
        ]
        return candidates

    def _is_already_converted(self, sra_file: str) -> bool:
        """
        检查SRA文件是否已转换（断点续传）|Check if SRA file already converted (checkpoint resume)

        通过检查输出目录中是否存在匹配的fastq文件来判断|Judge by checking for matching fastq files in output dir
        """
        base_name = Path(sra_file).stem

        for suffix in self.OUTPUT_SUFFIXES:
            # 单端: SRRxxx.fq.gz / 双端: SRRxxx_1.fq.gz, SRRxxx_2.fq.gz|Single/paired exact names
            if any(os.path.exists(p) for p in self._candidate_outputs(base_name, suffix)):
                return True
            # 过滤技术读段时的 _pass 变体: SRRxxx_pass_1.fq.gz|_pass variants (skip-technical)
            pass_pattern = os.path.join(
                self.config.output_dir, f"{base_name}_pass*{suffix}")
            if glob_module.glob(pass_pattern):
                return True

        return False

    def _rename_fastq_to_fq(self, sra_file: str) -> int:
        """
        将 .fastq(.gz) 统一重命名为 .fq(.gz)|Rename .fastq(.gz) to unified .fq(.gz)

        Returns:
            int: 重命名的文件数|Number of renamed files
        """
        base_name = Path(sra_file).stem
        output_dir = self.config.output_dir
        renamed_count = 0

        for fastq_suffix, fq_suffix in [('.fastq.gz', '.fq.gz'), ('.fastq', '.fq')]:
            # 精确名单 + 锚定 _pass 的 glob(避免前缀误伤)|Exact list + anchored _pass glob (no prefix collisions)
            sources = [p for p in self._candidate_outputs(base_name, fastq_suffix)
                       if os.path.exists(p)]
            sources.extend(glob_module.glob(
                os.path.join(output_dir, f"{base_name}_pass*{fastq_suffix}")))

            for src in sources:
                dst = src[:-len(fastq_suffix)] + fq_suffix
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

        # 根据工具类型构建参数|Build arguments based on tool type
        if self.config.use_parallel:
            args = self._build_parallel_fastq_dump_args(sra_file)
            tool_name = "parallel-fastq-dump (多线程加速)|(multi-threaded)"
        else:
            args = self._build_fastq_dump_args(sra_file)
            tool_name = "fastq-dump (单线程)|(single-threaded)"

        self.logger.info(f"使用工具|Using tool: {tool_name}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")

        # 构建完整命令(功能域环境自动包装)|Build full command (auto domain env wrapping)
        cmd = build_conda_command(self.config.tool_path, args)

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
