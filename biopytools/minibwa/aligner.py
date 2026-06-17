"""
序列比对模块|Sequence Alignment Module

封装 `minibwa map`，按规范13.2.1方案A/C构建管道：直接用二进制完整路径，
对conda环境内的samtools设置LD_LIBRARY_PATH避免库加载失败。
|Wraps `minibwa map`. Per spec 13.2.1 plan A/C, uses full binary paths directly
and sets LD_LIBRARY_PATH for conda env samtools.
"""

import os
import re
import shlex
from pathlib import Path

from .utils import build_conda_command


class MinibwaAligner:
    """Minibwa比对处理器|Minibwa Aligner"""

    def __init__(self, config, logger, cmd_runner, index_prefix: str):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.index_prefix = index_prefix

    def align_sample(self, sample_name: str, read1: str, read2: str) -> str:
        """
        比对单个样品|Align single sample

        Args:
            sample_name: 样品名|Sample name
            read1: R1文件|R1 file
            read2: R2文件，None表示单端|R2 file, None for single-end

        Returns:
            最终BAM路径，失败返回None|Final BAM path or None on failure
        """
        self.logger.info(f"开始比对样品|Starting alignment: {sample_name}")

        # 样品输出目录|Per-sample output directory
        sample_dir = self.config.output_path / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)

        # 输出文件命名遵循{sample}.{tool}.{ext}|Naming follows {sample}.{tool}.{ext}
        sorted_bam = sample_dir / f"{sample_name}.minibwa.bam"
        sorted_bai = sample_dir / f"{sample_name}.minibwa.bam.bai"

        # 断点续传|Checkpoint resume
        if self.config.resume and sorted_bam.exists() and sorted_bai.exists():
            self.logger.info(f"跳过已完成样品|Skipping completed sample: {sample_name}")
            return str(sorted_bam)

        # Step 1: 比对+排序（管道）|Step 1: Align + sort (pipeline)
        if not self._align_and_sort(sample_name, read1, read2, sorted_bam):
            return None

        # Step 2: 标记重复（可选）|Step 2: Mark duplicates (optional)
        final_bam = sorted_bam
        if self.config.markdup:
            final_bam = sample_dir / f"{sample_name}.minibwa.markdup.bam"
            if not self._mark_duplicates(sorted_bam, final_bam):
                return None
            # markdup产物替代原bam|markdup output replaces sorted bam
            sorted_bam.unlink(missing_ok=True)
            sorted_bam = final_bam

        # Step 3: 构建BAM索引|Step 3: Build BAM index
        if not self._index_bam(sorted_bam):
            return None

        self.logger.info(f"样品比对完成|Sample aligned: {sample_name}")
        return str(sorted_bam)

    def _align_and_sort(self, sample_name: str, read1: str, read2: str,
                        output_bam: Path) -> bool:
        """管道：minibwa map | samtools sort|Pipeline: minibwa map | samtools sort"""
        # minibwa map 参数|minibwa map args
        map_opts = self.config.get_map_options()

        # 构建minibwa命令段（含必要的环境变量）|Build minibwa cmd segment
        minibwa_seg = self._pipeline_segment(
            self.config.minibwa_path, ['map'] + map_opts + [self.index_prefix, read1]
        )
        if read2:
            # 追加R2到minibwa段（已shlex.join）|Append R2 to minibwa segment
            minibwa_seg = f"{minibwa_seg} {shlex.quote(read2)}"

        # 构建samtools sort命令段|Build samtools sort segment
        sort_opts = ['sort', '-@', str(self.config.threads), '-o', str(output_bam), '-']
        samtools_seg = self._pipeline_segment(self.config.samtools_path, sort_opts)

        pipeline = f"{minibwa_seg} | {samtools_seg}"

        return self.cmd_runner.run_pipeline(
            pipeline,
            f"minibwa比对+排序|minibwa alignment+sort: {sample_name}"
        )

    def _mark_duplicates(self, input_bam: Path, output_bam: Path) -> bool:
        """标记或移除重复|Mark or remove duplicates"""
        stats_dir = self.config.output_path / input_bam.parent.name
        stats_file = stats_dir / f"{input_bam.stem}.markdup_stats.txt"

        args = ['markdup', '-@', str(self.config.threads)]
        if self.config.remove_dup:
            args.append('-r')
        args.extend(['-f', str(stats_file), str(input_bam), str(output_bam)])

        cmd = build_conda_command(self.config.samtools_path, args)
        return self.cmd_runner.run(
            cmd,
            f"标记重复|Marking duplicates: {input_bam.name}"
        )

    def _index_bam(self, bam_file: Path) -> bool:
        """构建BAM索引|Build BAM index"""
        args = ['index', '-@', str(self.config.threads), str(bam_file)]
        cmd = build_conda_command(self.config.samtools_path, args)
        return self.cmd_runner.run(
            cmd,
            f"构建BAM索引|Building BAM index: {bam_file.name}"
        )

    def _pipeline_segment(self, binary_path: str, args: list) -> str:
        """
        构建管道中单个命令段（避免conda run | conda run）|Build pipeline segment

        按规范13.2.1方案A：conda env内的工具，设置LD_LIBRARY_PATH并直接调用完整路径
        |Per spec 13.2.1 plan A: for conda env tools, set LD_LIBRARY_PATH and call directly

        Args:
            binary_path: 二进制完整路径|Full binary path
            args: 参数列表（未转义）|Argument list (unquoted)

        Returns:
            可放入管道的shell命令段（已正确转义）|Shell-safe command segment
        """
        # 检测conda env并准备环境变量前缀|Detect conda env, prep env prefix
        env_prefix = ""
        match = re.search(r'(/envs/[^/]+)', binary_path)
        if match:
            env_root = match.group(1)
            env_lib = f"{env_root}/lib"
            current_ld = os.environ.get('LD_LIBRARY_PATH', '')
            env_prefix = f"export LD_LIBRARY_PATH={shlex.quote(env_lib)}:{shlex.quote(current_ld)}; "

        quoted_args = ' '.join(shlex.quote(str(a)) for a in args)
        return f"{env_prefix}{shlex.quote(binary_path)} {quoted_args}"
