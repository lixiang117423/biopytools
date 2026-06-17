"""
比对统计与覆盖度分析模块|Alignment Statistics and Coverage Analysis Module

包含两类后处理：
- AlignmentStatsGenerator: samtools flagstat + stats + 汇总报告
- CoverageAnalyzer: samtools depth + 滑窗覆盖度
|Includes two post-processors:
- AlignmentStatsGenerator: samtools flagstat + stats + summary
- CoverageAnalyzer: samtools depth + windowed coverage
"""

from pathlib import Path

from .utils import build_conda_command, format_number


class AlignmentStatsGenerator:
    """比对统计生成器|Alignment Statistics Generator"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def generate_stats(self, sample_name: str, bam_file: str) -> bool:
        """生成单个样品的flagstat和stats|Generate flagstat and stats for a sample"""
        self.logger.info(f"生成比对统计|Generating alignment stats: {sample_name}")

        sample_dir = Path(bam_file).parent

        # flagstat
        flagstat_file = sample_dir / f"{sample_name}.minibwa.flagstat.txt"
        if not self._run_flagstat(bam_file, flagstat_file):
            return False

        # stats
        stats_file = sample_dir / f"{sample_name}.minibwa.stats.txt"
        if not self._run_stats(bam_file, stats_file):
            return False

        self.logger.info(f"统计生成完成|Stats generated: {sample_name}")
        return True

    def _run_flagstat(self, bam_file: str, output_file: Path) -> bool:
        """运行samtools flagstat|Run samtools flagstat"""
        args = ['flagstat', str(bam_file)]
        cmd = build_conda_command(self.config.samtools_path, args)

        # flagstat结果写入文件需要shell重定向|Redirect to file via shell
        pipeline = f"{' '.join(cmd)} > {output_file}"
        return self.cmd_runner.run_pipeline(
            pipeline,
            f"samtools flagstat|samtools flagstat: {Path(bam_file).name}"
        )

    def _run_stats(self, bam_file: str, output_file: Path) -> bool:
        """运行samtools stats|Run samtools stats"""
        args = ['stats', '-@', str(self.config.threads), str(bam_file)]
        cmd = build_conda_command(self.config.samtools_path, args)

        pipeline = f"{' '.join(cmd)} > {output_file}"
        return self.cmd_runner.run_pipeline(
            pipeline,
            f"samtools stats|samtools stats: {Path(bam_file).name}"
        )

    def generate_summary_report(self, samples: list) -> bool:
        """生成多样本汇总报告|Generate multi-sample summary report"""
        self.logger.info("生成汇总报告|Generating summary report")

        report_file = self.config.pipeline_info_dir / "alignment_summary.txt"

        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("Minibwa比对分析汇总报告|Minibwa Alignment Summary Report\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"参考基因组|Reference genome: {self.config.genome}\n")
            f.write(f"输入目录|Input directory: {self.config.input_dir}\n")
            f.write(f"匹配模式|Pattern: {self.config.pattern}\n")
            f.write(f"比对模式|Alignment mode: {self.config.mode}\n")
            f.write(f"样品数量|Sample count: {len(samples)}\n\n")

            f.write("-" * 80 + "\n")
            f.write("关键参数|Key Parameters\n")
            f.write("-" * 80 + "\n")
            f.write(f"线程数|Threads: {self.config.threads}\n")
            f.write(f"预设|Preset: {self.config.preset}\n")
            f.write(f"最小种子|Min seed (-k): {self.config.min_seed}\n")
            f.write(f"标记重复|Mark duplicates: {'Yes' if self.config.markdup else 'No'}\n")
            f.write(f"移除重复|Remove duplicates: {'Yes' if self.config.remove_dup else 'No'}\n\n")

            f.write("-" * 80 + "\n")
            f.write("已处理样品|Processed Samples\n")
            f.write("-" * 80 + "\n")
            for i, sample in enumerate(samples, 1):
                f.write(f"  {i}. {sample}\n")
            f.write("\n" + "=" * 80 + "\n")

        self.logger.info(f"汇总报告已生成|Summary report: {report_file}")
        return True


class CoverageAnalyzer:
    """覆盖度分析器|Coverage Analyzer"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def analyze_sample_coverage(self, sample_name: str, bam_file: str) -> bool:
        """分析样品覆盖度|Analyze sample coverage"""
        self.logger.info(f"开始覆盖度分析|Starting coverage analysis: {sample_name}")

        sample_dir = Path(bam_file).parent

        # 位点覆盖度|Per-base depth
        depth_file = sample_dir / f"{sample_name}.minibwa.depth.txt"
        if not self._calculate_depth(bam_file, depth_file):
            return False

        # 滑窗覆盖度|Windowed coverage
        window_file = sample_dir / f"{sample_name}.minibwa.window.bed"
        if not self._calculate_windowed_coverage(depth_file, window_file):
            return False

        self.logger.info(f"覆盖度分析完成|Coverage analysis completed: {sample_name}")
        return True

    def _calculate_depth(self, bam_file: str, output_file: Path) -> bool:
        """samtools depth|Run samtools depth"""
        args = ['depth']
        if self.config.min_base_quality > 0:
            args.extend(['-q', str(self.config.min_base_quality)])
        if self.config.min_mapping_quality > 0:
            args.extend(['-Q', str(self.config.min_mapping_quality)])
        if self.config.max_depth > 0:
            args.extend(['-d', str(self.config.max_depth)])
        args.append(str(bam_file))

        cmd = build_conda_command(self.config.samtools_path, args)
        pipeline = f"{' '.join(cmd)} > {output_file}"

        return self.cmd_runner.run_pipeline(
            pipeline,
            f"samtools depth|samtools depth: {Path(bam_file).name}"
        )

    def _calculate_windowed_coverage(self, depth_file: Path,
                                     output_file: Path) -> bool:
        """计算滑窗覆盖度|Calculate windowed coverage"""
        self.logger.info(
            f"计算滑窗覆盖度|Calculating windowed coverage: "
            f"窗口|window={format_number(self.config.window_size)}bp, "
            f"步长|step={format_number(self.config.step_size)}bp"
        )

        try:
            windows = {}
            with open(depth_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 3:
                        continue
                    chrom = parts[0]
                    pos = int(parts[1])
                    depth = int(parts[2])

                    window_idx = (pos - 1) // self.config.step_size
                    window_start = window_idx * self.config.step_size + 1
                    window_end = window_start + self.config.window_size - 1
                    key = (chrom, window_start, window_end)

                    if key not in windows:
                        windows[key] = {'sum': 0, 'count': 0}
                    windows[key]['sum'] += depth
                    windows[key]['count'] += 1

            with open(output_file, 'w') as f:
                for (chrom, start, end), data in sorted(windows.items()):
                    avg = data['sum'] / data['count'] if data['count'] > 0 else 0
                    f.write(f"{chrom}\t{start}\t{end}\t{avg:.2f}\n")

            self.logger.info(f"滑窗覆盖度计算完成|Windowed coverage done: {output_file.name}")
            return True

        except Exception as e:
            self.logger.error(f"滑窗覆盖度计算失败|Windowed coverage failed: {e}")
            return False
