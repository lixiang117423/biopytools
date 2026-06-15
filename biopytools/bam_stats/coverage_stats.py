"""
覆盖度统计分析模块|Coverage Statistics Analysis Module
"""

import statistics
from pathlib import Path
from typing import Dict, List, Any
from .utils import CommandRunner


class CoverageStatsAnalyzer:
    """覆盖度统计分析器|Coverage Statistics Analyzer"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def get_coverage_stats(self, bam_file: str) -> Dict[str, Any]:
        """获取覆盖度统计（滑窗模式）|Get coverage statistics (sliding window)"""
        self.logger.info(
            f"计算覆盖度统计（{self.config.window_size // 1000}kb窗口，"
            f"{self.config.step_size // 1000}kb步长）|"
            f"Calculating coverage statistics "
            f"({self.config.window_size // 1000}kb window, "
            f"{self.config.step_size // 1000}kb step)"
        )

        chrom_sizes = self._get_chromosome_sizes(bam_file)
        if not chrom_sizes:
            return {}

        return self._calculate_windowed_coverage(bam_file, chrom_sizes)

    def _get_chromosome_sizes(self, bam_file: str) -> Dict[str, int]:
        """获取染色体大小信息|Get chromosome size information"""
        cmd = (
            f"{self.config.samtools_path} view -H {bam_file}"
            f"| grep '^@SQ' | awk '{{print $2,$3}}'"
            f"| sed 's/SN://g' | sed 's/LN://g'"
        )
        success, output = self.cmd_runner.run(
            cmd, "Get chromosome sizes"
        )

        if not success:
            return {}

        chrom_sizes = {}
        for line in output.strip().split('\n'):
            if line.strip():
                parts = line.strip().split()
                if len(parts) == 2:
                    chrom_sizes[parts[0]] = int(parts[1])

        return chrom_sizes

    def _calculate_windowed_coverage(
        self, bam_file: str, chrom_sizes: Dict[str, int]
    ) -> Dict[str, Any]:
        """计算滑窗覆盖度|Calculate windowed coverage"""
        import tempfile

        windows_file = Path(self.config.output_dir) / f"_{Path(bam_file).stem}_windows.bed"
        coverage_file = Path(self.config.output_dir) / f"_{Path(bam_file).stem}_coverage.bed"

        try:
            with open(windows_file, 'w') as f:
                for chrom, size in chrom_sizes.items():
                    for start in range(0, size, self.config.step_size):
                        end = min(start + self.config.window_size, size)
                        f.write(f"{chrom}\t{start}\t{end}\n")

            cmd = (
                f"{self.config.bedtools_path} coverage"
                f" -a {windows_file} -b {bam_file} -mean"
                f" > {coverage_file}"
            )
            success, _ = self.cmd_runner.run(
                cmd, "Calculate windowed coverage", use_threads=False
            )

            if not success:
                return {}

            return self._parse_windowed_coverage_file(coverage_file)

        finally:
            for f in [windows_file, coverage_file]:
                if f.exists():
                    f.unlink()

    def _parse_windowed_coverage_file(
        self, coverage_file: Path
    ) -> Dict[str, Any]:
        """解析滑窗覆盖度文件|Parse windowed coverage file"""
        coverages = []
        chromosome_coverage = {}
        total_windows = 0
        covered_windows = 0

        with open(coverage_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    chrom = parts[0]
                    coverage = float(parts[3])
                    coverages.append(coverage)
                    total_windows += 1
                    if coverage > 0:
                        covered_windows += 1

                    if chrom not in chromosome_coverage:
                        chromosome_coverage[chrom] = []
                    chromosome_coverage[chrom].append(coverage)

        if not coverages:
            return {}

        stats = {
            'total_windows': total_windows,
            'covered_windows': covered_windows,
            'coverage_rate': covered_windows / total_windows if total_windows > 0 else 0,
            'mean_coverage': round(statistics.mean(coverages), 4),
            'median_coverage': round(statistics.median(coverages), 4),
            'max_coverage': round(max(coverages), 4),
            'min_coverage': round(min(coverages), 4),
        }

        coverage_bins = self._calculate_coverage_bins(coverages)
        stats['coverage_distribution'] = coverage_bins

        if len(coverages) > 1:
            stats['coverage_std'] = round(statistics.stdev(coverages), 4)
            stats['coverage_cv'] = round(
                stats['coverage_std'] / stats['mean_coverage']
                if stats['mean_coverage'] > 0 else 0, 4
            )

        chrom_stats = {}
        for chrom, chrom_cov in chromosome_coverage.items():
            if chrom_cov:
                chrom_stats[chrom] = {
                    'windows': len(chrom_cov),
                    'mean_coverage': round(statistics.mean(chrom_cov), 4),
                    'covered_windows': sum(1 for c in chrom_cov if c > 0),
                    'coverage_rate': round(
                        sum(1 for c in chrom_cov if c > 0) / len(chrom_cov), 4
                    ),
                }

        stats['chromosome_coverage'] = chrom_stats
        return stats

    def _calculate_coverage_bins(
        self, coverages: List[float]
    ) -> Dict[str, int]:
        """计算覆盖度分布区间|Calculate coverage distribution bins"""
        if not coverages:
            return {}

        max_coverage = max(coverages)
        bin_size = max(1.0, max_coverage / self.config.coverage_bins)

        bins = {}
        for coverage in coverages:
            bin_index = int(coverage // bin_size)
            bin_key = f"{bin_index * bin_size:.1f}-{(bin_index + 1) * bin_size:.1f}X"
            bins[bin_key] = bins.get(bin_key, 0) + 1

        return bins
