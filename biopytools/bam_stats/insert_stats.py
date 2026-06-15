"""
插入片段统计分析模块|Insert Size Statistics Analysis Module
"""

import statistics
from pathlib import Path
from typing import Dict, List, Any
from .utils import CommandRunner


class InsertStatsAnalyzer:
    """插入片段统计分析器|Insert Size Statistics Analyzer"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def get_insert_stats(self, bam_file: str) -> Dict[str, Any]:
        """获取插入片段统计|Get insert size statistics"""
        self.logger.info(
            f"分析插入片段大小|Analyzing insert sizes: {Path(bam_file).name}"
        )

        cmd = (
            f"{self.config.samtools_path} view -f 2 -F 4 {bam_file}"
            f"| awk '{{if($9>0 && $9<={self.config.max_insert_size}) print $9}}'"
        )
        success, output = self.cmd_runner.run(
            cmd, "Extract insert sizes"
        )

        if not success:
            return {}

        insert_sizes = []
        for line in output.strip().split('\n'):
            if line.strip():
                try:
                    size = int(line.strip())
                    if 0 < size <= self.config.max_insert_size:
                        insert_sizes.append(size)
                except ValueError:
                    continue

        if not insert_sizes:
            return {'insert_stats_available': False}

        stats = {
            'insert_stats_available': True,
            'total_paired_reads': len(insert_sizes),
            'mean_insert_size': round(statistics.mean(insert_sizes), 2),
            'median_insert_size': round(statistics.median(insert_sizes), 2),
            'min_insert_size': min(insert_sizes),
            'max_insert_size': max(insert_sizes),
        }

        if len(insert_sizes) > 1:
            stats['insert_size_std'] = round(statistics.stdev(insert_sizes), 2)

        stats['insert_size_distribution'] = self._calculate_insert_bins(insert_sizes)
        return stats

    def _calculate_insert_bins(
        self, insert_sizes: List[int]
    ) -> Dict[str, int]:
        """计算插入大小分布区间|Calculate insert size distribution bins"""
        max_size = max(insert_sizes)
        bin_size = max(1, max_size // 50)

        bins = {}
        for size in insert_sizes:
            bin_key = (
                f"{(size // bin_size) * bin_size}-"
                f"{(size // bin_size + 1) * bin_size - 1}"
            )
            bins[bin_key] = bins.get(bin_key, 0) + 1

        return bins
