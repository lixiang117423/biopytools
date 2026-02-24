"""
覆盖度统计分析模块 |Coverage Statistics Analysis Module
"""

import statistics
from pathlib import Path
from typing import Dict, List, Any
from .utils import CommandRunner

class CoverageStatsAnalyzer:
    """覆盖度统计分析器 |Coverage Statistics Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def get_coverage_stats(self, bam_file: str) -> Dict[str, Any]:
        """获取覆盖度统计（滑窗模式）|Get coverage statistics (sliding window mode)"""
        self.logger.info(f" 计算覆盖度统计（滑窗: {self.config.window_size//1000}kb窗口，{self.config.step_size//1000}kb步长）| Calculating coverage statistics (sliding window: {self.config.window_size//1000}kb window, {self.config.step_size//1000}kb step)")
        
        # 使用bedtools makewindows创建滑窗，然后用bedtools coverage计算覆盖度
        # Use bedtools makewindows to create sliding windows, then bedtools coverage to calculate coverage
        
        # 首先获取BAM文件的染色体信息|First get chromosome information from BAM file
        chrom_sizes = self._get_chromosome_sizes(bam_file)
        if not chrom_sizes:
            return {}
        
        # 创建滑窗|Create sliding windows
        windows_stats = self._calculate_windowed_coverage(bam_file, chrom_sizes)
        
        return windows_stats
    
    def _get_chromosome_sizes(self, bam_file: str) -> Dict[str, int]:
        """获取染色体大小信息 |Get chromosome size information"""
        cmd = f"{self.config.samtools_path} view -H {bam_file}|grep '^@SQ'|awk '{{print $2,$3}}'|sed 's/SN://g'|sed 's/LN://g'"
        success, output = self.cmd_runner.run(cmd, "Get chromosome sizes", use_threads=True)
        
        if not success:
            return {}
        
        chrom_sizes = {}
        for line in output.strip().split('\n'):
            if line.strip():
                parts = line.strip().split()
                if len(parts) == 2:
                    chrom, size = parts[0], int(parts[1])
                    chrom_sizes[chrom] = size
        
        return chrom_sizes
    
    def _calculate_windowed_coverage(self, bam_file: str, chrom_sizes: Dict[str, int]) -> Dict[str, Any]:
        """计算滑窗覆盖度 |Calculate windowed coverage"""
        temp_windows_file = Path(self.config.output_dir) / f"{Path(bam_file).stem}_windows.bed"
        temp_coverage_file = Path(self.config.output_dir) / f"{Path(bam_file).stem}_coverage.bed"
        
        try:
            # 创建滑窗BED文件|Create sliding windows BED file
            with open(temp_windows_file, 'w') as f:
                for chrom, size in chrom_sizes.items():
                    for start in range(0, size, self.config.step_size):
                        end = min(start + self.config.window_size, size)
                        f.write(f"{chrom}\t{start}\t{end}\n")
            
            # 使用bedtools coverage计算每个窗口的覆盖度|Use bedtools coverage to calculate coverage for each window
            cmd = f"{self.config.bedtools_path} coverage -a {temp_windows_file} -b {bam_file} -mean > {temp_coverage_file}"
            success, _ = self.cmd_runner.run(cmd, "Calculate windowed coverage", use_threads=False)
            
            if not success:
                return {}
            
            # 解析覆盖度结果|Parse coverage results
            stats = self._parse_windowed_coverage_file(temp_coverage_file, chrom_sizes)
            
            return stats
            
        finally:
            # 清理临时文件|Clean up temporary files
            for temp_file in [temp_windows_file, temp_coverage_file]:
                if temp_file.exists():
                    temp_file.unlink()
    
    def _parse_windowed_coverage_file(self, coverage_file: Path, chrom_sizes: Dict[str, int]) -> Dict[str, Any]:
        """解析滑窗覆盖度文件 |Parse windowed coverage file"""
        coverages = []
        chromosome_coverage = {}
        total_windows = 0
        covered_windows = 0
        
        with open(coverage_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    chrom, start, end, coverage = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
                    coverages.append(coverage)
                    total_windows += 1
                    
                    if coverage > 0:
                        covered_windows += 1
                    
                    # 按染色体统计|Statistics by chromosome
                    if chrom not in chromosome_coverage:
                        chromosome_coverage[chrom] = []
                    chromosome_coverage[chrom].append(coverage)
        
        if not coverages:
            return {}
        
        # 计算统计指标|Calculate statistics
        stats = {
            'total_windows': total_windows,
            'covered_windows': covered_windows,
            'coverage_rate': covered_windows / total_windows if total_windows > 0 else 0,
            'mean_coverage': statistics.mean(coverages),
            'median_coverage': statistics.median(coverages),
            'max_coverage': max(coverages),
            'min_coverage': min(coverages),
            'window_size': self.config.window_size,
            'step_size': self.config.step_size,
        }
        
        # 计算覆盖度分布|Calculate coverage distribution
        coverage_bins = self._calculate_coverage_bins_windowed(coverages)
        stats['coverage_distribution'] = coverage_bins
        
        # 计算覆盖均匀性|Calculate coverage uniformity
        if len(coverages) > 1:
            stats['coverage_std'] = statistics.stdev(coverages)
            stats['coverage_cv'] = stats['coverage_std'] / stats['mean_coverage'] if stats['mean_coverage'] > 0 else 0
        
        # 染色体特异性统计|Chromosome-specific statistics
        chrom_stats = {}
        for chrom, chrom_coverages in chromosome_coverage.items():
            if chrom_coverages:
                chrom_stats[chrom] = {
                    'windows': len(chrom_coverages),
                    'mean_coverage': statistics.mean(chrom_coverages),
                    'covered_windows': sum(1 for c in chrom_coverages if c > 0),
                    'coverage_rate': sum(1 for c in chrom_coverages if c > 0) / len(chrom_coverages)
                }
        
        stats['chromosome_coverage'] = chrom_stats
        
        return stats
    
    def _calculate_coverage_bins_windowed(self, coverages: List[float]) -> Dict[str, int]:
        """计算滑窗覆盖度分布区间 |Calculate windowed coverage distribution bins"""
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
