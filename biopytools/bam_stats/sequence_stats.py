"""
序列特征统计分析模块 |Sequence Feature Statistics Analysis Module
"""

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Any
from .utils import CommandRunner

class SequenceStatsAnalyzer:
    """序列特征统计分析器 |Sequence Feature Statistics Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def get_sequence_stats(self, bam_file: str) -> Dict[str, Any]:
        """获取序列特征统计|Get sequence feature statistics"""
        self.logger.info(f" 分析序列特征|Analyzing sequence features: {bam_file}")
        
        # 获取读长分布|Get read length distribution
        read_length_stats = self._get_read_length_distribution(bam_file)
        
        # 获取GC含量分布|Get GC content distribution
        gc_stats = self._get_gc_content_distribution(bam_file)
        
        # 获取碱基质量分布|Get base quality distribution
        quality_stats = self._get_base_quality_distribution(bam_file)
        
        # 合并统计结果|Merge statistics
        stats = {}
        stats.update(read_length_stats)
        stats.update(gc_stats)
        stats.update(quality_stats)
        
        return stats
    
    def _get_read_length_distribution(self, bam_file: str) -> Dict[str, Any]:
        """获取读长分布 |Get read length distribution"""
        cmd = f"{self.config.samtools_path} view -F 4 {bam_file}|awk '{{print length($10)}}'|sort -n|uniq -c"
        success, output = self.cmd_runner.run(cmd, "Read length distribution", use_threads=True)
        
        if not success:
            return {}
        
        length_dist = {}
        total_reads = 0
        length_sum = 0
        
        for line in output.strip().split('\n'):
            if line.strip():
                parts = line.strip().split()
                if len(parts) == 2:
                    count, length = int(parts[0]), int(parts[1])
                    length_dist[length] = count
                    total_reads += count
                    length_sum += length * count
        
        stats = {
            'read_length_distribution': length_dist,
            'total_reads_for_length': total_reads,
            'average_read_length': length_sum / total_reads if total_reads > 0 else 0
        }
        
        if length_dist:
            stats['min_read_length'] = min(length_dist.keys())
            stats['max_read_length'] = max(length_dist.keys())
        
        return stats
    
    def _get_gc_content_distribution(self, bam_file: str) -> Dict[str, Any]:
        """获取GC含量分布 |Get GC content distribution"""
        # 使用samtools提取序列并计算GC含量（支持多线程）| Use samtools to extract sequences and calculate GC content (with multi-threading)
        cmd = f"{self.config.samtools_path} view -F 4 {bam_file}|cut -f10|head -10000"  # 采样分析
        success, output = self.cmd_runner.run(cmd, "GC content analysis", use_threads=True)
        
        if not success:
            return {}
        
        gc_contents = []
        
        for line in output.strip().split('\n'):
            if line.strip():
                seq = line.strip()
                if seq and seq != '*':
                    gc_count = seq.count('G') + seq.count('C')
                    gc_content = (gc_count / len(seq)) * 100 if len(seq) > 0 else 0
                    gc_contents.append(gc_content)
        
        if not gc_contents:
            return {}
        
        # 计算GC含量分布区间|Calculate GC content distribution bins
        gc_bins = defaultdict(int)
        for gc in gc_contents:
            bin_key = f"{int(gc//5)*5}-{int(gc//5)*5+4}%"
            gc_bins[bin_key] += 1
        
        return {
            'gc_content_distribution': dict(gc_bins),
            'average_gc_content': sum(gc_contents) / len(gc_contents),
            'sampled_reads_for_gc': len(gc_contents)
        }
    
    def _get_base_quality_distribution(self, bam_file: str) -> Dict[str, Any]:
        """获取碱基质量分布|Get base quality distribution"""
        cmd = f"{self.config.samtools_path} view -F 4 {bam_file}|cut -f11|head -1000"  # 采样分析
        success, output = self.cmd_runner.run(cmd, "Base quality analysis")
        
        if not success:
            return {}
        
        quality_counts = defaultdict(int)
        total_bases = 0
        
        for line in output.strip().split('\n'):
            if line.strip():
                qual_str = line.strip()
                if qual_str and qual_str != '*':
                    for qual_char in qual_str:
                        qual_score = ord(qual_char) - 33  # Phred+33编码
                        quality_counts[qual_score] += 1
                        total_bases += 1
        
        if not quality_counts:
            return {}
        
        # 计算平均质量|Calculate average quality
        weighted_sum = sum(qual * count for qual, count in quality_counts.items())
        average_quality = weighted_sum / total_bases if total_bases > 0 else 0
        
        return {
            'base_quality_distribution': dict(quality_counts),
            'average_base_quality': average_quality,
            'total_bases_analyzed': total_bases
        }
