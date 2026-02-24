"""
比对统计分析模块 |Alignment Statistics Analysis Module
"""

import json
from pathlib import Path
from typing import Dict, List, Any
from .utils import CommandRunner

class AlignmentStatsAnalyzer:
    """比对统计分析器 |Alignment Statistics Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def get_basic_alignment_stats(self, bam_file: str) -> Dict[str, Any]:
        """获取基本比对统计 |Get basic alignment statistics"""
        self.logger.info(f" 获取比对统计|Getting alignment statistics: {bam_file}")
        
        # samtools flagstat (支持多线程)|samtools flagstat (with multi-threading)
        cmd = f"{self.config.samtools_path} flagstat {bam_file}"
        success, output = self.cmd_runner.run(cmd, f"flagstat for {Path(bam_file).name}", use_threads=True)
        
        if not success:
            return {}
        
        stats = self._parse_flagstat_output(output)
        
        # 添加MAPQ统计|Add MAPQ statistics
        mapq_stats = self._get_mapq_distribution(bam_file)
        stats.update(mapq_stats)
        
        return stats
    
    def _parse_flagstat_output(self, output: str) -> Dict[str, Any]:
        """解析flagstat输出|Parse flagstat output"""
        stats = {}
        lines = output.strip().split('\n')
        
        for line in lines:
            if 'in total' in line:
                stats['total_reads'] = int(line.split()[0])
            elif 'secondary' in line:
                stats['secondary_reads'] = int(line.split()[0])
            elif 'supplementary' in line:
                stats['supplementary_reads'] = int(line.split()[0])
            elif 'duplicates' in line:
                stats['duplicate_reads'] = int(line.split()[0])
            elif 'mapped (' in line and 'primary' not in line:
                stats['mapped_reads'] = int(line.split()[0])
                # 提取比对率|Extract mapping rate
                percentage_part = line.split('(')[1].split('%')[0]
                stats['mapping_rate'] = float(percentage_part)
            elif 'paired in sequencing' in line:
                stats['paired_reads'] = int(line.split()[0])
            elif 'properly paired' in line:
                stats['properly_paired'] = int(line.split()[0])
                # 提取正确配对率|Extract proper pair rate
                percentage_part = line.split('(')[1].split('%')[0]
                stats['proper_pair_rate'] = float(percentage_part)
            elif 'with itself and mate mapped' in line:
                stats['both_mapped'] = int(line.split()[0])
            elif 'singletons' in line:
                stats['singletons'] = int(line.split()[0])
        
        # 计算未比对读段数|Calculate unmapped reads
        if 'total_reads' in stats and 'mapped_reads' in stats:
            stats['unmapped_reads'] = stats['total_reads'] - stats['mapped_reads']
        
        return stats
    
    def _get_mapq_distribution(self, bam_file: str) -> Dict[str, Any]:
        """获取MAPQ分布 |Get MAPQ distribution"""
        cmd = f"{self.config.samtools_path} view -F 4 {bam_file}|cut -f5|sort -n|uniq -c"
        success, output = self.cmd_runner.run(cmd, "MAPQ distribution", use_threads=True)
        
        if not success:
            return {}
        
        mapq_dist = {}
        unique_mapped = 0
        multi_mapped = 0
        
        for line in output.strip().split('\n'):
            if line.strip():
                parts = line.strip().split()
                if len(parts) == 2:
                    count, mapq = int(parts[0]), int(parts[1])
                    mapq_dist[mapq] = count
                    
                    if mapq >= self.config.min_mapq:
                        unique_mapped += count
                    elif mapq > 0:
                        multi_mapped += count
        
        return {
            'mapq_distribution': mapq_dist,
            'unique_mapped': unique_mapped,
            'multi_mapped': multi_mapped
        }
