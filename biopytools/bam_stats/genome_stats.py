"""
基因组级别统计模块 🌍|Genome-level Statistics Module
"""

import json
import statistics
from pathlib import Path
from typing import Dict, List, Any
from collections import defaultdict
from datetime import datetime

class GenomeStatsGenerator:
    """基因组统计生成器 🌍|Genome Statistics Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_genome_stats(self, sample_stats: List[Dict[str, Any]]) -> Dict[str, Any]:
        """生成基因组级别统计|Generate genome-level statistics"""
        self.logger.info("🌍 生成基因组级别统计|Generating genome-level statistics")
        
        if not sample_stats:
            return {}
        
        genome_stats = {
            'analysis_info': {
                'total_samples': len(sample_stats),
                'analysis_timestamp': datetime.now().isoformat(),
                'samples_analyzed': [stats.get('sample_name', 'unknown') for stats in sample_stats]
            }
        }
        
        # 比对统计汇总|Alignment statistics summary
        genome_stats['alignment_summary'] = self._aggregate_alignment_stats(sample_stats)
        
        # 覆盖度统计汇总|Coverage statistics summary
        genome_stats['coverage_summary'] = self._aggregate_coverage_stats(sample_stats)
        
        # 序列特征汇总|Sequence feature summary
        genome_stats['sequence_summary'] = self._aggregate_sequence_stats(sample_stats)
        
        # 插入片段汇总|Insert size summary
        genome_stats['insert_size_summary'] = self._aggregate_insert_stats(sample_stats)
        
        # 重复序列汇总|Duplicate summary
        genome_stats['duplicate_summary'] = self._aggregate_duplicate_stats(sample_stats)
        
        # 染色体级别汇总|Chromosome-level summary
        genome_stats['chromosome_summary'] = self._aggregate_chromosome_stats(sample_stats)
        
        return genome_stats
    
    def _aggregate_alignment_stats(self, sample_stats: List[Dict[str, Any]]) -> Dict[str, Any]:
        """汇总比对统计|Aggregate alignment statistics"""
        total_reads = []
        mapping_rates = []
        proper_pair_rates = []
        unique_mapped_counts = []
        
        for stats in sample_stats:
            alignment = stats.get('alignment_stats', {})
            if alignment:
                if 'total_reads' in alignment:
                    total_reads.append(alignment['total_reads'])
                if 'mapping_rate' in alignment:
                    mapping_rates.append(alignment['mapping_rate'])
                if 'proper_pair_rate' in alignment:
                    proper_pair_rates.append(alignment['proper_pair_rate'])
                if 'unique_mapped' in alignment:
                    unique_mapped_counts.append(alignment['unique_mapped'])
        
        summary = {}
        if total_reads:
            summary['total_reads_across_samples'] = sum(total_reads)
            summary['average_reads_per_sample'] = statistics.mean(total_reads)
            summary['median_reads_per_sample'] = statistics.median(total_reads)
        
        if mapping_rates:
            summary['average_mapping_rate'] = statistics.mean(mapping_rates)
            summary['median_mapping_rate'] = statistics.median(mapping_rates)
            summary['min_mapping_rate'] = min(mapping_rates)
            summary['max_mapping_rate'] = max(mapping_rates)
        
        if proper_pair_rates:
            summary['average_proper_pair_rate'] = statistics.mean(proper_pair_rates)
        
        if unique_mapped_counts:
            summary['total_unique_mapped'] = sum(unique_mapped_counts)
        
        return summary
    
    def _aggregate_coverage_stats(self, sample_stats: List[Dict[str, Any]]) -> Dict[str, Any]:
        """汇总覆盖度统计|Aggregate coverage statistics"""
        mean_coverages = []
        coverage_rates = []
        
        for stats in sample_stats:
            coverage = stats.get('coverage_stats', {})
            if coverage:
                if 'mean_coverage' in coverage:
                    mean_coverages.append(coverage['mean_coverage'])
                if 'coverage_rate' in coverage:
                    coverage_rates.append(coverage['coverage_rate'])
        
        summary = {}
        if mean_coverages:
            summary['average_coverage_across_samples'] = statistics.mean(mean_coverages)
            summary['median_coverage_across_samples'] = statistics.median(mean_coverages)
            summary['min_coverage'] = min(mean_coverages)
            summary['max_coverage'] = max(mean_coverages)
        
        if coverage_rates:
            summary['average_coverage_rate'] = statistics.mean(coverage_rates)
        
        return summary
    
    def _aggregate_sequence_stats(self, sample_stats: List[Dict[str, Any]]) -> Dict[str, Any]:
        """汇总序列特征统计|Aggregate sequence feature statistics"""
        read_lengths = []
        gc_contents = []
        base_qualities = []
        
        for stats in sample_stats:
            sequence = stats.get('sequence_stats', {})
            if sequence:
                if 'average_read_length' in sequence:
                    read_lengths.append(sequence['average_read_length'])
                if 'average_gc_content' in sequence:
                    gc_contents.append(sequence['average_gc_content'])
                if 'average_base_quality' in sequence:
                    base_qualities.append(sequence['average_base_quality'])
        
        summary = {}
        if read_lengths:
            summary['average_read_length'] = statistics.mean(read_lengths)
        
        if gc_contents:
            summary['average_gc_content'] = statistics.mean(gc_contents)
            summary['gc_content_range'] = [min(gc_contents), max(gc_contents)]
        
        if base_qualities:
            summary['average_base_quality'] = statistics.mean(base_qualities)
        
        return summary
    
    def _aggregate_insert_stats(self, sample_stats: List[Dict[str, Any]]) -> Dict[str, Any]:
        """汇总插入片段统计|Aggregate insert size statistics"""
        insert_sizes = []
        
        for stats in sample_stats:
            insert_stats = stats.get('insert_stats', {})
            if insert_stats and insert_stats.get('insert_stats_available', False):
                if 'mean_insert_size' in insert_stats:
                    insert_sizes.append(insert_stats['mean_insert_size'])
        
        summary = {}
        if insert_sizes:
            summary['average_insert_size'] = statistics.mean(insert_sizes)
            summary['insert_size_range'] = [min(insert_sizes), max(insert_sizes)]
            summary['samples_with_insert_data'] = len(insert_sizes)
        else:
            summary['note'] = 'No paired-end data available for insert size analysis'
        
        return summary
    
    def _aggregate_duplicate_stats(self, sample_stats: List[Dict[str, Any]]) -> Dict[str, Any]:
        """汇总重复序列统计|Aggregate duplicate statistics"""
        duplicate_rates = []
        
        for stats in sample_stats:
            dup_stats = stats.get('duplicate_stats', {})
            if dup_stats:
                if 'duplicate_rate' in dup_stats:
                    duplicate_rates.append(dup_stats['duplicate_rate'])
                elif 'estimated_duplicate_rate' in dup_stats:
                    duplicate_rates.append(dup_stats['estimated_duplicate_rate'])
        
        summary = {}
        if duplicate_rates:
            summary['average_duplicate_rate'] = statistics.mean(duplicate_rates)
            summary['duplicate_rate_range'] = [min(duplicate_rates), max(duplicate_rates)]
        
        return summary
    
    def _aggregate_chromosome_stats(self, sample_stats: List[Dict[str, Any]]) -> Dict[str, Any]:
        """汇总染色体统计|Aggregate chromosome statistics"""
        chromosome_data = defaultdict(list)
        
        for stats in sample_stats:
            coverage = stats.get('coverage_stats', {})
            chrom_coverage = coverage.get('chromosome_coverage', {})
            
            for chrom, chrom_stats in chrom_coverage.items():
                if 'mean_coverage' in chrom_stats:
                    chromosome_data[chrom].append(chrom_stats['mean_coverage'])
        
        summary = {}
        for chrom, coverages in chromosome_data.items():
            if coverages:
                summary[chrom] = {
                    'average_coverage': statistics.mean(coverages),
                    'coverage_range': [min(coverages), max(coverages)],
                    'samples_analyzed': len(coverages)
                }
        
        return summary
