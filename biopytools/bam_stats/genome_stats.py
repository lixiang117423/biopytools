"""
基因组级别统计模块|Genome-level Statistics Module
"""

import json
import statistics
from pathlib import Path
from typing import Dict, List, Any
from collections import defaultdict
from datetime import datetime


class GenomeStatsGenerator:
    """基因组统计生成器|Genome Statistics Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_genome_stats(
        self, sample_stats: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """生成基因组级别统计|Generate genome-level statistics"""
        self.logger.info(
            "生成基因组级别统计|Generating genome-level statistics"
        )

        if not sample_stats:
            return {}

        return {
            'analysis_info': {
                'total_samples': len(sample_stats),
                'analysis_timestamp': datetime.now().isoformat(),
                'samples_analyzed': [
                    s.get('sample_name', 'unknown') for s in sample_stats
                ],
            },
            'alignment_summary': self._aggregate_alignment_stats(sample_stats),
            'coverage_summary': self._aggregate_coverage_stats(sample_stats),
            'sequence_summary': self._aggregate_sequence_stats(sample_stats),
            'insert_size_summary': self._aggregate_insert_stats(sample_stats),
            'duplicate_summary': self._aggregate_duplicate_stats(sample_stats),
            'variation_summary': self._aggregate_variation_stats(sample_stats),
            'chromosome_summary': self._aggregate_chromosome_stats(sample_stats),
        }

    def _aggregate_alignment_stats(
        self, sample_stats: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """汇总比对统计|Aggregate alignment statistics"""
        total_reads = []
        mapping_rates = []
        proper_pair_rates = []

        for stats in sample_stats:
            alignment = stats.get('alignment_stats', {})
            if alignment:
                if 'total_reads' in alignment:
                    total_reads.append(alignment['total_reads'])
                if 'mapping_rate' in alignment:
                    mapping_rates.append(alignment['mapping_rate'])
                if 'proper_pair_rate' in alignment:
                    proper_pair_rates.append(alignment['proper_pair_rate'])

        summary = {}
        if total_reads:
            summary['average_reads_per_sample'] = round(
                statistics.mean(total_reads), 0
            )
            summary['median_reads_per_sample'] = round(
                statistics.median(total_reads), 0
            )

        if mapping_rates:
            summary['average_mapping_rate'] = round(
                statistics.mean(mapping_rates), 2
            )
            summary['min_mapping_rate'] = round(min(mapping_rates), 2)
            summary['max_mapping_rate'] = round(max(mapping_rates), 2)

        if proper_pair_rates:
            summary['average_proper_pair_rate'] = round(
                statistics.mean(proper_pair_rates), 2
            )

        return summary

    def _aggregate_coverage_stats(
        self, sample_stats: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
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
            summary['average_coverage'] = round(
                statistics.mean(mean_coverages), 2
            )
            summary['min_coverage'] = round(min(mean_coverages), 2)
            summary['max_coverage'] = round(max(mean_coverages), 2)

        if coverage_rates:
            summary['average_coverage_rate'] = round(
                statistics.mean(coverage_rates), 4
            )

        return summary

    def _aggregate_sequence_stats(
        self, sample_stats: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """汇总序列特征统计|Aggregate sequence feature statistics"""
        read_lengths = []
        gc_contents = []

        for stats in sample_stats:
            sequence = stats.get('sequence_stats', {})
            if sequence:
                if 'average_read_length' in sequence:
                    read_lengths.append(sequence['average_read_length'])
                if 'average_gc_content' in sequence:
                    gc_contents.append(sequence['average_gc_content'])

        summary = {}
        if read_lengths:
            summary['average_read_length'] = round(
                statistics.mean(read_lengths), 2
            )

        if gc_contents:
            summary['average_gc_content'] = round(
                statistics.mean(gc_contents), 2
            )
            summary['gc_content_range'] = [
                round(min(gc_contents), 2),
                round(max(gc_contents), 2),
            ]

        return summary

    def _aggregate_insert_stats(
        self, sample_stats: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """汇总插入片段统计|Aggregate insert size statistics"""
        insert_sizes = []

        for stats in sample_stats:
            ins = stats.get('insert_stats', {})
            if ins and ins.get('insert_stats_available', False):
                if 'mean_insert_size' in ins:
                    insert_sizes.append(ins['mean_insert_size'])

        summary = {}
        if insert_sizes:
            summary['average_insert_size'] = round(
                statistics.mean(insert_sizes), 2
            )
            summary['insert_size_range'] = [
                round(min(insert_sizes), 2),
                round(max(insert_sizes), 2),
            ]
            summary['samples_with_insert_data'] = len(insert_sizes)
        else:
            summary['note'] = 'No paired-end data available'

        return summary

    def _aggregate_duplicate_stats(
        self, sample_stats: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """汇总重复序列统计|Aggregate duplicate statistics"""
        duplicate_rates = []

        for stats in sample_stats:
            dup = stats.get('duplicate_stats', {})
            if dup:
                if 'duplicate_rate' in dup:
                    duplicate_rates.append(dup['duplicate_rate'])

        summary = {}
        if duplicate_rates:
            summary['average_duplicate_rate'] = round(
                statistics.mean(duplicate_rates), 2
            )
            summary['duplicate_rate_range'] = [
                round(min(duplicate_rates), 2),
                round(max(duplicate_rates), 2),
            ]

        return summary

    def _aggregate_variation_stats(
        self, sample_stats: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """汇总变异统计|Aggregate variation statistics"""
        mismatch_rates = []
        softclip_rates = []
        gc_biases = []

        for stats in sample_stats:
            var = stats.get('variation_stats', {})
            if var:
                if var.get('mismatch_rate') is not None:
                    mismatch_rates.append(var['mismatch_rate'])
                if var.get('softclip_rate') is not None:
                    softclip_rates.append(var['softclip_rate'])
                if var.get('gc_bias') is not None:
                    gc_biases.append(var['gc_bias'])

        summary = {}
        if mismatch_rates:
            summary['average_mismatch_rate'] = round(
                statistics.mean(mismatch_rates), 4
            )
            summary['max_mismatch_rate'] = round(max(mismatch_rates), 4)

        if softclip_rates:
            summary['average_softclip_rate'] = round(
                statistics.mean(softclip_rates), 4
            )

        if gc_biases:
            summary['average_gc_bias'] = round(
                statistics.mean(gc_biases), 2
            )

        return summary

    def _aggregate_chromosome_stats(
        self, sample_stats: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
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
                    'average_coverage': round(
                        statistics.mean(coverages), 4
                    ),
                    'coverage_range': [
                        round(min(coverages), 4),
                        round(max(coverages), 4),
                    ],
                    'samples_analyzed': len(coverages),
                }

        return summary
