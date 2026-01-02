"""
BLAST Statistics Analysis Module
"""

import os
from pathlib import Path
from datetime import datetime
from typing import Dict


class StatisticsGenerator:
    """Generate BLAST analysis statistics reports"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_statistics_report(self, summary_file: str) -> str:
        """
        Generate statistics report

        Args:
            summary_file: BLAST summary results file path

        Returns:
            Path to generated statistics file
        """
        self.logger.info("Generating BLAST statistics report")

        stats_file = os.path.join(self.config.output, "blast_statistics.txt")

        try:
            # Load summary data
            stats = self._load_summary_stats(summary_file)

            # Write statistics report
            with open(stats_file, 'w', encoding='utf-8') as f:
                self._write_basic_info(f)
                self._write_summary_stats(f, stats)
                self._write_sample_stats(f, stats)
                self._write_evalue_distribution(f, stats)
                self._write_identity_distribution(f, stats)
                self._write_coverage_distribution(f, stats)

            self.logger.info(f"Statistics report generated: {stats_file}")
            return stats_file

        except Exception as e:
            self.logger.error(f"Failed to generate statistics report: {e}")
            return ""

    def _load_summary_stats(self, summary_file: str) -> Dict:
        """Load statistics from summary file"""
        stats = {
            'total_alignments': 0,
            'samples_count': 0,
            'unique_queries': set(),
            'unique_subjects': set(),
            'identities': [],
            'evalues': [],
            'coverages': [],
            'sample_stats': {}
        }

        if not os.path.exists(summary_file) or os.path.getsize(summary_file) == 0:
            return stats

        with open(summary_file, 'r', encoding='utf-8') as f:
            # Skip header
            next(f, None)

            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                # 新格式：16列 - Sample qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen coverage qseq sseq
                if len(parts) < 16:
                    continue

                try:
                    sample_name = parts[0]
                    query_id = parts[1]
                    subject_id = parts[2]
                    identity = float(parts[3])
                    evalue = parts[11]
                    coverage = float(parts[14])

                    stats['total_alignments'] += 1
                    stats['unique_queries'].add(query_id)
                    stats['unique_subjects'].add(subject_id)
                    stats['identities'].append(identity)
                    stats['evalues'].append(evalue)
                    stats['coverages'].append(coverage)

                    # Per-sample statistics
                    if sample_name not in stats['sample_stats']:
                        stats['sample_stats'][sample_name] = {
                            'count': 0,
                            'identities': [],
                            'evalues': [],
                            'coverages': []
                        }

                    stats['sample_stats'][sample_name]['count'] += 1
                    stats['sample_stats'][sample_name]['identities'].append(identity)
                    stats['sample_stats'][sample_name]['evalues'].append(evalue)
                    stats['sample_stats'][sample_name]['coverages'].append(coverage)

                except (ValueError, IndexError) as e:
                    self.logger.debug(f"跳过格式错误的行 | Skipping malformed line: {e}")
                    continue

        stats['unique_queries'] = len(stats['unique_queries'])
        stats['unique_subjects'] = len(stats['unique_subjects'])

        return stats

    def _write_basic_info(self, f):
        """Write basic information"""
        f.write("BLAST Alignment Statistics Report\n")
        f.write("=" * 80 + "\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"BLAST Program: {self.config.blast_type.upper()}\n")
        f.write(f"Target Database: {self.config.reference}\n")
        f.write(f"Threads: {self.config.threads}\n")
        f.write(f"E-value Threshold: {self.config.evalue}\n")
        f.write(f"Min Identity: {self.config.min_identity}%\n")
        f.write(f"Min Coverage: {self.config.min_coverage}%\n")
        f.write("\n")

    def _write_summary_stats(self, f, stats: Dict):
        """Write summary statistics"""
        f.write("Summary Statistics\n")
        f.write("-" * 40 + "\n")

        if stats['total_alignments'] == 0:
            f.write("No alignments found\n\n")
            return

        f.write(f"Total Alignments: {stats['total_alignments']}\n")
        f.write(f"Samples Count: {stats['samples_count']}\n")
        f.write(f"Unique Queries: {stats['unique_queries']}\n")
        f.write(f"Unique Subjects: {stats['unique_subjects']}\n")

        if stats['identities']:
            avg_identity = sum(stats['identities']) / len(stats['identities'])
            f.write(f"Average Identity: {avg_identity:.2f}%\n")

        if stats['coverages']:
            avg_coverage = sum(stats['coverages']) / len(stats['coverages'])
            f.write(f"Average Coverage: {avg_coverage:.2f}%\n")

        f.write("\n")

    def _write_sample_stats(self, f, stats: Dict):
        """Write per-sample statistics"""
        f.write("Sample Statistics\n")
        f.write("-" * 40 + "\n")

        if not stats['sample_stats']:
            f.write("No sample data\n\n")
            return

        f.write("Sample\tAlignments\tAvg Identity\tAvg Coverage\n")
        for sample_name, sample_stats in sorted(stats['sample_stats'].items()):
            count = sample_stats['count']
            avg_id = sum(sample_stats['identities']) / len(sample_stats['identities']) if sample_stats['identities'] else 0
            avg_cov = sum(sample_stats['coverages']) / len(sample_stats['coverages']) if sample_stats['coverages'] else 0
            f.write(f"{sample_name}\t{count}\t{avg_id:.2f}%\t{avg_cov:.2f}%\n")

        f.write("\n")

    def _write_evalue_distribution(self, f, stats: Dict):
        """Write E-value distribution"""
        f.write("E-value Distribution\n")
        f.write("-" * 40 + "\n")

        if not stats['evalues']:
            f.write("No E-value data\n\n")
            return

        def parse_evalue(e):
            """Parse e-value string to float"""
            try:
                if e.startswith('e') or e.startswith('E'):
                    return float(f"1{e}")
                return float(e)
            except:
                return 1.0

        evalue_ranges = [
            ("<=1e-50", lambda x: x <= 1e-50),
            ("1e-50 to 1e-30", lambda x: 1e-50 < x <= 1e-30),
            ("1e-30 to 1e-20", lambda x: 1e-30 < x <= 1e-20),
            ("1e-20 to 1e-10", lambda x: 1e-20 < x <= 1e-10),
            ("1e-10 to 1e-5", lambda x: 1e-10 < x <= 1e-5),
            (">1e-5", lambda x: x > 1e-5)
        ]

        parsed_evalues = [parse_evalue(e) for e in stats['evalues']]

        f.write("E-value Range\tCount\n")
        for range_name, condition in evalue_ranges:
            count = sum(1 for e in parsed_evalues if condition(e))
            f.write(f"{range_name}\t{count}\n")

        f.write("\n")

    def _write_identity_distribution(self, f, stats: Dict):
        """Write identity distribution"""
        f.write("Identity Distribution\n")
        f.write("-" * 40 + "\n")

        if not stats['identities']:
            f.write("No identity data\n\n")
            return

        identity_ranges = [
            (">=95%", lambda x: x >= 95),
            ("90-95%", lambda x: 90 <= x < 95),
            ("80-90%", lambda x: 80 <= x < 90),
            ("70-80%", lambda x: 70 <= x < 80),
            ("<70%", lambda x: x < 70)
        ]

        f.write("Identity Range\tCount\n")
        for range_name, condition in identity_ranges:
            count = sum(1 for i in stats['identities'] if condition(i))
            f.write(f"{range_name}\t{count}\n")

        f.write("\n")

    def _write_coverage_distribution(self, f, stats: Dict):
        """Write coverage distribution"""
        f.write("Coverage Distribution\n")
        f.write("-" * 40 + "\n")

        if not stats['coverages']:
            f.write("No coverage data\n\n")
            return

        coverage_ranges = [
            (">=90%", lambda x: x >= 90),
            ("80-90%", lambda x: 80 <= x < 90),
            ("70-80%", lambda x: 70 <= x < 80),
            ("50-70%", lambda x: 50 <= x < 70),
            ("<50%", lambda x: x < 50)
        ]

        f.write("Coverage Range\tCount\n")
        for range_name, condition in coverage_ranges:
            count = sum(1 for c in stats['coverages'] if condition(c))
            f.write(f"{range_name}\t{count}\n")

        f.write("\n")
