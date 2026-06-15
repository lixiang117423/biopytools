"""质控统计处理器|Quality Control Statistics Processor"""

import os
import json
from ..utils import CommandRunner, format_number


class QCProcessor:
    """质控统计处理器|Quality Control Statistics Processor"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run(self) -> bool:
        """运行质控统计|Run QC statistics"""
        self.logger.info("开始质控统计|Starting QC statistics")

        qc_results = {}

        # 1. k-mer计数统计|K-mer count statistics
        qc_results['kmer_counts'] = self._analyze_kmer_counts()

        # 2. 矩阵统计|Matrix statistics
        qc_results['matrices'] = self._analyze_matrices()

        # 3. 过滤统计|Filter statistics
        qc_results['filtered'] = self._analyze_filtered()

        # 4. 关联统计|Association statistics
        qc_results['association'] = self._analyze_association()

        # 保存QC报告|Save QC report
        report_file = os.path.join(self.config.dirs['qc'], 'qc_report.json')
        with open(report_file, 'w') as f:
            json.dump(qc_results, f, indent=2)

        self.logger.info(f"QC报告已保存|QC report saved: {report_file}")

        return True

    def _analyze_kmer_counts(self):
        """分析k-mer计数|Analyze k-mer counts"""
        counts_dir = self.config.dirs['counts']
        if not os.path.exists(counts_dir):
            return {}

        from pathlib import Path
        kmer_files = list(Path(counts_dir).glob('*.bin'))

        stats = {
            'total_files': len(kmer_files),
            'total_size_bytes': sum(f.stat().st_size for f in kmer_files),
            'avg_size_bytes': 0
        }

        if kmer_files:
            stats['avg_size_bytes'] = stats['total_size_bytes'] / len(kmer_files)

        self.logger.info(f"k-mer文件数|K-mer file count: {stats['total_files']}")
        self.logger.info(f"总大小|Total size: {format_number(stats['total_size_bytes'])} bytes")

        return stats

    def _analyze_matrices(self):
        """分析矩阵|Analyze matrices"""
        matrices_dir = self.config.dirs['matrices']
        if not os.path.exists(matrices_dir):
            return {}

        from pathlib import Path
        matrix_files = list(Path(matrices_dir).glob('*'))

        return {
            'total_files': len(matrix_files),
            'files': [f.name for f in matrix_files[:10]]
        }

    def _analyze_filtered(self):
        """分析过滤结果|Analyze filtered results"""
        filtered_dir = self.config.dirs['filtered']
        if not os.path.exists(filtered_dir):
            return {}

        from pathlib import Path
        filtered_files = list(Path(filtered_dir).glob('*'))

        return {
            'total_files': len(filtered_files)
        }

    def _analyze_association(self):
        """分析关联结果|Analyze association results"""
        asso_dir = self.config.dirs['association']
        if not os.path.exists(asso_dir):
            return {}

        from pathlib import Path
        asso_files = list(Path(asso_dir).glob('*'))

        return {
            'total_files': len(asso_files),
            'files': [f.name for f in asso_files]
        }
