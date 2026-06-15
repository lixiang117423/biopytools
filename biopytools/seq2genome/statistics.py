"""
统计报告生成模块|Statistics Report Generation Module
"""

import os
from typing import List, Dict
from collections import defaultdict
import numpy as np
from .parser import PAFRecord


class StatisticsGenerator:
    """统计报告生成器|Statistics Report Generator"""

    def __init__(self, config, logger):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def generate_statistics(self, records: List[PAFRecord], output_file: str = None):
        """生成统计报告|Generate statistics report

        Args:
            records: PAF记录列表|List of PAF records
            output_file: 输出文件路径|Output file path
        """
        if output_file is None:
            output_file = os.path.join(self.config.output_dir, "alignment_statistics.txt")

        self.logger.info(f"生成统计报告|Generating statistics report: {output_file}")

        with open(output_file, 'w') as f:
            # 标题|Title
            f.write("=" * 80 + "\n")
            f.write("蛋白质到基因组比对统计报告|Protein to Genome Alignment Statistics Report\n")
            f.write("=" * 80 + "\n\n")

            # 总体统计|Overall statistics
            self._write_overall_stats(f, records)

            # 查询序列统计|Query statistics
            self._write_query_stats(f, records)

            # 目标序列统计|Target statistics
            self._write_target_stats(f, records)

            # 质量统计|Quality statistics
            self._write_quality_stats(f, records)

            # 覆盖度统计|Coverage statistics
            self._write_coverage_stats(f, records)

        self.logger.info(f"统计报告已保存|Statistics report saved: {output_file}")

    def _write_overall_stats(self, f, records: List[PAFRecord]):
        """写入总体统计|Write overall statistics"""
        f.write("1. 总体统计|Overall Statistics\n")
        f.write("-" * 80 + "\n")
        f.write(f"总比对记录数|Total alignments: {len(records)}\n")
        f.write(f"唯一查询序列数|Unique queries: {len(set(r.query_name for r in records))}\n")
        f.write(f"唯一目标序列数|Unique targets: {len(set(r.target_name for r in records))}\n")
        f.write(f"正向链|Forward strand (+): {sum(1 for r in records if r.strand == '+')}\n")
        f.write(f"反向链|Reverse strand (-): {sum(1 for r in records if r.strand == '-')}\n")
        f.write("\n")

    def _write_query_stats(self, f, records: List[PAFRecord]):
        """写入查询序列统计|Write query statistics"""
        f.write("2. 查询序列统计|Query Statistics\n")
        f.write("-" * 80 + "\n")

        # 按查询序列分组|Group by query
        query_groups = defaultdict(list)
        for record in records:
            query_groups[record.query_name].append(record)

        f.write(f"查询序列总数|Total query sequences: {len(query_groups)}\n")

        # 每个查询序列的比对数|Alignments per query
        alignment_counts = [len(alignments) for alignments in query_groups.values()]
        f.write(f"平均比对数/查询|Avg alignments/query: {np.mean(alignment_counts):.2f}\n")
        f.write(f"最大比对数|Max alignments: {max(alignment_counts)}\n")
        f.write(f"最小比对数|Min alignments: {min(alignment_counts)}\n")

        # 多重比对统计|Multiple alignment statistics
        multi_aligned = sum(1 for count in alignment_counts if count > 1)
        f.write(f"多重比对查询数|Queries with multiple alignments: {multi_aligned} "
                f"({multi_aligned/len(query_groups)*100:.1f}%)\n")
        f.write("\n")

    def _write_target_stats(self, f, records: List[PAFRecord]):
        """写入目标序列统计|Write target statistics"""
        f.write("3. 目标序列统计|Target Statistics\n")
        f.write("-" * 80 + "\n")

        # 按目标序列分组|Group by target
        target_groups = defaultdict(list)
        for record in records:
            target_groups[record.target_name].append(record)

        f.write(f"目标序列总数|Total target sequences: {len(target_groups)}\n")

        # 每个目标序列的比对数|Alignments per target
        alignment_counts = [len(alignments) for alignments in target_groups.values()]
        f.write(f"平均比对数/目标|Avg alignments/target: {np.mean(alignment_counts):.2f}\n")
        f.write(f"最大比对数|Max alignments: {max(alignment_counts)}\n")
        f.write(f"最小比对数|Min alignments: {min(alignment_counts)}\n")
        f.write("\n")

    def _write_quality_stats(self, f, records: List[PAFRecord]):
        """写入质量统计|Write quality statistics"""
        f.write("4. 比对质量统计|Alignment Quality Statistics\n")
        f.write("-" * 80 + "\n")

        # Mapping quality
        mq_values = [r.mapping_quality for r in records]
        f.write(f"Mapping Quality:\n")
        f.write(f"  平均|Average: {np.mean(mq_values):.2f}\n")
        f.write(f"  中位数|Median: {np.median(mq_values):.2f}\n")
        f.write(f"  最小|Min: {min(mq_values)}\n")
        f.write(f"  最大|Max: {max(mq_values)}\n")

        # Identity
        identity_values = [r.identity for r in records]
        f.write(f"\n序列一致性|Sequence Identity (%):\n")
        f.write(f"  平均|Average: {np.mean(identity_values):.2f}\n")
        f.write(f"  中位数|Median: {np.median(identity_values):.2f}\n")
        f.write(f"  最小|Min: {min(identity_values):.2f}\n")
        f.write(f"  最大|Max: {max(identity_values):.2f}\n")

        # High-quality alignments (MQ >= 30, Identity >= 80)
        high_quality = sum(1 for r in records if r.mapping_quality >= 30 and r.identity >= 80)
        f.write(f"\n高质量比对数（MQ≥30且Identity≥80%）|High-quality alignments (MQ≥30 & ID≥80%): "
                f"{high_quality} ({high_quality/len(records)*100:.1f}%)\n")
        f.write("\n")

    def _write_coverage_stats(self, f, records: List[PAFRecord]):
        """写入覆盖度统计|Write coverage statistics"""
        f.write("5. 覆盖度统计|Coverage Statistics\n")
        f.write("-" * 80 + "\n")

        # Query coverage
        query_cov_values = [r.query_coverage for r in records]
        f.write(f"查询序列覆盖度|Query Coverage (%):\n")
        f.write(f"  平均|Average: {np.mean(query_cov_values):.2f}\n")
        f.write(f"  中位数|Median: {np.median(query_cov_values):.2f}\n")
        f.write(f"  最小|Min: {min(query_cov_values):.2f}\n")
        f.write(f"  最大|Max: {max(query_cov_values):.2f}\n")

        # Target coverage
        target_cov_values = [r.target_coverage for r in records]
        f.write(f"\n目标序列覆盖度|Target Coverage (%):\n")
        f.write(f"  平均|Average: {np.mean(target_cov_values):.2f}\n")
        f.write(f"  中位数|Median: {np.median(target_cov_values):.2f}\n")
        f.write(f"  最小|Min: {min(target_cov_values):.2f}\n")
        f.write(f"  最大|Max: {max(target_cov_values):.2f}\n")
        f.write("\n")
