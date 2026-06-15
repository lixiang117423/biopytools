"""
样品级别统计模块|Sample-level Statistics Module
"""

from typing import Dict, Any
from .alignment_stats import AlignmentStatsAnalyzer
from .coverage_stats import CoverageStatsAnalyzer
from .sequence_stats import SequenceStatsAnalyzer
from .insert_stats import InsertStatsAnalyzer
from .duplicate_stats import DuplicateStatsAnalyzer
from .variation_stats import VariationStatsAnalyzer
from .utils import get_sample_name


class SampleStatsGenerator:
    """样品统计生成器|Sample Statistics Generator"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

        self.alignment_analyzer = AlignmentStatsAnalyzer(
            config, logger, cmd_runner
        )
        self.coverage_analyzer = CoverageStatsAnalyzer(
            config, logger, cmd_runner
        )
        self.sequence_analyzer = SequenceStatsAnalyzer(
            config, logger, cmd_runner
        )
        self.insert_analyzer = InsertStatsAnalyzer(
            config, logger, cmd_runner
        )
        self.duplicate_analyzer = DuplicateStatsAnalyzer(
            config, logger, cmd_runner
        )
        self.variation_analyzer = VariationStatsAnalyzer(
            config, logger, cmd_runner
        )

    def analyze_sample(self, bam_file: str) -> Dict[str, Any]:
        """分析单个样品|Analyze single sample"""
        sample_name = get_sample_name(bam_file)
        self.logger.info(
            f"开始分析样品|Starting sample analysis: {sample_name}"
        )

        sample_stats = {
            'sample_name': sample_name,
            'bam_file': bam_file,
        }

        if not self.config.skip_alignment:
            self.logger.info("分析比对统计|Analyzing alignment statistics")
            sample_stats['alignment_stats'] = (
                self.alignment_analyzer.get_basic_alignment_stats(bam_file)
            )

        if not self.config.skip_coverage:
            self.logger.info("分析覆盖度统计|Analyzing coverage statistics")
            sample_stats['coverage_stats'] = (
                self.coverage_analyzer.get_coverage_stats(bam_file)
            )

        if not self.config.skip_sequence:
            self.logger.info("分析序列特征|Analyzing sequence features")
            sample_stats['sequence_stats'] = (
                self.sequence_analyzer.get_sequence_stats(bam_file)
            )

        if not self.config.skip_insert:
            self.logger.info("分析插入片段|Analyzing insert sizes")
            sample_stats['insert_stats'] = (
                self.insert_analyzer.get_insert_stats(bam_file)
            )

        if not self.config.skip_duplicate:
            self.logger.info("分析重复序列|Analyzing duplicates")
            sample_stats['duplicate_stats'] = (
                self.duplicate_analyzer.get_duplicate_stats(bam_file)
            )

        if not self.config.skip_variation:
            self.logger.info(
                "分析变异与特异性|Analyzing variation & specificity"
            )
            sample_stats['variation_stats'] = (
                self.variation_analyzer.get_variation_stats(bam_file)
            )

        self.logger.info(
            f"样品分析完成|Sample analysis completed: {sample_name}"
        )
        return sample_stats
