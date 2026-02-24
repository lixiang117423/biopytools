"""
样品级别统计模块|Sample-level Statistics Module
"""

import json
from pathlib import Path
from typing import Dict, Any
from .alignment_stats import AlignmentStatsAnalyzer
from .coverage_stats import CoverageStatsAnalyzer
from .sequence_stats import SequenceStatsAnalyzer
from .insert_stats import InsertStatsAnalyzer
from .duplicate_stats import DuplicateStatsAnalyzer
from .utils import get_sample_name

class SampleStatsGenerator:
    """样品统计生成器|Sample Statistics Generator"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

        # 初始化各个分析器|Initialize analyzers
        self.alignment_analyzer = AlignmentStatsAnalyzer(config, logger, cmd_runner)
        self.coverage_analyzer = CoverageStatsAnalyzer(config, logger, cmd_runner)
        self.sequence_analyzer = SequenceStatsAnalyzer(config, logger, cmd_runner)
        self.insert_analyzer = InsertStatsAnalyzer(config, logger, cmd_runner)
        self.duplicate_analyzer = DuplicateStatsAnalyzer(config, logger, cmd_runner)

    def analyze_sample(self, bam_file: str) -> Dict[str, Any]:
        """分析单个样品|Analyze single sample"""
        sample_name = get_sample_name(bam_file)
        self.logger.info(f"开始分析样品|Starting sample analysis: {sample_name}")

        sample_stats = {
            'sample_name': sample_name,
            'bam_file': bam_file,
            'analysis_modules': {}
        }

        # 比对统计|Alignment statistics
        if not self.config.skip_alignment_stats:
            self.logger.info("分析比对统计|Analyzing alignment statistics")
            sample_stats['alignment_stats'] = self.alignment_analyzer.get_basic_alignment_stats(bam_file)
            sample_stats['analysis_modules']['alignment'] = True
        else:
            sample_stats['analysis_modules']['alignment'] = False

        # 覆盖度统计|Coverage statistics
        if not self.config.skip_coverage_stats:
            self.logger.info("分析覆盖度统计|Analyzing coverage statistics")
            sample_stats['coverage_stats'] = self.coverage_analyzer.get_coverage_stats(bam_file)
            sample_stats['analysis_modules']['coverage'] = True
        else:
            sample_stats['analysis_modules']['coverage'] = False

        # 序列特征统计|Sequence feature statistics
        if not self.config.skip_sequence_stats:
            self.logger.info("分析序列特征|Analyzing sequence features")
            sample_stats['sequence_stats'] = self.sequence_analyzer.get_sequence_stats(bam_file)
            sample_stats['analysis_modules']['sequence'] = True
        else:
            sample_stats['analysis_modules']['sequence'] = False

        # 插入片段统计|Insert size statistics
        if not self.config.skip_insert_stats:
            self.logger.info("分析插入片段|Analyzing insert sizes")
            sample_stats['insert_stats'] = self.insert_analyzer.get_insert_stats(bam_file)
            sample_stats['analysis_modules']['insert'] = True
        else:
            sample_stats['analysis_modules']['insert'] = False

        # 重复序列统计|Duplicate statistics
        if not self.config.skip_duplicate_stats:
            self.logger.info("分析重复序列|Analyzing duplicates")
            sample_stats['duplicate_stats'] = self.duplicate_analyzer.get_duplicate_stats(bam_file)
            sample_stats['analysis_modules']['duplicate'] = True
        else:
            sample_stats['analysis_modules']['duplicate'] = False

        self.logger.info(f"样品分析完成|Sample analysis completed: {sample_name}")
        return sample_stats
