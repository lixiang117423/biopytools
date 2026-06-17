"""
HiC-Pro质量控制评估核心计算模块|HiC-Pro QC Assessment Core Calculation Module
"""

import os
from pathlib import Path
from typing import Dict, Any, Optional


class HiCProQCCalculator:
    """HiC-Pro质量控制评估计算器|HiC-Pro QC Assessment Calculator"""

    def __init__(self, config, logger):
        """初始化计算器|Initialize calculator

        Args:
            config: HiCProQCConfig配置对象|HiCProQCConfig object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger
        self.stats_data = {}

    def parse_mmapstat(self, mmapstat_file: str) -> Dict[str, Any]:
        """解析合并的mapping统计文件|Parse merged mapping statistics file

        Args:
            mmapstat_file: .mmapstat文件路径|.mmapstat file path

        Returns:
            Dict[str, Any]: Mapping统计数据|Mapping statistics
        """
        self.logger.info(f"解析mapping统计文件|Parsing mapping statistics: {mmapstat_file}")

        if not os.path.exists(mmapstat_file):
            self.logger.warning(f"未找到mapping统计文件|Mapping statistics file not found: {mmapstat_file}")
            return {}

        stats = {}
        try:
            with open(mmapstat_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split('\t')
                    if len(parts) >= 2:
                        key = parts[0]
                        value = parts[1]

                        # 转换为数字|Convert to number
                        try:
                            value = int(value)
                        except ValueError:
                            try:
                                value = float(value)
                            except ValueError:
                                pass

                        stats[key] = value

            self.logger.info(f"Mapped statistics loaded: {len(stats)} fields")
            return stats

        except Exception as e:
            self.logger.error(f"解析mapping统计文件失败|Failed to parse mapping statistics: {str(e)}")
            return {}

    def parse_mpairstat(self, mpairstat_file: str) -> Dict[str, Any]:
        """解析合并的pairing统计文件|Parse merged pairing statistics file

        Args:
            mpairstat_file: .mpairstat文件路径|.mpairstat file path

        Returns:
            Dict[str, Any]: Pairing统计数据|Pairing statistics
        """
        self.logger.info(f"解析pairing统计文件|Parsing pairing statistics: {mpairstat_file}")

        if not os.path.exists(mpairstat_file):
            self.logger.warning(f"未找到pairing统计文件|Pairing statistics file not found: {mpairstat_file}")
            return {}

        stats = {}
        try:
            with open(mpairstat_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split('\t')
                    if len(parts) >= 2:
                        key = parts[0]
                        value = parts[1]

                        # 转换为数字|Convert to number
                        try:
                            value = int(value)
                        except ValueError:
                            try:
                                value = float(value)
                            except ValueError:
                                pass

                        stats[key] = value

            self.logger.info(f"Pairing statistics loaded: {len(stats)} fields")
            return stats

        except Exception as e:
            self.logger.error(f"解析pairing统计文件失败|Failed to parse pairing statistics: {str(e)}")
            return {}

    def parse_mrsstat(self, mrsstat_file: str) -> Dict[str, Any]:
        """解析合并的RS统计文件（Valid pairs过滤统计）|Parse merged RS statistics file

        Args:
            mrsstat_file: .mRSstat文件路径|.mRSstat file path

        Returns:
            Dict[str, Any]: RS统计数据|RS statistics
        """
        self.logger.info(f"解析RS统计文件|Parsing RS statistics: {mrsstat_file}")

        if not os.path.exists(mrsstat_file):
            self.logger.warning(f"未找到RS统计文件|RS statistics file not found: {mrsstat_file}")
            return {}

        stats = {}
        try:
            with open(mrsstat_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split('\t')
                    if len(parts) >= 2:
                        key = parts[0]
                        value = parts[1]

                        # 转换为数字|Convert to number
                        try:
                            value = int(value)
                        except ValueError:
                            try:
                                value = float(value)
                            except ValueError:
                                pass

                        stats[key] = value

            self.logger.info(f"RS statistics loaded: {len(stats)} fields")
            return stats

        except Exception as e:
            self.logger.error(f"解析RS统计文件失败|Failed to parse RS statistics: {str(e)}")
            return {}

    def calculate_metrics(self, mapstat: Dict, pairstat: Dict, rsstat: Dict) -> Dict[str, Any]:
        """计算质量指标|Calculate quality metrics

        Args:
            mapstat: Mapping统计|Mapping statistics
            pairstat: Pairing统计|Pairing statistics
            rsstat: RS统计|RS statistics

        Returns:
            Dict[str, Any]: 计算出的质量指标|Calculated quality metrics
        """
        self.logger.info("计算质量指标|Calculating quality metrics")

        metrics = {}

        # 1. Mapping指标|Mapping metrics
        if mapstat:
            # 提取R1和R2的统计|Extract R1 and R2 statistics
            total_r1 = mapstat.get('total_R1', 0)
            total_r2 = mapstat.get('total_R2', 0)
            mapped_r1 = mapstat.get('mapped_R1', 0)
            mapped_r2 = mapstat.get('mapped_R2', 0)

            total_reads = total_r1 + total_r2
            total_mapped = mapped_r1 + mapped_r2

            if total_reads > 0:
                metrics['total_reads'] = total_reads
                metrics['total_mapped'] = total_mapped
                metrics['mapping_rate'] = (total_mapped / total_reads) * 100

            # 计算unique mapping rate（假设global即为unique）|Calculate unique mapping rate
            global_r1 = mapstat.get('global_R1', 0)
            global_r2 = mapstat.get('global_R2', 0)
            total_global = global_r1 + global_r2

            if total_reads > 0:
                metrics['total_global_unique'] = total_global
                metrics['unique_mapping_rate'] = (total_global / total_reads) * 100

        # 2. Valid pairs指标|Valid pairs metrics
        if rsstat:
            valid_pairs = rsstat.get('Valid_interaction_pairs', 0)
            dangling_ends = rsstat.get('Dangling_end_pairs', 0)
            self_circle = rsstat.get('Self_Cycle_pairs', 0)
            religation = rsstat.get('Religation_pairs', 0)

            # 计算total pairs（valid + invalid）|Calculate total pairs
            total_pairs = (
                valid_pairs + dangling_ends + self_circle + religation +
                rsstat.get('Single-end_pairs', 0) +
                rsstat.get('Filtered_pairs', 0) +
                rsstat.get('Dumped_pairs', 0)
            )

            if total_pairs > 0:
                metrics['total_pairs'] = total_pairs
                metrics['valid_pairs'] = valid_pairs
                metrics['valid_pairs_rate'] = (valid_pairs / total_pairs) * 100

                metrics['dangling_ends'] = dangling_ends
                metrics['dangling_ends_rate'] = (dangling_ends / total_pairs) * 100

                metrics['self_ligation'] = self_circle
                metrics['self_ligation_rate'] = (self_circle / total_pairs) * 100

                metrics['religation'] = religation
                metrics['religation_rate'] = (religation / total_pairs) * 100

        # 3. Cis/Trans指标|Cis/Trans metrics
        # 需要从valid pairs文件中计算（这里暂时跳过，需要额外处理）
        # 暂时使用默认值
        metrics['cis_trans_ratio'] = 0.0

        # 4. PCR重复指标|PCR duplication metrics
        # HiC-Pro会在去重前后生成统计，暂时设置为0
        metrics['duplication_rate'] = 0.0

        return metrics

    def assess_quality(self, metrics: Dict[str, Any]) -> Dict[str, Any]:
        """评估数据质量|Assess data quality

        Args:
            metrics: 质量指标|Quality metrics

        Returns:
            Dict[str, Any]: 评估结果|Assessment results
        """
        self.logger.info("评估数据质量|Assessing data quality")

        thresholds = self.config.get_thresholds()
        results = {
            'metrics': metrics,
            'thresholds': thresholds,
            'assessment': {},
            'passed': True
        }

        # 定义评估指标|Define assessment metrics
        assessment_metrics = [
            # (metric_key, threshold_key, comparison_type)
            ('mapping_rate', 'min_mapping_rate', 'ge'),
            ('unique_mapping_rate', 'min_unique_rate', 'ge'),
            ('valid_pairs_rate', 'min_valid_pairs_rate', 'ge'),
            ('dangling_ends_rate', 'max_dangling_ends_rate', 'le'),
            ('self_ligation_rate', 'max_self_ligation_rate', 'le'),
            ('religation_rate', 'max_religation_rate', 'le'),
            ('duplication_rate', 'max_duplication_rate', 'le'),
            ('cis_trans_ratio', 'min_cis_trans_ratio', 'ge'),
        ]

        for metric_key, threshold_key, comparison in assessment_metrics:
            value = metrics.get(metric_key, 0)
            threshold = thresholds.get(threshold_key, 0)

            # 判断是否通过|Check if passed
            if comparison == 'ge':
                passed = value >= threshold
                if 'rate' in metric_key or 'ratio' in metric_key:
                    desc = f"{value:.2f}% (阈值|threshold >= {threshold}%)"
                else:
                    desc = f"{value:.2f} (阈值|threshold >= {threshold})"
            else:  # 'le'
                passed = value <= threshold
                if 'rate' in metric_key or 'ratio' in metric_key:
                    desc = f"{value:.2f}% (阈值|threshold <= {threshold}%)"
                else:
                    desc = f"{value:.2f} (阈值|threshold <= {threshold})"

            results['assessment'][metric_key] = {
                'value': value,
                'threshold': threshold,
                'passed': passed,
                'description': desc
            }

            if not passed:
                results['passed'] = False

        return results

    def calculate(self) -> Dict[str, Any]:
        """执行完整的质量评估流程|Execute complete quality assessment pipeline

        Returns:
            Dict[str, Any]: 评估结果|Assessment results
        """
        # 查找统计文件|Find statistics files
        sample_name = self.config.sample_name
        bowtie_dir = Path(self.config.hicpro_dir) / 'bowtie_results' / 'bwt2'
        hic_data_dir = Path(self.config.hicpro_dir) / 'hic_results' / 'data'

        # 解析mmapstat|Parse mmapstat
        mmapstat_file = bowtie_dir / f'{sample_name}.mmapstat'
        mapstat = self.parse_mmapstat(str(mmapstat_file))

        # 解析mpairstat|Parse mpairstat
        mpairstat_file = bowtie_dir / f'{sample_name}.mpairstat'
        pairstat = self.parse_mpairstat(str(mpairstat_file))

        # 解析mRSstat|Parse mRSstat
        mrsstat_file = hic_data_dir / f'{sample_name}.mRSstat'
        rsstat = self.parse_mrsstat(str(mrsstat_file))

        # 计算质量指标|Calculate quality metrics
        metrics = self.calculate_metrics(mapstat, pairstat, rsstat)

        # 评估质量|Assess quality
        results = self.assess_quality(metrics)

        return results

    def format_number(self, num: int) -> str:
        """格式化数字，添加千位分隔符|Format number with thousands separator

        Args:
            num: 要格式化的数字|Number to format

        Returns:
            str: 格式化后的字符串|Formatted string
        """
        return f"{num:,}"
