"""
NCBI分类学注释统计模块|NCBI Taxonomy Annotation Statistics Module
"""

import os
from collections import Counter, defaultdict
from typing import Dict, List, Tuple
from .config import NCBITaxoConfig
from .utils import parse_lineage, format_number


class TaxonomyStatsCalculator:
    """分类学统计计算器|Taxonomy Statistics Calculator"""

    def __init__(self, config: NCBITaxoConfig, logger):
        self.config = config
        self.logger = logger

    def calculate_statistics(
        self,
        accession_lineage: Dict[str, str],
        blast_hits: Dict[str, int]
    ) -> Dict[str, List[Tuple]]:
        """计算统计信息|Calculate statistics

        Args:
            accession_lineage: accession到lineage的映射
            blast_hits: accession到hit次数的映射

        Returns:
            Dict[str, List[Tuple]]: 每个层级的统计结果
                {
                    'genus': [(name, count, percentage), ...],
                    'species': [(name, count, percentage), ...]
                }
        """
        self.logger.info("开始计算统计信息|Starting statistics calculation")

        # 检查是否有有效数据|Check if there is valid data
        if not accession_lineage:
            self.logger.warning("没有有效的lineage数据，跳过统计|No valid lineage data, skipping statistics")
            return {}

        stats = {}

        # 对每个统计层级进行统计|Calculate statistics for each level
        for level in self.config.stats_by:
            if self.config.stats_target in ['unique_accessions', 'both']:
                # 统计唯一accession|Count unique accessions
                level_stats_accessions = self._count_by_level(accession_lineage, level, count_type='accessions')
                stats[f'{level}_accessions'] = level_stats_accessions

            if self.config.stats_target in ['blast_hits', 'both']:
                # 统计blast hits|Count blast hits
                if blast_hits:
                    level_stats_hits = self._count_by_level(accession_lineage, level, blast_hits, count_type='hits')
                    stats[f'{level}_hits'] = level_stats_hits

        self.logger.info("统计信息计算完成|Statistics calculation completed")
        return stats

    def _count_by_level(
        self,
        accession_lineage: Dict[str, str],
        level: str,
        hit_counts: Dict[str, int] = None,
        count_type: str = 'accessions'
    ) -> List[Tuple[str, int, float]]:
        """按指定层级统计|Count by specified level

        Args:
            accession_lineage: accession到lineage的映射
            level: 统计层级 (genus/species等)
            hit_counts: accession到hit次数的映射（用于统计blast hits）
            count_type: 统计类型 ('accessions' 或 'hits')

        Returns:
            List[Tuple[str, int, float]]: [(名称, 数量, 百分比), ...] 按数量降序排列
        """
        counter = Counter()

        for accession, lineage in accession_lineage.items():
            # 从lineage中提取指定层级的名称
            name = parse_lineage(lineage, level)

            if count_type == 'accessions':
                # 统计唯一accession数量
                counter[name] += 1
            elif count_type == 'hits' and hit_counts:
                # 统计blast hit次数
                hits = hit_counts.get(accession, 0)
                counter[name] += hits

        # 计算总数和百分比
        total = sum(counter.values())

        # 转换为列表并按数量降序排序
        result = []
        for name, count in counter.most_common():
            percentage = (count / total * 100) if total > 0 else 0
            result.append((name, count, percentage))

        return result

    def format_output(
        self,
        stats: Dict[str, List[Tuple]],
        accession_lineage: Dict[str, str],
        blast_hits: Dict[str, int]
    ) -> str:
        """格式化统计输出|Format statistics output

        返回格式3: 汇总表格
        Level    Level_Name    Count    Percentage
        Genus    Phytophthora  743      45.23%
        ...
        """
        lines = []
        lines.append("=" * 80)
        lines.append("NCBI分类学注释统计结果|NCBI Taxonomy Annotation Statistics")
        lines.append("=" * 80)
        lines.append("")

        # 添加总体统计信息|Add overall statistics
        total_accessions = len(accession_lineage)
        total_hits = sum(blast_hits.values()) if blast_hits else 0

        lines.append("总体统计|Overall Statistics:")
        lines.append("-" * 80)
        lines.append(f"总accession数量|Total unique accessions: {total_accessions}")

        if blast_hits:
            lines.append(f"总BLAST hit次数|Total BLAST hits: {format_number(total_hits)}")
            if total_accessions > 0:
                avg_hits = total_hits / total_accessions
                lines.append(f"平均每accession hit次数|Average hits per accession: {avg_hits:.2f}")

        lines.append("")
        lines.append("=" * 80)
        lines.append("")

        # 输出各层级统计|Output statistics for each level
        for level in self.config.stats_by:
            lines.append(f"{level.upper()} 统计|{level.upper()} Statistics")
            lines.append("=" * 80)
            lines.append("")

            # 输出accessions统计|Output accessions statistics
            if f'{level}_accessions' in stats:
                lines.append(f"  [1] 基于唯一Accessions数量统计|By Unique Accessions Count")
                lines.append("  " + "-" * 76)
                lines.append(f"  {'Level':<10}{'Level_Name':<48}{'Count':>12}{'Percentage':>12}")
                lines.append("  " + "-" * 76)

                stats_accessions = stats[f'{level}_accessions']

                # 分离>=1%和<1%的项|Separate items >=1% and <1%
                major_items = []
                other_count = 0
                other_percentage = 0.0

                for name, count, percentage in stats_accessions:
                    if percentage >= 1.0:
                        major_items.append((name, count, percentage))
                    else:
                        other_count += count
                        other_percentage += percentage

                # 输出主要项（>=1%）|Output major items (>=1%)
                level_label = level.capitalize()
                for name, count, percentage in major_items:
                    display_name = name[:47] + '...' if len(name) > 50 else name
                    lines.append(f"  {level_label:<10}{display_name:<48}{count:>12}{percentage:>11.1f}%")

                # 输出Other项（<1%的汇总）|Output Other item (summary of <1%)
                if other_count > 0:
                    lines.append(f"  {level_label:<10}{'Other (< 1% each)':<48}{other_count:>12}{other_percentage:>11.1f}%")

                lines.append("")

            # 输出blast hits统计|Output blast hits statistics
            if f'{level}_hits' in stats:
                lines.append(f"  [2] 基于BLAST Hits次数统计|By BLAST Hits Count")
                lines.append("  " + "-" * 76)
                lines.append(f"  {'Level':<10}{'Level_Name':<48}{'Count':>12}{'Percentage':>12}")
                lines.append("  " + "-" * 76)

                stats_hits = stats[f'{level}_hits']

                # 分离>=1%和<1%的项|Separate items >=1% and <1%
                major_items = []
                other_count = 0
                other_percentage = 0.0

                for name, count, percentage in stats_hits:
                    if percentage >= 1.0:
                        major_items.append((name, count, percentage))
                    else:
                        other_count += count
                        other_percentage += percentage

                # 输出主要项（>=1%）|Output major items (>=1%)
                level_label = level.capitalize()
                for name, count, percentage in major_items:
                    display_name = name[:47] + '...' if len(name) > 50 else name
                    lines.append(f"  {level_label:<10}{display_name:<48}{format_number(count):>12}{percentage:>11.1f}%")

                # 输出Other项（<1%的汇总）|Output Other item (summary of <1%)
                if other_count > 0:
                    lines.append(f"  {level_label:<10}{'Other (< 1% each)':<48}{format_number(other_count):>12}{other_percentage:>11.1f}%")

            lines.append("")
            lines.append("=" * 80)
            lines.append("")

        return "\n".join(lines)

    def write_statistics(
        self,
        stats: Dict[str, List[Tuple]],
        accession_lineage: Dict[str, str],
        blast_hits: Dict[str, int]
    ) -> str:
        """写入统计文件|Write statistics file

        Returns:
            str: 主要统计文件路径
        """
        self.logger.info("开始写入统计文件|Starting to write statistics file")

        # 确定输出文件路径|Determine output file path
        if self.config.stats_output == 'csv':
            stats_file = f"{self.config.output_prefix}.statistics.csv"
            content = self._format_csv(stats, accession_lineage, blast_hits)
        else:
            stats_file = f"{self.config.output_prefix}.statistics.txt"
            content = self.format_output(stats, accession_lineage, blast_hits)

        # 写入主要统计文件|Write main statistics file
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write(content)

        self.logger.info(f"统计文件已保存|Statistics file saved: {stats_file}")

        return stats_file

    def _format_csv(
        self,
        stats: Dict[str, List[Tuple]],
        accession_lineage: Dict[str, str],
        blast_hits: Dict[str, int]
    ) -> str:
        """格式化为CSV输出|Format as CSV output"""
        lines = []

        # CSV header
        lines.append("Statistic_Type,Level,Level_Name,Count,Percentage")

        # 总体统计|Overall statistics
        total_accessions = len(accession_lineage)
        lines.append(f"OVERALL,TOTAL,Unique_Accessions,{total_accessions},")

        if blast_hits:
            total_hits = sum(blast_hits.values())
            lines.append(f"OVERALL,TOTAL,BLAST_HITS,{total_hits},")
            if total_accessions > 0:
                avg_hits = total_hits / total_accessions
                lines.append(f"OVERALL,TOTAL,Average_Hits_Per_Accession,{avg_hits:.2f},")

        # 各层级统计|Statistics by level
        for level in self.config.stats_by:
            # Accessions统计
            if f'{level}_accessions' in stats:
                for name, count, percentage in stats[f'{level}_accessions']:
                    lines.append(f"UNIQUE_ACCESSIONS,{level.capitalize()},\"{name}\",{count},{percentage:.2f}")

            # BLAST hits统计
            if f'{level}_hits' in stats:
                for name, count, percentage in stats[f'{level}_hits']:
                    lines.append(f"BLAST_HITS,{level.capitalize()},\"{name}\",{count},{percentage:.2f}")

        return "\n".join(lines)

    def generate_summary_report(
        self,
        accession_lineage: Dict[str, str],
        blast_hits: Dict[str, int],
        stats_file: str
    ):
        """生成汇总报告|Generate summary report"""
        self.logger.info("生成汇总报告|Generating summary report")

        total_accessions = len(accession_lineage)
        total_hits = sum(blast_hits.values()) if blast_hits else 0

        self.logger.info("=" * 60)
        self.logger.info("处理汇总|Processing Summary")
        self.logger.info("=" * 60)
        self.logger.info(f"总accession数量|Total accessions: {total_accessions}")

        if blast_hits:
            self.logger.info(f"总BLAST hit次数|Total BLAST hits: {format_number(total_hits)}")
            self.logger.info(f"平均每accession hit次数|Average hits per accession: {total_hits/total_accessions:.2f}")

        self.logger.info(f"统计文件已保存|Statistics saved: {stats_file}")
        self.logger.info("=" * 60)
