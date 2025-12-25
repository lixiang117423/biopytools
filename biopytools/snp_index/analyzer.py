"""
SNP Index结果分析模块 | SNP Index Result Analysis Module
"""

import csv
import logging
import numpy as np
from collections import Counter, defaultdict
from typing import List, Dict, Tuple, Optional
from .config import SNPIndexConfig
from .utils import format_number, format_large_number, check_file_exists


class SNPIndexAnalyzer:
    """SNP Index结果分析器 | SNP Index Result Analyzer"""

    def __init__(self, result_file: str, config: Optional[SNPIndexConfig] = None):
        """
        初始化分析器 | Initialize analyzer

        Args:
            result_file: 结果文件路径 | Result file path
            config: 配置对象 | Configuration object
        """
        self.result_file = result_file
        self.config = config or SNPIndexConfig()
        self.logger = logging.getLogger(self.__class__.__name__)

        # 数据存储 | Data storage
        self.data = []
        self.sample_names = []
        self.statistics = {}
        self.extreme_sites = []
        self.potential_regions = []

    def load_data(self) -> bool:
        """
        加载结果数据 | Load result data

        Returns:
            bool: 加载是否成功 | Whether loading succeeded
        """
        self.logger.info("加载SNP index结果数据 | Loading SNP index result data")

        try:
            check_file_exists(self.result_file, "结果文件 | Result file")

            with open(self.result_file, 'r', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\t')
                header = next(reader)

                # 解析样本名称 | Parse sample names
                self.sample_names = self._parse_sample_names(header)
                self.logger.info(f"检测到样本 | Detected samples: {self.sample_names}")

                # 读取数据行 | Read data rows
                for row in reader:
                    if len(row) >= 11:
                        try:
                            parsed_row = self._parse_data_row(row)
                            if parsed_row:
                                self.data.append(parsed_row)
                        except (ValueError, IndexError) as e:
                            self.logger.warning(f"解析数据行时出错，跳过 | Error parsing data row, skipping: {str(e)}")

            self.logger.info(f"成功加载 | Successfully loaded: {len(self.data):,} 条记录")
            return len(self.data) > 0

        except Exception as e:
            self.logger.error(f"加载数据时出错 | Error loading data: {str(e)}")
            return False

    def _parse_sample_names(self, header: List[str]) -> List[str]:
        """
        从表头解析样本名称 | Parse sample names from header

        Args:
            header: 表头列表 | Header list

        Returns:
            list: 样本名称列表 | Sample names list
        """
        sample_names = []
        for col in header:
            if '_SNP_index' in col:
                sample_name = col.replace('_SNP_index', '')
                sample_names.append(sample_name)
                if len(sample_names) >= 2:  # 只需要前两个 | Only need first two
                    break
        return sample_names[:2] if sample_names else ['Sample1', 'Sample2']

    def _parse_data_row(self, row: List[str]) -> Optional[Dict]:
        """
        解析数据行 | Parse data row

        Args:
            row: 数据行列表 | Data row list

        Returns:
            dict: 解析后的数据字典 | Parsed data dictionary
        """
        if len(row) < 11:
            return None

        return {
            'chromosome': row[0],
            'position': int(row[1]),
            'reference': row[2],
            'alternative': row[3],
            'sample1_ref_depth': int(row[4]),
            'sample1_alt_depth': int(row[5]),
            'sample1_snp_index': float(row[6]),
            'sample2_ref_depth': int(row[7]),
            'sample2_alt_depth': int(row[8]),
            'sample2_snp_index': float(row[9]),
            'delta_snp_index': float(row[10])
        }

    def analyze(self) -> bool:
        """
        执行完整分析 | Perform complete analysis

        Returns:
            bool: 分析是否成功 | Whether analysis succeeded
        """
        self.logger.info("=" * 60)
        self.logger.info("开始SNP index结果分析 | Starting SNP index result analysis")
        self.logger.info("=" * 60)

        # 加载数据 | Load data
        if not self.load_data():
            return False

        try:
            # 计算基本统计 | Calculate basic statistics
            self._calculate_basic_statistics()

            # 查找极端位点 | Find extreme sites
            self._find_extreme_sites()

            # 查找潜在目标区域 | Find potential target regions
            self._find_potential_regions()

            # 输出分析结果 | Output analysis results
            self._print_analysis_results()

            self.logger.info("SNP index结果分析完成 | SNP index result analysis completed")
            return True

        except Exception as e:
            self.logger.error(f"分析过程中发生错误 | Error during analysis: {str(e)}", exc_info=True)
            return False

    def _calculate_basic_statistics(self) -> None:
        """计算基本统计信息 | Calculate basic statistics"""
        self.logger.info("计算基本统计信息 | Calculating basic statistics")

        if not self.data:
            return

        # 提取数值数组 | Extract numeric arrays
        delta_snp_indices = [d['delta_snp_index'] for d in self.data]
        snp_indices1 = [d['sample1_snp_index'] for d in self.data]
        snp_indices2 = [d['sample2_snp_index'] for d in self.data]

        # 计算统计量 | Calculate statistics
        self.statistics = {
            'total_snps': len(self.data),
            'delta_snp_index': {
                'mean': np.mean(delta_snp_indices),
                'median': np.median(delta_snp_indices),
                'std': np.std(delta_snp_indices),
                'min': np.min(delta_snp_indices),
                'max': np.max(delta_snp_indices)
            },
            'sample1_snp_index': {
                'mean': np.mean(snp_indices1),
                'median': np.median(snp_indices1)
            },
            'sample2_snp_index': {
                'mean': np.mean(snp_indices2),
                'median': np.median(snp_indices2)
            }
        }

        # 染色体分布统计 | Chromosome distribution statistics
        chromosomes = [d['chromosome'] for d in self.data]
        chr_counts = Counter(chromosomes)
        self.statistics['chromosome_distribution'] = dict(chr_counts.most_common())

    def _find_extreme_sites(self) -> None:
        """查找极端ΔSNP index位点 | Find extreme ΔSNP index sites"""
        self.logger.info(f"查找极端ΔSNP index位点 (|ΔSNP index| > {self.config.extreme_threshold}) | "
                        f"Finding extreme ΔSNP index sites (|ΔSNP index| > {self.config.extreme_threshold})")

        threshold = self.config.extreme_threshold
        self.extreme_sites = [
            d for d in self.data
            if abs(d['delta_snp_index']) > threshold
        ]

        self.logger.info(f"找到 | Found: {len(self.extreme_sites)} 个极端位点")

    def _find_potential_regions(self) -> None:
        """查找潜在的目标区域 | Find potential target regions"""
        self.logger.info(f"查找潜在目标区域 (连续{self.config.min_region_snps}个以上 |ΔSNP index| > {self.config.region_threshold}) | "
                        f"Finding potential target regions (>{self.config.min_region_snps} consecutive sites with |ΔSNP index| > {self.config.region_threshold})")

        # 按染色体分组并排序 | Group by chromosome and sort
        chr_data = defaultdict(list)
        for d in self.data:
            chr_data[d['chromosome']].append(d)

        # 对每个染色体找区域 | Find regions for each chromosome
        self.potential_regions = []
        for chrom, sites in chr_data.items():
            sites.sort(key=lambda x: x['position'])
            regions = self._find_regions_in_chromosome(sites)
            self.potential_regions.extend(regions)

        self.logger.info(f"找到 | Found: {len(self.potential_regions)} 个潜在目标区域")

    def _find_regions_in_chromosome(self, sites: List[Dict]) -> List[Dict]:
        """
        在单个染色体中查找连续区域 | Find consecutive regions in a single chromosome

        Args:
            sites: 位点列表 | Site list

        Returns:
            list: 区域列表 | Region list
        """
        regions = []
        current_region = None

        threshold = self.config.region_threshold
        min_snps = self.config.min_region_snps
        max_gap = self.config.max_region_gap

        for site in sites:
            if abs(site['delta_snp_index']) > threshold:
                if current_region is None:
                    # 开始新区域 | Start new region
                    current_region = {
                        'chromosome': site['chromosome'],
                        'start': site['position'],
                        'end': site['position'],
                        'count': 1,
                        'mean_delta': site['delta_snp_index'],
                        'sites': [site]
                    }
                elif (site['chromosome'] == current_region['chromosome'] and
                      site['position'] - current_region['end'] <= max_gap):
                    # 扩展当前区域 | Extend current region
                    current_region['end'] = site['position']
                    current_region['count'] += 1
                    current_region['sites'].append(site)
                    # 更新平均ΔSNP index | Update mean ΔSNP index
                    current_region['mean_delta'] = np.mean([s['delta_snp_index'] for s in current_region['sites']])
                else:
                    # 保存当前区域，开始新区域 | Save current region, start new region
                    if current_region['count'] >= min_snps:
                        regions.append(current_region)
                    current_region = {
                        'chromosome': site['chromosome'],
                        'start': site['position'],
                        'end': site['position'],
                        'count': 1,
                        'mean_delta': site['delta_snp_index'],
                        'sites': [site]
                    }
            else:
                # 保存当前区域 | Save current region
                if current_region and current_region['count'] >= min_snps:
                    regions.append(current_region)
                current_region = None

        # 检查最后一个区域 | Check last region
        if current_region and current_region['count'] >= min_snps:
            regions.append(current_region)

        return regions

    def _print_analysis_results(self) -> None:
        """打印分析结果 | Print analysis results"""
        self.logger.info("=" * 60)
        self.logger.info("SNP Index分析结果 | SNP Index Analysis Results")
        self.logger.info("=" * 60)

        # 基本统计 | Basic statistics
        self.logger.info(f"总SNP数量 | Total SNPs: {format_large_number(self.statistics['total_snps'])}")

        delta_stats = self.statistics['delta_snp_index']
        self.logger.info(f"ΔSNP index统计 | ΔSNP index statistics:")
        self.logger.info(f"  平均值 | Mean: {format_number(delta_stats['mean'])}")
        self.logger.info(f"  中位数 | Median: {format_number(delta_stats['median'])}")
        self.logger.info(f"  标准差 | Std: {format_number(delta_stats['std'])}")
        self.logger.info(f"  最小值 | Min: {format_number(delta_stats['min'])}")
        self.logger.info(f"  最大值 | Max: {format_number(delta_stats['max'])}")

        if len(self.sample_names) >= 2:
            sample1_stats = self.statistics['sample1_snp_index']
            sample2_stats = self.statistics['sample2_snp_index']
            self.logger.info(f"{self.sample_names[0]} SNP index统计 | {self.sample_names[0]} SNP index statistics:")
            self.logger.info(f"  平均值 | Mean: {format_number(sample1_stats['mean'])}")
            self.logger.info(f"  中位数 | Median: {format_number(sample1_stats['median'])}")

            self.logger.info(f"{self.sample_names[1]} SNP index统计 | {self.sample_names[1]} SNP index statistics:")
            self.logger.info(f"  平均值 | Mean: {format_number(sample2_stats['mean'])}")
            self.logger.info(f"  中位数 | Median: {format_number(sample2_stats['median'])}")

        # 染色体分布 | Chromosome distribution
        self.logger.info("\n染色体分布 | Chromosome distribution:")
        for chrom, count in list(self.statistics['chromosome_distribution'].items())[:10]:
            self.logger.info(f"  {chrom}: {format_large_number(count)} SNPs")

        # 极端位点 | Extreme sites
        if self.extreme_sites:
            self.logger.info(f"\n极端ΔSNP index位点 | Extreme ΔSNP index sites (|ΔSNP index| > {self.config.extreme_threshold}):")
            self.logger.info(f"  总数 | Total: {len(self.extreme_sites)}")
            # 显示前10个 | Show top 10
            for site in self.extreme_sites[:10]:
                self.logger.info(f"  {site['chromosome']}:{site['position']} "
                               f"Ref:{site['reference']} Alt:{site['alternative']} "
                               f"ΔSNP_index={format_number(site['delta_snp_index'])}")
            if len(self.extreme_sites) > 10:
                self.logger.info(f"  ... 还有 {len(self.extreme_sites) - 10} 个极端位点")

        # 潜在目标区域 | Potential target regions
        if self.potential_regions:
            self.logger.info(f"\n潜在目标区域 | Potential target regions:")
            self.logger.info(f"  总数 | Total: {len(self.potential_regions)}")
            for region in self.potential_regions:
                length = region['end'] - region['start']
                self.logger.info(f"  {region['chromosome']}: {region['start']}-{region['end']} "
                               f"({length:,} bp, {region['count']} SNPs, "
                               f"mean ΔSNP index={format_number(region['mean_delta'])})")

    def get_statistics(self) -> Dict:
        """
        获取统计结果 | Get statistics results

        Returns:
            dict: 统计结果字典 | Statistics results dictionary
        """
        return {
            'basic_statistics': self.statistics,
            'extreme_sites': self.extreme_sites,
            'potential_regions': self.potential_regions,
            'sample_names': self.sample_names
        }

    def export_results(self, output_prefix: str) -> bool:
        """
        导出分析结果 | Export analysis results

        Args:
            output_prefix: 输出文件前缀 | Output file prefix

        Returns:
            bool: 导出是否成功 | Whether export succeeded
        """
        try:
            # 导出极端位点 | Export extreme sites
            if self.extreme_sites:
                extreme_file = f"{output_prefix}_extreme_sites.tsv"
                with open(extreme_file, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow(['Chromosome', 'Position', 'Reference', 'Alternative',
                                   'ΔSNP_index', 'Sample1_SNP_index', 'Sample2_SNP_index'])
                    for site in self.extreme_sites:
                        writer.writerow([
                            site['chromosome'], site['position'],
                            site['reference'], site['alternative'],
                            format_number(site['delta_snp_index']),
                            format_number(site['sample1_snp_index']),
                            format_number(site['sample2_snp_index'])
                        ])
                self.logger.info(f"极端位点已导出 | Extreme sites exported: {extreme_file}")

            # 导出潜在区域 | Export potential regions
            if self.potential_regions:
                region_file = f"{output_prefix}_potential_regions.tsv"
                with open(region_file, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow(['Chromosome', 'Start', 'End', 'Length', 'SNP_Count', 'Mean_ΔSNP_index'])
                    for region in self.potential_regions:
                        length = region['end'] - region['start']
                        writer.writerow([
                            region['chromosome'], region['start'], region['end'],
                            length, region['count'], format_number(region['mean_delta'])
                        ])
                self.logger.info(f"潜在区域已导出 | Potential regions exported: {region_file}")

            return True

        except Exception as e:
            self.logger.error(f"导出结果时出错 | Error exporting results: {str(e)}")
            return False