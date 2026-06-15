"""
Pi结果合并模块|Pi Results Merger Module
合并pi计算结果到一个汇总表
Merge pi calculation results into a summary table
"""

from pathlib import Path
from typing import Dict, List, Tuple

from .config import PiConfig
from .utils import PiRow, format_number


class PiResultsMerger:
    """Pi结果合并器|Pi Results Merger"""

    def __init__(self, config: PiConfig, logger):
        self.config = config
        self.logger = logger

    def merge_and_write(self, vcftools_results: List[PiRow]) -> str:
        """
        合并结果并写入TSV文件|Merge results and write to TSV file

        Args:
            vcftools_results: vcftools pi结果列表|vcftools pi results

        Returns:
            输出文件路径|Output file path
        """
        output_file = str(self.config.output_path / 'pi_merged.tsv')

        is_windowed = self.config.window_size is not None

        if not vcftools_results:
            self.logger.error("没有任何计算结果|No calculation results available")
            return output_file

        self.logger.info("开始合并结果|Starting result merge")

        if is_windowed:
            merged = self._merge_windowed(vcftools_results)
        else:
            merged = self._merge_genome_wide(vcftools_results)

        # 写入TSV文件|Write TSV file
        self._write_tsv(merged, output_file, is_windowed)

        self.logger.info(
            f"合并完成，共{len(merged)}行结果写入|"
            f"Merge completed, {len(merged)} rows written to: {output_file}"
        )

        return output_file

    def _merge_windowed(self, vcftools_results: List[PiRow]) -> List[dict]:
        """
        合并窗口模式结果|Merge windowed results

        按键: (population, chromosome, window_start, window_end) 去重合并
        Key: (population, chromosome, window_start, window_end)
        """
        merged_data: Dict[Tuple, dict] = {}

        for row in vcftools_results:
            key = (row.population, row.chromosome, row.window_start, row.window_end)
            if key not in merged_data:
                merged_data[key] = {
                    'population': row.population,
                    'chromosome': row.chromosome,
                    'window_start': row.window_start,
                    'window_end': row.window_end,
                    'pi': row.pi_value,
                    'n_sites': row.n_sites,
                }
            else:
                merged_data[key]['pi'] = row.pi_value
                merged_data[key]['n_sites'] = row.n_sites

        # 排序|Sort by population, chromosome, window_start
        sorted_results = sorted(
            merged_data.values(),
            key=lambda x: (x['population'], x['chromosome'], x['window_start'] or 0)
        )

        return sorted_results

    def _merge_genome_wide(self, vcftools_results: List[PiRow]) -> List[dict]:
        """
        合并全基因组模式结果|Merge genome-wide results

        按键: (population, chromosome) 去重合并
        Key: (population, chromosome)
        """
        merged_data: Dict[Tuple, dict] = {}

        for row in vcftools_results:
            key = (row.population, row.chromosome)
            if key not in merged_data:
                merged_data[key] = {
                    'population': row.population,
                    'chromosome': row.chromosome,
                    'pi': row.pi_value,
                    'n_sites': row.n_sites,
                }
            else:
                merged_data[key]['pi'] = row.pi_value
                merged_data[key]['n_sites'] = row.n_sites

        # 排序|Sort by population, chromosome
        sorted_results = sorted(
            merged_data.values(),
            key=lambda x: (x['population'], x['chromosome'])
        )

        return sorted_results

    def _write_tsv(self, merged: List[dict], output_file: str, is_windowed: bool):
        """
        写入TSV文件|Write TSV file
        """
        with open(output_file, 'w') as f:
            # 写入header|Write header
            if is_windowed:
                header_parts = ['population', 'chromosome', 'window_start', 'window_end']
            else:
                header_parts = ['population', 'chromosome']

            header_parts.append('pi')
            header_parts.append('n_sites')

            f.write('\t'.join(header_parts) + '\n')

            # 写入数据行|Write data rows
            for row in merged:
                if is_windowed:
                    parts = [
                        row['population'],
                        row['chromosome'],
                        str(row['window_start']),
                        str(row['window_end']),
                    ]
                else:
                    parts = [
                        row['population'],
                        row['chromosome'],
                    ]

                pi_val = format_number(row['pi']) if row['pi'] is not None else 'NA'
                parts.append(pi_val)
                parts.append(str(row['n_sites']))

                f.write('\t'.join(parts) + '\n')
