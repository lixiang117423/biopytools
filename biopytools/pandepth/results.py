"""
PanDepth结果合并模块|PanDepth Results Merging Module
"""

import os
import gzip
import pandas as pd
from pathlib import Path


class PanDepthResultsMerger:
    """PanDepth结果合并器|PanDepth Results Merger"""

    def __init__(self, logger):
        self.logger = logger

    def _read_gzipped_file(self, file_path: str) -> list:
        """读取gzip压缩的统计文件|Read gzipped statistics file

        Args:
            file_path: 文件路径|File path

        Returns:
            list: 文件行列表|List of file lines
        """
        try:
            with gzip.open(file_path, 'rt') as f:
                return f.readlines()
        except Exception as e:
            self.logger.error(f"读取文件失败|Failed to read file: {file_path}")
            self.logger.error(f"错误|Error: {e}")
            return None

    def _find_stat_files(self, output_dir: str, file_extension: str) -> list:
        """查找统计文件|Find statistics files

        Args:
            output_dir: 输出目录|Output directory
            file_extension: 文件扩展名|File extension (e.g., '.chr.stat.gz')

        Returns:
            list: 文件路径列表|List of file paths
        """
        stat_files = []
        for file in Path(output_dir).glob(f'*{file_extension}'):
            stat_files.append(str(file))
        return sorted(stat_files)

    def _parse_chr_stat_file(self, lines: list, sample_name: str) -> dict:
        """解析染色体统计文件|Parse chromosome statistics file

        Args:
            lines: 文件行列表|List of file lines
            sample_name: 样本名称|Sample name

        Returns:
            dict: 染色体统计数据|Chromosome statistics data
        """
        data = {}

        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            cols = line.split('\t')
            if len(cols) >= 6:
                chr_name = cols[0]
                coverage_pct = float(cols[4])
                mean_depth = float(cols[5])

                if chr_name not in data:
                    data[chr_name] = {}

                data[chr_name][f'{sample_name}_coverage'] = coverage_pct
                data[chr_name][f'{sample_name}_mean_depth'] = mean_depth

        return data

    def _parse_gene_stat_file(self, lines: list, sample_name: str) -> dict:
        """解析基因统计文件|Parse gene statistics file

        Args:
            lines: 文件行列表|List of file lines
            sample_name: 样本名称|Sample name

        Returns:
            dict: 基因统计数据|Gene statistics data
        """
        data = {}

        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            cols = line.split('\t')
            if len(cols) >= 8:
                gene_id = cols[3]
                coverage_pct = float(cols[6])
                mean_depth = float(cols[7])

                if gene_id not in data:
                    data[gene_id] = {}

                data[gene_id][f'{sample_name}_coverage'] = coverage_pct
                data[gene_id][f'{sample_name}_mean_depth'] = mean_depth

        return data

    def merge_chr_statistics(self, output_dir: str, output_file: str = 'merged_chr_statistics.txt') -> bool:
        """合并染色体统计结果|Merge chromosome statistics results

        Args:
            output_dir: 输出目录|Output directory
            output_file: 输出文件名|Output file name

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("开始合并染色体统计结果|Starting chromosome statistics merging")

        # 查找所有染色体统计文件|Find all chromosome statistics files
        stat_files = self._find_stat_files(output_dir, '.chr.stat.gz')

        if not stat_files:
            self.logger.warning("未找到染色体统计文件|No chromosome statistics files found")
            return False

        self.logger.info(f"找到 {len(stat_files)} 个染色体统计文件|Found {len(stat_files)} chromosome statistics files")

        # 解析所有文件|Parse all files
        merged_data = {}

        for stat_file in stat_files:
            sample_name = Path(stat_file).stem.replace('.chr', '')
            self.logger.debug(f"处理样本|Processing sample: {sample_name}")

            lines = self._read_gzipped_file(stat_file)
            if lines is None:
                continue

            sample_data = self._parse_chr_stat_file(lines, sample_name)

            # 合并到总数据|Merge to total data
            for chr_name, metrics in sample_data.items():
                if chr_name not in merged_data:
                    merged_data[chr_name] = metrics
                else:
                    merged_data[chr_name].update(metrics)

        # 转换为DataFrame并保存|Convert to DataFrame and save
        if not merged_data:
            self.logger.error("没有有效数据可合并|No valid data to merge")
            return False

        # 创建DataFrame|Create DataFrame
        df = pd.DataFrame.from_dict(merged_data, orient='index')

        # 重置索引，将染色体名作为一列|Reset index, make chromosome name a column
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'Chr'}, inplace=True)

        # 对列进行排序，先coverage后mean_depth，按样本名排序|Sort columns
        coverage_cols = [col for col in df.columns if col.endswith('_coverage')]
        depth_cols = [col for col in df.columns if col.endswith('_mean_depth')]

        coverage_cols_sorted = sorted(coverage_cols)
        depth_cols_sorted = sorted(depth_cols)

        # 重新排列列顺序|Reorder columns
        new_column_order = ['Chr'] + coverage_cols_sorted + depth_cols_sorted
        df = df[new_column_order]

        # 保存结果|Save results
        output_path = os.path.join(output_dir, output_file)
        df.to_csv(output_path, sep='\t', index=False, float_format='%.2f')

        self.logger.info(f"染色体统计结果已保存|Chromosome statistics saved: {output_path}")
        self.logger.info(f"合并了 {len(merged_data)} 个染色体，{len(stat_files)} 个样本|Merged {len(merged_data)} chromosomes, {len(stat_files)} samples")

        return True

    def merge_gene_statistics(self, output_dir: str, output_file: str = 'merged_gene_statistics.txt') -> bool:
        """合并基因统计结果|Merge gene statistics results

        Args:
            output_dir: 输出目录|Output directory
            output_file: 输出文件名|Output file name

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("开始合并基因统计结果|Starting gene statistics merging")

        # 查找所有基因统计文件|Find all gene statistics files
        stat_files = self._find_stat_files(output_dir, '.gene.stat.gz')

        if not stat_files:
            self.logger.warning("未找到基因统计文件|No gene statistics files found")
            return False

        self.logger.info(f"找到 {len(stat_files)} 个基因统计文件|Found {len(stat_files)} gene statistics files")

        # 解析所有文件|Parse all files
        merged_data = {}

        for stat_file in stat_files:
            sample_name = Path(stat_file).stem.replace('.gene', '')
            self.logger.debug(f"处理样本|Processing sample: {sample_name}")

            lines = self._read_gzipped_file(stat_file)
            if lines is None:
                continue

            sample_data = self._parse_gene_stat_file(lines, sample_name)

            # 合并到总数据|Merge to total data
            for gene_id, metrics in sample_data.items():
                if gene_id not in merged_data:
                    merged_data[gene_id] = metrics
                else:
                    merged_data[gene_id].update(metrics)

        # 转换为DataFrame并保存|Convert to DataFrame and save
        if not merged_data:
            self.logger.error("没有有效数据可合并|No valid data to merge")
            return False

        # 创建DataFrame|Create DataFrame
        df = pd.DataFrame.from_dict(merged_data, orient='index')

        # 重置索引，将基因ID作为一列|Reset index, make gene ID a column
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'GeneID'}, inplace=True)

        # 对列进行排序|Sort columns
        coverage_cols = [col for col in df.columns if col.endswith('_coverage')]
        depth_cols = [col for col in df.columns if col.endswith('_mean_depth')]

        coverage_cols_sorted = sorted(coverage_cols)
        depth_cols_sorted = sorted(depth_cols)

        # 重新排列列顺序|Reorder columns
        new_column_order = ['GeneID'] + coverage_cols_sorted + depth_cols_sorted
        df = df[new_column_order]

        # 保存结果|Save results
        output_path = os.path.join(output_dir, output_file)
        df.to_csv(output_path, sep='\t', index=False, float_format='%.2f')

        self.logger.info(f"基因统计结果已保存|Gene statistics saved: {output_path}")
        self.logger.info(f"合并了 {len(merged_data)} 个基因，{len(stat_files)} 个样本|Merged {len(merged_data)} genes, {len(stat_files)} samples")

        return True

    def merge_results(self, output_dir: str, has_gff: bool = False) -> bool:
        """自动合并统计结果|Automatically merge statistics results

        Args:
            output_dir: 输出目录|Output directory
            has_gff: 是否使用了GFF文件|Whether GFF file was used

        Returns:
            bool: 是否成功|Whether successful
        """
        if has_gff:
            # 基因覆盖度分析|Gene coverage analysis
            return self.merge_gene_statistics(output_dir)
        else:
            # 染色体覆盖度分析|Chromosome coverage analysis
            return self.merge_chr_statistics(output_dir)
