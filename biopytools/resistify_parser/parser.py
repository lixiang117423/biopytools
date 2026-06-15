"""
Resistify结果解析模块|Resistify Results Parser Module
"""

import pandas as pd
from typing import List, Dict
from pathlib import Path


class ResistifyParser:
    """Resistify结果解析器|Resistify Results Parser"""

    def __init__(self, logger, config):
        self.logger = logger
        self.config = config

    def parse_results(self) -> pd.DataFrame:
        """
        解析results.tsv主文件|Parse results.tsv main file

        Returns:
            pd.DataFrame: 解析结果|Parsed results
        """
        self.logger.info("解析results.tsv|Parsing results.tsv")

        results_file = Path(self.config.input_dir) / 'results.tsv'
        df = pd.read_csv(results_file, sep='\t')

        self.logger.info(f"解析完成，共{len(df)}条NLR基因|Parsing completed, {len(df)} NLR genes")
        return df

    def parse_domains(self) -> pd.DataFrame:
        """
        解析domains.tsv文件|Parse domains.tsv file

        Returns:
            pd.DataFrame: Domain信息|Domain information
        """
        self.logger.info("解析domains.tsv|Parsing domains.tsv")

        domains_file = Path(self.config.input_dir) / 'domains.tsv'
        df = pd.read_csv(domains_file, sep='\t')

        self.logger.info(f"解析完成，共{len(df)}条domain记录|Parsing completed, {len(df)} domain records")
        return df

    def parse_annotations(self) -> pd.DataFrame:
        """
        解析annotations.tsv文件|Parse annotations.tsv file

        Returns:
            pd.DataFrame: 注释信息|Annotation information
        """
        self.logger.info("解析annotations.tsv|Parsing annotations.tsv")

        annotations_file = Path(self.config.input_dir) / 'annotations.tsv'
        df = pd.read_csv(annotations_file, sep='\t')

        self.logger.info(f"解析完成，共{len(df)}条注释记录|Parsing completed, {len(df)} annotation records")
        return df

    def merge_data(self, results_df: pd.DataFrame, domains_df: pd.DataFrame = None) -> pd.DataFrame:
        """
        合并results和domains数据|Merge results and domains data

        Args:
            results_df: results.tsv的DataFrame|DataFrame from results.tsv
            domains_df: domains.tsv的DataFrame|DataFrame from domains.tsv

        Returns:
            pd.DataFrame: 合并后的数据|Merged data
        """
        self.logger.info("整合数据|Integrating data")

        # 创建domain特征列|Create domain feature columns
        results_df['has_TIR'] = False
        results_df['has_NB_ARC'] = False
        results_df['has_CC'] = False
        results_df['has_LRR'] = False
        results_df['has_RPW8'] = False

        if domains_df is not None and len(domains_df) > 0:
            # 为每个gene标记domain存在|Mark domain existence for each gene
            for seq_id in results_df['Sequence'].unique():
                seq_domains = domains_df[domains_df['Sequence'] == seq_id]['Domain'].values

                if seq_id in results_df['Sequence'].values:
                    idx = results_df[results_df['Sequence'] == seq_id].index[0]
                    results_df.at[idx, 'has_TIR'] = 'TIR' in seq_domains
                    results_df.at[idx, 'has_NB_ARC'] = 'NB-ARC' in seq_domains
                    results_df.at[idx, 'has_CC'] = 'CC' in seq_domains
                    results_df.at[idx, 'has_LRR'] = 'LRR' in seq_domains
                    results_df.at[idx, 'has_RPW8'] = 'RPW8' in seq_domains

        self.logger.info("数据整合完成|Data integration completed")
        return results_df

    def filter_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        过滤数据|Filter data

        Args:
            df: 输入数据框|Input DataFrame

        Returns:
            pd.DataFrame: 过滤后的数据框|Filtered DataFrame
        """
        filtered_df = df.copy()
        original_count = len(filtered_df)

        # 按分类筛选|Filter by classification
        if self.config.filter_classification:
            pattern = self.config.filter_classification
            filtered_df = filtered_df[filtered_df['Classification'].str.contains(pattern, na=False)]
            self.logger.info(f"分类筛选|Classification filter ({pattern}): {original_count} -> {len(filtered_df)}")

        # 按长度筛选|Filter by length
        if self.config.min_length:
            filtered_df = filtered_df[filtered_df['Length'] >= self.config.min_length]
            self.logger.info(f"最小长度筛选|Min length filter ({self.config.min_length}): {len(filtered_df)}")

        if self.config.max_length:
            filtered_df = filtered_df[filtered_df['Length'] <= self.config.max_length]
            self.logger.info(f"最大长度筛选|Max length filter ({self.config.max_length}): {len(filtered_df)}")

        if self.config.min_lrr_length:
            filtered_df = filtered_df[filtered_df['LRR_Length'] >= self.config.min_lrr_length]
            self.logger.info(f"最小LRR长度筛选|Min LRR length filter ({self.config.min_lrr_length}): {len(filtered_df)}")

        self.logger.info(f"过滤完成|Filtering completed: {original_count} -> {len(filtered_df)}")
        return filtered_df

    def prepare_summary_table(self, results_df: pd.DataFrame) -> pd.DataFrame:
        """
        准备汇总表格|Prepare summary table

        Args:
            results_df: results.tsv的DataFrame|DataFrame from results.tsv

        Returns:
            pd.DataFrame: 汇总表格|Summary table
        """
        self.logger.info("准备汇总表格|Preparing summary table")

        # 选择主要列|Select main columns
        summary_df = results_df[[
            'Sequence', 'Length', 'LRR_Length', 'Domains', 'Classification',
            'NBARC_motifs', 'has_TIR', 'has_NB_ARC', 'has_CC', 'has_LRR'
        ]].copy()

        return summary_df
