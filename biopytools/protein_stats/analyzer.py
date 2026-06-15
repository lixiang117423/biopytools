"""
蛋白质性质分析模块|Protein Properties Analysis Module
"""

from Bio import SeqIO
from Bio.SeqUtils import ProtParam
from typing import List, Dict, Any
import pandas as pd


class ProteinStatsAnalyzer:
    """蛋白质性质统计器|Protein Statistics Analyzer"""

    def __init__(self, logger, config):
        self.logger = logger
        self.config = config

    def analyze(self) -> pd.DataFrame:
        """
        分析蛋白质序列|Analyze protein sequences

        Returns:
            pd.DataFrame: 分析结果|Analysis results
        """
        self.logger.info(f"开始分析蛋白质序列|Analyzing protein sequences: {self.config.protein_fasta}")

        results = []
        count = 0

        ambiguous_aa = set('BJOUXZbjouxz')

        for record in SeqIO.parse(self.config.protein_fasta, "fasta"):
            count += 1
            result = {'id': record.id}
            seq = str(record.seq)

            # 清理模糊氨基酸|Remove ambiguous amino acids for ProtParam
            clean_seq = ''.join(c for c in seq if c not in ambiguous_aa)

            # 基本性质|Basic properties
            if self.config.calculate_length:
                result['length'] = len(seq)

            # 分子量和等电点|Molecular weight and pI
            if self.config.calculate_mw or self.config.calculate_pi:
                analyser = ProtParam.ProteinAnalysis(clean_seq)

                if self.config.calculate_mw:
                    result['mw_da'] = round(analyser.molecular_weight(), 2)

                if self.config.calculate_pi:
                    result['pi'] = round(analyser.isoelectric_point(), 2)

                    # 判断酸碱性|Determine acid/base property
                    pI = result['pi']
                    if pI < 7.0:
                        result['acid_base'] = 'acidic'
                        result['酸碱性'] = '酸性'
                    elif pI > 7.0:
                        result['acid_base'] = 'basic'
                        result['酸碱性'] = '碱性'
                    else:
                        result['acid_base'] = 'neutral'
                        result['酸碱性'] = '中性'

            # 氨基酸组成|Amino acid composition
            if self.config.calculate_aa_composition:
                analyser = ProtParam.ProteinAnalysis(clean_seq)
                aa_comp = analyser.get_amino_acids_percent()
                for aa, percent in aa_comp.items():
                    result[f'aa_{aa}'] = round(percent, 2)

            # 不稳定指数|Instability index
            if self.config.calculate_instability_index:
                analyser = ProtParam.ProteinAnalysis(clean_seq)
                result['instability_index'] = round(analyser.instability_index(), 2)

            # 脂肪指数|Gravy (hydropathy)
            if self.config.calculate_gravy:
                analyser = ProtParam.ProteinAnalysis(clean_seq)
                result['gravy'] = round(analyser.gravy(), 3)

            # 芳香性|Aromaticity
            if self.config.calculate_aromacity:
                analyser = ProtParam.ProteinAnalysis(clean_seq)
                result['aromaticity'] = round(analyser.aromaticity(), 3)

            results.append(result)

            if count % 1000 == 0:
                self.logger.debug(f"已处理|Processed: {count} 条序列|sequences")

        self.logger.info(f"分析完成，共{count}条序列|Analysis completed, {count} sequences")

        # 转换为DataFrame|Convert to DataFrame
        df = pd.DataFrame(results)

        return df

    def save_results(self, df: pd.DataFrame):
        """
        保存结果到文件|Save results to file

        Args:
            df: 结果DataFrame|Result DataFrame
        """
        self.logger.info(f"保存结果到|Saving results to: {self.config.output_file}")

        if self.config.output_format == 'tsv':
            df.to_csv(self.config.output_file, sep='\t', index=False)
        elif self.config.output_format == 'csv':
            df.to_csv(self.config.output_file, index=False)
        elif self.config.output_format == 'excel':
            df.to_excel(self.config.output_file, index=False, engine='openpyxl')

        self.logger.info(f"结果已保存|Results saved successfully")

    def print_summary(self, df: pd.DataFrame):
        """
        打印统计摘要|Print statistical summary

        Args:
            df: 结果DataFrame|Result DataFrame
        """
        self.logger.info("-" * 60)
        self.logger.info("统计摘要|Statistical Summary")
        self.logger.info("-" * 60)
        self.logger.info(f"总序列数|Total sequences: {len(df)}")

        if 'length' in df.columns:
            self.logger.info(f"平均长度|Average length: {df['length'].mean():.2f}")
            self.logger.info(f"长度范围|Length range: {df['length'].min()} - {df['length'].max()}")

        if 'mw_da' in df.columns:
            self.logger.info(f"平均分子量|Average MW: {df['mw_da'].mean():.2f} Da")

        if 'pi' in df.columns:
            self.logger.info(f"平均pI|Average pI: {df['pi'].mean():.2f}")

        if 'acid_base' in df.columns:
            acid_count = (df['acid_base'] == 'acidic').sum()
            basic_count = (df['acid_base'] == 'basic').sum()
            neutral_count = (df['acid_base'] == 'neutral').sum()
            self.logger.info(f"酸性蛋白|Acidic proteins (pI<7): {acid_count} ({acid_count/len(df)*100:.1f}%)")
            self.logger.info(f"碱性蛋白|Basic proteins (pI>7): {basic_count} ({basic_count/len(df)*100:.1f}%)")
            self.logger.info(f"中性蛋白|Neutral proteins (pI=7): {neutral_count} ({neutral_count/len(df)*100:.1f}%)")

        self.logger.info("-" * 60)
