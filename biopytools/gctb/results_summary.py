"""
GCTB结果汇总模块|GCTB Results Summary Module
"""

import os
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional
import logging


class ResultsSummarizer:
    """GCTB结果汇总器|GCTB Results Summarizer"""

    def __init__(self, logger: logging.Logger, output_dir: Path):
        self.logger = logger
        self.output_dir = Path(output_dir)

    def summarize_batch_results(self,
                                trait_names: List[str],
                                analysis_dirs: Dict[str, Path],
                                bayes_type: str = "S") -> bool:
        """
        汇总批量分析结果|Summarize batch analysis results

        Args:
            trait_names: 表型名称列表|List of trait names
            analysis_dirs: 表型到分析目录的映射|Mapping from trait to analysis directory
            bayes_type: 贝叶斯模型类型|Bayesian model type

        Returns:
            bool: 汇总是否成功|Whether summarization succeeded
        """
        self.logger.info("="*60)
        self.logger.info("开始汇总分析结果|Starting results summarization")
        self.logger.info("="*60)

        # 创建汇总目录|Create summary directory
        summary_dir = self.output_dir / "summary"
        summary_dir.mkdir(parents=True, exist_ok=True)

        try:
            # 1. 汇总参数表|Summary parameter table
            param_df = self._summarize_parameters(trait_names, analysis_dirs, bayes_type)
            param_file = summary_dir / "all_traits_parameters.tsv"
            param_df.to_csv(param_file, sep='\t', index=False, float_format='%.6f')
            self.logger.info(f"✓ 参数汇总表|Parameter summary: {param_file}")

            # 2. 提取top SNPs|Extract top SNPs
            top_snps_df = self._extract_top_snps(trait_names, analysis_dirs, bayes_type, top_n=50)
            top_snps_file = summary_dir / "all_traits_top_snps.tsv"
            top_snps_df.to_csv(top_snps_file, sep='\t', index=False, float_format='%.6f')
            self.logger.info(f"✓ Top SNPs表|Top SNPs table: {top_snps_file}")

            # 3. 合并所有SNPs（可选，可能很大）|Merge all SNPs (optional, may be large)
            all_snps_df = self._merge_all_snps(trait_names, analysis_dirs, bayes_type)
            all_snps_file = summary_dir / "all_traits_all_snps.tsv"
            all_snps_df.to_csv(all_snps_file, sep='\t', index=False, float_format='%.6f')
            self.logger.info(f"✓ 完整SNPs表|All SNPs table: {all_snps_file} ({len(all_snps_df)} rows)")

            # 4. 生成Markdown报告|Generate Markdown report
            report_file = summary_dir / "analysis_report.md"
            self._generate_markdown_report(param_df, top_snps_df, report_file, bayes_type)
            self.logger.info(f"✓ Markdown报告|Markdown report: {report_file}")

            self.logger.info("="*60)
            self.logger.info("结果汇总完成|Results summarization completed")
            self.logger.info("="*60)

            return True

        except Exception as e:
            self.logger.error(f"结果汇总失败|Results summarization failed: {e}")
            return False

    def _summarize_parameters(self,
                             trait_names: List[str],
                             analysis_dirs: Dict[str, Path],
                             bayes_type: str) -> pd.DataFrame:
        """汇总参数表|Summarize parameters"""
        results = []

        for trait in trait_names:
            trait_dir = analysis_dirs.get(trait)
            if not trait_dir:
                self.logger.warning(f"未找到表型目录|Trait directory not found: {trait}")
                continue

            par_file = trait_dir / f"bayes{bayes_type.lower()}.parRes"
            if not par_file.exists():
                self.logger.warning(f"参数文件不存在|Parameter file not found: {par_file}")
                continue

            try:
                # 读取参数文件|Read parameter file
                par_df = pd.read_csv(par_file, sep='\s+', header=0)

                # 提取关键参数|Extract key parameters
                param_dict = {'Trait': trait}

                for _, row in par_df.iterrows():
                    param_name = row['Parameter']
                    param_dict[f"{param_name}_mean"] = row['Mean']
                    if 'SD' in row:
                        param_dict[f"{param_name}_sd"] = row['SD']

                results.append(param_dict)

            except Exception as e:
                self.logger.error(f"读取参数文件失败|Failed to read parameter file for {trait}: {e}")

        return pd.DataFrame(results)

    def _extract_top_snps(self,
                         trait_names: List[str],
                         analysis_dirs: Dict[str, Path],
                         bayes_type: str,
                         top_n: int = 50) -> pd.DataFrame:
        """提取每个表型的top SNPs|Extract top SNPs for each trait"""
        all_top_snps = []

        for trait in trait_names:
            trait_dir = analysis_dirs.get(trait)
            if not trait_dir:
                continue

            snp_file = trait_dir / f"bayes{bayes_type.lower()}.snpRes"
            if not snp_file.exists():
                self.logger.warning(f"SNP结果文件不存在|SNP result file not found: {snp_file}")
                continue

            try:
                # 读取SNP结果|Read SNP results
                snp_df = pd.read_csv(snp_file, sep='\s+')

                # 按PIP排序，取top N|Sort by PIP and take top N
                top_snps = snp_df.nlargest(top_n, 'PIP').copy()
                top_snps.insert(0, 'Trait', trait)

                all_top_snps.append(top_snps)

            except Exception as e:
                self.logger.error(f"读取SNP文件失败|Failed to read SNP file for {trait}: {e}")

        if all_top_snps:
            return pd.concat(all_top_snps, ignore_index=True)
        return pd.DataFrame()

    def _merge_all_snps(self,
                       trait_names: List[str],
                       analysis_dirs: Dict[str, Path],
                       bayes_type: str) -> pd.DataFrame:
        """合并所有表型的所有SNPs|Merge all SNPs from all traits"""
        all_snps = []

        for trait in trait_names:
            trait_dir = analysis_dirs.get(trait)
            if not trait_dir:
                continue

            snp_file = trait_dir / f"bayes{bayes_type.lower()}.snpRes"
            if not snp_file.exists():
                continue

            try:
                # 读取SNP结果|Read SNP results
                snp_df = pd.read_csv(snp_file, sep='\s+')
                snp_df.insert(0, 'Trait', trait)

                all_snps.append(snp_df)

            except Exception as e:
                self.logger.error(f"读取SNP文件失败|Failed to read SNP file for {trait}: {e}")

        if all_snps:
            return pd.concat(all_snps, ignore_index=True)
        return pd.DataFrame()

    def _generate_markdown_report(self,
                                  param_df: pd.DataFrame,
                                  top_snps_df: pd.DataFrame,
                                  report_file: Path,
                                  bayes_type: str):
        """生成Markdown报告|Generate Markdown report"""
        with open(report_file, 'w', encoding='utf-8') as f:
            # 标题|Title
            f.write("# GCTB批量分析结果汇总|GCTB Batch Analysis Results Summary\n\n")
            f.write(f"**模型|Model**: Bayes{bayes_type}\n\n")
            f.write("---\n\n")

            # 1. 参数汇总表|Parameter summary table
            f.write("## 1. 表型参数汇总|Phenotype Parameters Summary\n\n")
            f.write("| 表型|Trait | 遗传力|h2 | 有效SNP数|NnzSnp | 遗传方差|GenVar | 残差方差|ResVar |\n")
            f.write("|---------|---------|---------------|---------------|---------------|\n")

            for _, row in param_df.iterrows():
                trait = row['Trait']
                h2 = row.get('hsq_mean', 'N/A')
                nnz = row.get('NnzSnp_mean', 'N/A')
                gen_var = row.get('GenVar_mean', 'N/A')
                res_var = row.get('ResVar_mean', 'N/A')

                f.write(f"| {trait} | {h2:.4f} | {nnz:.0f} | {gen_var:.4f} | {res_var:.4f} |\n")

            f.write("\n---\n\n")

            # 2. 每个表型的top SNPs|Top SNPs for each trait
            f.write("## 2. 各表型Top候选SNPs|Top Candidate SNPs by Trait\n\n")

            for trait in param_df['Trait']:
                trait_snps = top_snps_df[top_snps_df['Trait'] == trait].head(10)

                if len(trait_snps) > 0:
                    f.write(f"### {trait}\n\n")
                    f.write("| 排名|Rank | SNP名称|Name | 染色体|Chr | 位置|Pos | 效应|Effect | PIP |\n")
                    f.write("|-------|-------|-------|-------|-------|\n")

                    for idx, (_, row) in enumerate(trait_snps.iterrows(), 1):
                        snp_name = row['Name']
                        chrom = row['Chrom']
                        pos = row['Position']
                        effect = row['A1Effect']
                        pip = row['PIP']

                        f.write(f"| {idx} | {snp_name} | {chrom} | {pos:.0f} | {effect:.6f} | {pip:.6f} |\n")

                    f.write("\n")

            f.write("---\n\n")
            f.write(f"*报告生成时间|Report generated: {pd.Timestamp.now()}*\n")

        self.logger.info(f"Markdown报告已生成|Markdown report generated: {report_file}")
