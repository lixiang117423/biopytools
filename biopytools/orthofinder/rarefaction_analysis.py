"""
稀释曲线分析模块|Rarefaction Curve Analysis Module
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
from itertools import combinations
import random


class RarefactionAnalyzer:
    """稀释曲线分析器|Rarefaction Curve Analyzer"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.classification_data = None
        self.genome_names = []
        self.total_genomes = 0

    def load_classification_results(self, classification_file: Path):
        """加载分类结果数据|Load classification results data"""
        self.logger.info(f"加载分类结果用于稀释分析|Loading classification results for rarefaction analysis: {classification_file}")

        try:
            self.classification_data = pd.read_csv(classification_file, sep='\t')

            self.genome_names = list(self.classification_data['Genome_Name'].unique())
            self.total_genomes = len(self.genome_names)

            self.logger.info(f"稀释分析数据已加载|Rarefaction analysis data loaded: {self.total_genomes} genomes")
            self.logger.info(f"基因组列表|Genome list: {', '.join(self.genome_names[:5])}{'...' if self.total_genomes > 5 else ''}")

        except Exception as e:
            self.logger.error(f"加载分类结果失败|Failed to load classification results: {e}")
            raise

    def run_rarefaction_analysis(self) -> Dict[str, List]:
        """运行稀释曲线分析|Run rarefaction curve analysis"""
        self.logger.info("开始稀释曲线分析|Starting rarefaction curve analysis")
        self.logger.info(f"迭代次数|Iterations: {self.config.rarefaction_iterations}")

        rarefaction_results = {
            'sample_sizes': [],
            'iterations': [],
            'pan_counts': [],
            'core_counts': []
        }

        for sample_size in range(1, self.total_genomes + 1):
            self.logger.info(f"分析样本大小|Analyzing sample size: {sample_size}/{self.total_genomes}")

            if sample_size == self.total_genomes:
                genome_combinations = [self.genome_names]
            else:
                from math import comb
                total_combinations = comb(self.total_genomes, sample_size)

                if total_combinations <= self.config.rarefaction_iterations:
                    genome_combinations = list(combinations(self.genome_names, sample_size))
                else:
                    genome_combinations = []
                    for _ in range(self.config.rarefaction_iterations):
                        selected_genomes = random.sample(self.genome_names, sample_size)
                        genome_combinations.append(selected_genomes)

            for iteration, selected_genomes in enumerate(genome_combinations, 1):
                pan_count, core_count = self._calculate_pan_core_counts(selected_genomes)

                rarefaction_results['sample_sizes'].append(sample_size)
                rarefaction_results['iterations'].append(iteration)
                rarefaction_results['pan_counts'].append(pan_count)
                rarefaction_results['core_counts'].append(core_count)

            current_pan = [rarefaction_results['pan_counts'][i] for i in range(len(rarefaction_results['sample_sizes']))
                          if rarefaction_results['sample_sizes'][i] == sample_size]
            current_core = [rarefaction_results['core_counts'][i] for i in range(len(rarefaction_results['sample_sizes']))
                           if rarefaction_results['sample_sizes'][i] == sample_size]

            self.logger.info(f"  样本大小 {sample_size}: Pan均值={np.mean(current_pan):.1f}, Core均值={np.mean(current_core):.1f}")

        self.logger.info("稀释曲线分析完成|Rarefaction curve analysis completed")
        return rarefaction_results

    def _calculate_pan_core_counts(self, selected_genomes: List[str]) -> Tuple[int, int]:
        """计算选定基因组的pan和core同源群数量|Calculate pan and core orthogroup counts for selected genomes"""
        selected_data = self.classification_data[self.classification_data['Genome_Name'].isin(selected_genomes)]
        og_genome_counts = selected_data.groupby('Orthogroup_ID')['Genome_Name'].nunique()

        pan_count = len(og_genome_counts)
        core_count = sum(1 for count in og_genome_counts if count == len(selected_genomes))

        return pan_count, core_count

    def save_rarefaction_results(self, rarefaction_results: Dict, output_dir: Path) -> Path:
        """保存稀释分析结果|Save rarefaction analysis results"""
        self.logger.info("保存稀释曲线分析结果|Saving rarefaction curve analysis results")

        rarefaction_dir = Path(output_dir) / "03_rarefaction_analysis"
        rarefaction_dir.mkdir(parents=True, exist_ok=True)

        df_detailed = pd.DataFrame({
            'Sample_Size': rarefaction_results['sample_sizes'],
            'Iteration': rarefaction_results['iterations'],
            'Core_Count': rarefaction_results['core_counts'],
            'Pan_Count': rarefaction_results['pan_counts']
        })

        detailed_file = rarefaction_dir / "rarefaction_curve_detailed.tsv"
        df_detailed.to_csv(detailed_file, sep='\t', index=False)

        summary_data = []
        for sample_size in range(1, self.total_genomes + 1):
            size_data = df_detailed[df_detailed['Sample_Size'] == sample_size]

            summary_data.append({
                'Sample_Size': sample_size,
                'Pan_Mean': size_data['Pan_Count'].mean(),
                'Pan_Std': size_data['Pan_Count'].std(),
                'Core_Mean': size_data['Core_Count'].mean(),
                'Core_Std': size_data['Core_Count'].std(),
                'Iterations': len(size_data)
            })

        df_summary = pd.DataFrame(summary_data)

        summary_file = rarefaction_dir / "rarefaction_curve_summary.tsv"
        df_summary.to_csv(summary_file, sep='\t', index=False, float_format='%.2f')

        self.logger.info(f"稀释曲线详细数据已保存|Rarefaction curve detailed data saved: {detailed_file}")
        self.logger.info(f"稀释曲线摘要数据已保存|Rarefaction curve summary data saved: {summary_file}")

        return summary_file
