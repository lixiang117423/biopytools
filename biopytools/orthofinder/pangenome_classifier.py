"""
泛基因组分类器模块|Pangenome Classifier Module
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Set
from collections import defaultdict


class PangenomeClassifier:
    """泛基因组分类器|Pangenome Classifier"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.orthogroups_data = None
        self.gene_count_data = None
        self.genome_names = []
        self.total_genomes = 0

    def load_orthofinder_results(self, orthogroups_file: Path, gene_count_file: Path):
        """加载OrthoFinder结果|Load OrthoFinder results"""
        self.logger.info("加载OrthoFinder结果文件|Loading OrthoFinder results files")

        try:
            self.gene_count_data = pd.read_csv(gene_count_file, sep='\t', index_col=0)
            if 'Total' in self.gene_count_data.columns:
                self.gene_count_data = self.gene_count_data.drop('Total', axis=1)

            self.genome_names = list(self.gene_count_data.columns)
            self.total_genomes = len(self.genome_names)
            self.logger.info(f"成功加载基因计数文件|Successfully loaded gene count file: {self.total_genomes} genomes")
            self.logger.info(f"基因组名称|Genome names: {', '.join(self.genome_names[:5])}{'...' if self.total_genomes > 5 else ''}")
        except Exception as e:
            self.logger.error(f"加载基因计数文件失败|Failed to load gene count file: {e}")
            raise

        try:
            orthogroups_tsv_file = orthogroups_file.parent / "Orthogroups.tsv"
            if not orthogroups_tsv_file.exists():
                raise FileNotFoundError(f"未找到Orthogroups.tsv文件|Orthogroups.tsv file not found: {orthogroups_tsv_file}")

            self.orthogroups_data = pd.read_csv(orthogroups_tsv_file, sep='\t', dtype=str)
            self.logger.info(f"成功加载同源基因群详细文件|Successfully loaded orthogroups detailed file: {len(self.orthogroups_data)} groups")

            orthogroups_columns = list(self.orthogroups_data.columns[1:])
            if 'Total' in orthogroups_columns:
                orthogroups_columns.remove('Total')

            self.logger.info(f"Orthogroups.tsv基因组列|Orthogroups.tsv genome columns: {orthogroups_columns[:5]}{'...' if len(orthogroups_columns) > 5 else ''}")

            if set(self.genome_names) != set(orthogroups_columns):
                self.logger.warning("两个文件的基因组名称不完全一致|Genome names in two files are not identical")
                common_genomes = list(set(self.genome_names) & set(orthogroups_columns))
                self.genome_names = common_genomes
                self.total_genomes = len(common_genomes)
                self.logger.info(f"使用共同基因组|Using common genomes: {self.total_genomes} genomes")

        except Exception as e:
            self.logger.error(f"加载同源基因群详细文件失败|Failed to load orthogroups detailed file: {e}")
            raise

    def classify_pangenome(self) -> Dict[str, List[Tuple[str, List[str], List[str]]]]:
        """分类泛基因组|Classify pangenome"""
        self.logger.info("开始泛基因组分类|Starting pangenome classification")

        self.logger.info(f"分类参数|Classification parameters:")
        self.logger.info(f"  Softcore缺失阈值|Softcore missing threshold: <={self.config.softcore_missing_threshold}")
        self.logger.info(f"  Dispensable缺失阈值|Dispensable missing threshold: >{self.config.dispensable_missing_threshold}")

        classification_results = {
            'core': [],
            'softcore': [],
            'dispensable': [],
            'private': [],
            'single_copy': []
        }

        for og_id, row in self.gene_count_data.iterrows():
            presence_count = sum(1 for count in row if count > 0)
            missing_count = self.total_genomes - presence_count

            gene_ids, genome_names = self._get_gene_details(og_id)
            gene_info = (og_id, gene_ids, genome_names)

            if missing_count == 0:
                classification_results['core'].append(gene_info)
                if self._is_single_copy_orthogroup(og_id):
                    classification_results['single_copy'].append(gene_info)
            elif missing_count <= self.config.softcore_missing_threshold:
                classification_results['softcore'].append(gene_info)
            elif presence_count == 1:
                classification_results['private'].append(gene_info)
            else:
                classification_results['dispensable'].append(gene_info)

        self.logger.info("泛基因组分类统计|Pangenome classification statistics:")
        main_categories = ['core', 'softcore', 'dispensable', 'private']
        total_orthogroups = sum(len(classification_results[c]) for c in main_categories)

        for category in main_categories:
            results = classification_results[category]
            og_count = len(results)
            gene_count = sum(len(gene_info[1]) for gene_info in results)
            percentage = (og_count / total_orthogroups * 100) if total_orthogroups > 0 else 0

            category_name = {
                'core': '核心基因',
                'softcore': '软核心基因',
                'dispensable': '非必需基因',
                'private': '私有基因'
            }[category]

            self.logger.info(f"  {category_name}|{category.capitalize()}: {og_count} orthogroups ({percentage:.2f}%), {gene_count} genes")

        sc_count = len(classification_results['single_copy'])
        self.logger.info(f"  单拷贝基因|Single copy: {sc_count} orthogroups (核心基因的子集|subset of core)")

        return classification_results

    def _get_gene_details(self, og_id: str) -> Tuple[List[str], List[str]]:
        """从Orthogroups.tsv获取具体基因信息|Get gene details from Orthogroups.tsv"""
        try:
            og_row = self.orthogroups_data[self.orthogroups_data.iloc[:, 0] == og_id]

            if og_row.empty:
                return [], []

            row = og_row.iloc[0]

            all_genes = []
            all_genomes = []

            for genome_name in self.genome_names:
                if genome_name in self.orthogroups_data.columns:
                    col_idx = self.orthogroups_data.columns.get_loc(genome_name)
                    cell_value = row.iloc[col_idx]

                    if pd.notna(cell_value) and str(cell_value).strip() and str(cell_value) != 'nan':
                        gene_string = str(cell_value).strip()
                        gene_list = [g.strip() for g in gene_string.split(',') if g.strip()]
                        for gene_id in gene_list:
                            all_genes.append(gene_id)
                            all_genomes.append(genome_name)

            return list(all_genes), list(all_genomes)

        except Exception as e:
            self.logger.error(f"处理同源群 {og_id} 时出错|Error processing orthogroup {og_id}: {e}")
            return [], []

    def save_classification_results(self, classification_results: Dict, output_dir: Path):
        """保存分类结果|Save classification results"""
        self.logger.info("保存泛基因组分类结果|Saving pangenome classification results")

        classification_dir = Path(output_dir) / "02_pangenome_classification"
        classification_dir.mkdir(parents=True, exist_ok=True)

        output_file = classification_dir / "pangenome_gene_families.txt"

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("Category\tOrthogroup_ID\tGene_ID\tGenome_Name\n")

            for category in ['core', 'softcore', 'dispensable', 'private']:
                for og_id, gene_ids, genome_names in classification_results[category]:
                    for gene_id, genome_name in zip(gene_ids, genome_names):
                        f.write(f"{category}\t{og_id}\t{gene_id}\t{genome_name}\n")

        self.logger.info(f"分类结果已保存|Classification results saved: {output_file}")

        self._save_classification_summary(classification_results, classification_dir)

        return output_file

    def _save_classification_summary(self, classification_results: Dict, output_dir: Path):
        """保存分类统计摘要|Save classification summary"""
        summary_file = output_dir / "pangenome_classification_summary.txt"

        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("泛基因组分类统计摘要|Pangenome Classification Summary\n")
            f.write("=" * 60 + "\n\n")

            f.write(f"总基因组数量|Total genomes: {self.total_genomes}\n")
            f.write(f"Softcore缺失阈值|Softcore missing threshold: <={self.config.softcore_missing_threshold}\n")
            f.write(f"Dispensable缺失阈值|Dispensable missing threshold: >{self.config.dispensable_missing_threshold}\n\n")

            f.write("分类定义|Classification Definitions:\n")
            f.write("-" * 40 + "\n")
            f.write("Core: 存在于所有基因组|Present in all genomes\n")
            f.write(f"Softcore: 最多缺失{self.config.softcore_missing_threshold}个基因组|Missing in <={self.config.softcore_missing_threshold} genomes\n")
            f.write(f"Dispensable: 缺失超过{self.config.dispensable_missing_threshold}个基因组|Missing in >{self.config.dispensable_missing_threshold} genomes\n")
            f.write("Private: 仅存在于单个基因组|Present in only one genome\n\n")

            f.write("分类结果|Classification Results:\n")
            f.write("-" * 40 + "\n")

            total_orthogroups = sum(len(classification_results[c]) for c in ['core', 'softcore', 'dispensable', 'private'])
            total_genes = sum(sum(len(gene_info[1]) for gene_info in classification_results[c])
                            for c in ['core', 'softcore', 'dispensable', 'private'])

            for category in ['core', 'softcore', 'dispensable', 'private']:
                og_count = len(classification_results[category])
                gene_count = sum(len(gene_info[1]) for gene_info in classification_results[category])

                og_percentage = (og_count / total_orthogroups * 100) if total_orthogroups > 0 else 0
                gene_percentage = (gene_count / total_genes * 100) if total_genes > 0 else 0

                category_name = {
                    'core': '核心基因|Core genes',
                    'softcore': '软核心基因|Softcore genes',
                    'dispensable': '非必需基因|Dispensable genes',
                    'private': '私有基因|Private genes'
                }[category]

                f.write(f"{category_name}:\n")
                f.write(f"  同源基因群数量|Orthogroups: {og_count:,} ({og_percentage:.2f}%)\n")
                f.write(f"  基因数量|Genes: {gene_count:,} ({gene_percentage:.2f}%)\n\n")

            f.write(f"总计|Total:\n")
            f.write(f"  同源基因群总数|Total orthogroups: {total_orthogroups:,}\n")
            f.write(f"  基因总数|Total genes: {total_genes:,}\n")

        self.logger.info(f"分类摘要已保存|Classification summary saved: {summary_file}")

    def get_genome_gene_counts(self, classification_results: Dict) -> Dict[str, Dict[str, int]]:
        """获取每个基因组的基因数量统计|Get gene count statistics for each genome"""
        genome_stats = {}

        for genome_name in self.genome_names:
            genome_stats[genome_name] = {
                'core': 0,
                'softcore': 0,
                'dispensable': 0,
                'private': 0,
                'single_copy': 0,
                'total': 0
            }

        for category in ['core', 'softcore', 'dispensable', 'private']:
            for gene_info in classification_results[category]:
                try:
                    if len(gene_info) != 3:
                        continue

                    og_id, gene_ids, genome_names = gene_info

                    if not isinstance(gene_ids, (list, tuple)):
                        continue
                    if not isinstance(genome_names, (list, tuple)):
                        continue
                    if len(gene_ids) != len(genome_names):
                        continue

                    is_sc = og_id in {g[0] for g in classification_results['single_copy']}

                    for gene_id, genome_name in zip(gene_ids, genome_names):
                        if genome_name in genome_stats:
                            genome_stats[genome_name][category] += 1
                            genome_stats[genome_name]['total'] += 1
                            if is_sc:
                                genome_stats[genome_name]['single_copy'] += 1

                except Exception:
                    continue

        return genome_stats

    def calculate_frequency_distribution(self, classification_results: Dict) -> Dict[str, Dict[int, int]]:
        """计算每个分类中同源群的出现频率分布|Calculate frequency distribution of orthogroups in each category"""
        frequency_distributions = {}

        for category in ['core', 'softcore', 'dispensable', 'private']:
            frequency_dist = defaultdict(int)

            for og_id, gene_ids, genome_names in classification_results[category]:
                unique_genomes = set(genome_names)
                frequency = len(unique_genomes)
                frequency_dist[frequency] += 1

            frequency_distributions[category] = dict(frequency_dist)

        return frequency_distributions

    def _is_single_copy_orthogroup(self, og_id: str) -> bool:
        """检查是否为单拷贝同源群|Check if single copy orthogroup"""
        if og_id not in self.gene_count_data.index:
            return False

        row = self.gene_count_data.loc[og_id]
        return all(count == 1 for count in row)
