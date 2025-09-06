"""
泛基因组分类器模块 | Pangenome Classifier Module
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Set
from collections import defaultdict

class PangenomeClassifier:
    """泛基因组分类器 | Pangenome Classifier"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.orthogroups_data = None
        self.gene_count_data = None
        self.genome_names = []
        self.total_genomes = 0
        
    def load_orthofinder_results(self, orthogroups_file: Path, gene_count_file: Path):
        """加载OrthoFinder结果 | Load OrthoFinder results"""
        self.logger.info("加载OrthoFinder结果文件 | Loading OrthoFinder results files")
        
        # 加载基因计数文件用于分类 | Load gene count file for classification
        try:
            self.gene_count_data = pd.read_csv(gene_count_file, sep='\t', index_col=0)
            # 去除Total列 | Remove Total column
            if 'Total' in self.gene_count_data.columns:
                self.gene_count_data = self.gene_count_data.drop('Total', axis=1)
            
            self.genome_names = list(self.gene_count_data.columns)
            self.total_genomes = len(self.genome_names)
            self.logger.info(f"成功加载基因计数文件 | Successfully loaded gene count file: {self.total_genomes} genomes")
            self.logger.info(f"基因组名称 | Genome names: {', '.join(self.genome_names[:5])}{'...' if self.total_genomes > 5 else ''}")
        except Exception as e:
            self.logger.error(f"加载基因计数文件失败 | Failed to load gene count file: {e}")
            raise
        
        # 加载同源基因群详细信息文件（tsv格式）| Load orthogroups detailed info file (tsv format)
        try:
            # 查找正确的tsv文件
            orthogroups_tsv_file = orthogroups_file.parent / "Orthogroups.tsv"
            if not orthogroups_tsv_file.exists():
                raise FileNotFoundError(f"未找到Orthogroups.tsv文件 | Orthogroups.tsv file not found: {orthogroups_tsv_file}")
            
            self.orthogroups_data = pd.read_csv(orthogroups_tsv_file, sep='\t', dtype=str)
            self.logger.info(f"成功加载同源基因群详细文件 | Successfully loaded orthogroups detailed file: {len(self.orthogroups_data)} groups")
            
            # 验证列名一致性
            orthogroups_columns = list(self.orthogroups_data.columns[1:])  # 排除第一列Orthogroup
            if 'Total' in orthogroups_columns:
                orthogroups_columns.remove('Total')
            
            self.logger.info(f"Orthogroups.tsv基因组列 | Orthogroups.tsv genome columns: {orthogroups_columns[:5]}{'...' if len(orthogroups_columns) > 5 else ''}")
            
            # 确保两个文件的基因组名称一致
            if set(self.genome_names) != set(orthogroups_columns):
                self.logger.warning("两个文件的基因组名称不完全一致 | Genome names in two files are not identical")
                # 使用交集
                common_genomes = list(set(self.genome_names) & set(orthogroups_columns))
                self.genome_names = common_genomes
                self.total_genomes = len(common_genomes)
                self.logger.info(f"使用共同基因组 | Using common genomes: {self.total_genomes} genomes")
            
        except Exception as e:
            self.logger.error(f"加载同源基因群详细文件失败 | Failed to load orthogroups detailed file: {e}")
            raise
    
    def classify_pangenome(self) -> Dict[str, List[Tuple[str, List[str], List[str]]]]:
        """分类泛基因组 | Classify pangenome"""
        self.logger.info("开始泛基因组分类 | Starting pangenome classification")
        
        # 输出分类参数 | Output classification parameters
        self.logger.info(f"分类参数 | Classification parameters:")
        self.logger.info(f"  Softcore缺失阈值 | Softcore missing threshold: <={self.config.softcore_missing_threshold}")
        self.logger.info(f"  Dispensable缺失阈值 | Dispensable missing threshold: >{self.config.dispensable_missing_threshold}")
        

        # 分类结果字典 | Classification results dictionary
        classification_results = {
            'core': [],          # 核心基因 | Core genes
            'softcore': [],      # 软核心基因 | Softcore genes  
            'dispensable': [],   # 非必需基因 | Dispensable genes
            'private': [],       # 私有基因 | Private genes
            'single_copy': []    # 单拷贝基因 | Single copy genes
        }
        
        # 基于基因计数数据进行分类 | Classify based on gene count data
        for og_id, row in self.gene_count_data.iterrows():
            # 计算存在该同源群的基因组数量 | Count genomes with this orthogroup
            presence_count = sum(1 for count in row if count > 0)
            missing_count = self.total_genomes - presence_count
            
            # 从Orthogroups.txt获取具体基因信息 | Get specific gene info from Orthogroups.txt
            gene_ids, genome_names = self._get_gene_details(og_id)
            
            # 分类逻辑 | Classification logic
            gene_info = (og_id, gene_ids, genome_names)
            
            if missing_count == 0:
                # 检查是否为单拷贝基因 | Check if single copy gene
                if self._is_single_copy_orthogroup(og_id):
                    classification_results['single_copy'].append(gene_info)
                else:
                    classification_results['core'].append(gene_info)
            elif missing_count <= self.config.softcore_missing_threshold:
                classification_results['softcore'].append(gene_info)
            elif presence_count == 1:
                classification_results['private'].append(gene_info)
            else:
                classification_results['dispensable'].append(gene_info)
        
        # 输出分类统计 | Output classification statistics
        self.logger.info("泛基因组分类统计 | Pangenome classification statistics:")
        total_orthogroups = sum(len(results) for results in classification_results.values())
        
        for category, results in classification_results.items():
            og_count = len(results)
            gene_count = sum(len(gene_info[1]) for gene_info in results)
            percentage = (og_count / total_orthogroups * 100) if total_orthogroups > 0 else 0
            
            category_name = {
                'core': '核心基因',
                'softcore': '软核心基因',
                'dispensable': '非必需基因',
                'private': '私有基因',
                'single_copy': '单拷贝基因'
            }[category]
            
            self.logger.info(f"   {category_name} | {category.capitalize()}: {og_count} orthogroups ({percentage:.2f}%), {gene_count} genes")
        
        return classification_results
        
    #     return all_genes, all_genomes
    # def _get_gene_details(self, og_id: str) -> Tuple[List[str], List[str]]:
    #     """从Orthogroups.tsv获取具体基因信息 | Get gene details from Orthogroups.tsv"""
    #     # 查找对应的同源群行 | Find corresponding orthogroup row
    #     og_row = self.orthogroups_data[self.orthogroups_data.iloc[:, 0] == og_id]
        
    #     if og_row.empty:
    #         return [], []
        
    #     row = og_row.iloc[0]
        
    #     # 收集所有基因ID和对应基因组 | Collect all gene IDs and corresponding genomes
    #     all_genes = []
    #     all_genomes = []
        
    #     # 遍历每个基因组列 | Iterate through each genome column
    #     for genome_name in self.genome_names:
    #         if genome_name in self.orthogroups_data.columns:
    #             col_idx = self.orthogroups_data.columns.get_loc(genome_name)
    #             cell_value = row.iloc[col_idx]
                
    #             if pd.notna(cell_value) and str(cell_value).strip():
    #                 gene_string = str(cell_value).strip()
    #                 if gene_string and gene_string != 'nan':
    #                     # 用逗号分割基因，并去除空格 | Split genes by comma and remove spaces
    #                     gene_list = [g.strip() for g in gene_string.split(',') if g.strip()]
    #                     for gene_id in gene_list:
    #                         all_genes.append(gene_id)
    #                         all_genomes.append(genome_name)
        
    #     return all_genes, all_genomes
    def _get_gene_details(self, og_id: str) -> Tuple[List[str], List[str]]:
        """从Orthogroups.tsv获取具体基因信息 | Get gene details from Orthogroups.tsv"""
        try:
            # 查找对应的同源群行 | Find corresponding orthogroup row
            og_row = self.orthogroups_data[self.orthogroups_data.iloc[:, 0] == og_id]
            
            if og_row.empty:
                return [], []
            
            row = og_row.iloc[0]
            
            # 收集所有基因ID和对应基因组 | Collect all gene IDs and corresponding genomes
            all_genes = []
            all_genomes = []
            
            # 遍历每个基因组列 | Iterate through each genome column
            for genome_name in self.genome_names:
                if genome_name in self.orthogroups_data.columns:
                    col_idx = self.orthogroups_data.columns.get_loc(genome_name)
                    cell_value = row.iloc[col_idx]
                    
                    if pd.notna(cell_value) and str(cell_value).strip() and str(cell_value) != 'nan':
                        gene_string = str(cell_value).strip()
                        # 用逗号分割基因，并去除空格 | Split genes by comma and remove spaces
                        gene_list = [g.strip() for g in gene_string.split(',') if g.strip()]
                        for gene_id in gene_list:
                            all_genes.append(gene_id)
                            all_genomes.append(genome_name)
            
            # 确保返回的是列表 | Ensure returns are lists
            return list(all_genes), list(all_genomes)
            
        except Exception as e:
            self.logger.error(f"处理同源群 {og_id} 时出错: {e}")
            return [], []
    
    def save_classification_results(self, classification_results: Dict, output_dir: Path):
        """保存分类结果 | Save classification results"""
        self.logger.info("保存泛基因组分类结果 | Saving pangenome classification results")
        
        # 创建输出文件 | Create output file
        output_file = output_dir / "pangenome_gene_families.txt"
        
        with open(output_file, 'w', encoding='utf-8') as f:
            # 写入表头 | Write header
            f.write("Category\tOrthogroup_ID\tGene_ID\tGenome_Name\n")
            
            # 写入各类别的基因信息 | Write gene information for each category
            # for category in ['core', 'softcore', 'dispensable', 'private']:
            for category in ['core', 'softcore', 'dispensable', 'private', 'single_copy']:
                for og_id, gene_ids, genome_names in classification_results[category]:
                    for gene_id, genome_name in zip(gene_ids, genome_names):
                        f.write(f"{category}\t{og_id}\t{gene_id}\t{genome_name}\n")
        
        self.logger.info(f"分类结果已保存 | Classification results saved: {output_file}")
        
        # 保存分类统计摘要 | Save classification summary
        self._save_classification_summary(classification_results, output_dir)
        
        return output_file
    
    def _save_classification_summary(self, classification_results: Dict, output_dir: Path):
        """保存分类统计摘要 | Save classification summary"""
        summary_file = output_dir / "pangenome_classification_summary.txt"
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("泛基因组分类统计摘要 | Pangenome Classification Summary\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"总基因组数量 | Total genomes: {self.total_genomes}\n")
            f.write(f"Softcore缺失阈值 | Softcore missing threshold: <={self.config.softcore_missing_threshold}\n")
            f.write(f"Dispensable缺失阈值 | Dispensable missing threshold: >{self.config.dispensable_missing_threshold}\n\n")
            
            f.write("分类定义 | Classification Definitions:\n")
            f.write("-" * 40 + "\n")
            f.write("Core: 存在于所有基因组 | Present in all genomes\n")
            f.write(f"Softcore: 最多缺失{self.config.softcore_missing_threshold}个基因组 | Missing in ≤{self.config.softcore_missing_threshold} genomes\n")
            f.write(f"Dispensable: 缺失超过{self.config.dispensable_missing_threshold}个基因组 | Missing in >{self.config.dispensable_missing_threshold} genomes\n")
            f.write("Private: 仅存在于单个基因组 | Present in only one genome\n\n")
            
            f.write("分类结果 | Classification Results:\n")
            f.write("-" * 40 + "\n")
            
            total_orthogroups = sum(len(results) for results in classification_results.values())
            total_genes = sum(sum(len(gene_info[1]) for gene_info in results) 
                            for results in classification_results.values())
            
            for category in ['core', 'softcore', 'dispensable', 'private']:
                og_count = len(classification_results[category])
                gene_count = sum(len(gene_info[1]) for gene_info in classification_results[category])
                
                og_percentage = (og_count / total_orthogroups * 100) if total_orthogroups > 0 else 0
                gene_percentage = (gene_count / total_genes * 100) if total_genes > 0 else 0
                
                category_name = {
                    'core': '核心基因 | Core genes',
                    'softcore': '软核心基因 | Softcore genes', 
                    'dispensable': '非必需基因 | Dispensable genes',
                    'private': '私有基因 | Private genes'
                }[category]
                
                f.write(f"{category_name}:\n")
                f.write(f"  同源基因群数量 | Orthogroups: {og_count:,} ({og_percentage:.2f}%)\n")
                f.write(f"  基因数量 | Genes: {gene_count:,} ({gene_percentage:.2f}%)\n\n")
            
            f.write(f"总计 | Total:\n")
            f.write(f"  同源基因群总数 | Total orthogroups: {total_orthogroups:,}\n")
            f.write(f"  基因总数 | Total genes: {total_genes:,}\n")
        
        self.logger.info(f"分类摘要已保存 | Classification summary saved: {summary_file}")
    
    # def get_genome_gene_counts(self, classification_results: Dict) -> Dict[str, Dict[str, int]]:
    #     """获取每个基因组的基因数量统计 | Get gene count statistics for each genome"""
    #     genome_stats = {}
        
    #     for genome_name in self.genome_names:
    #         genome_stats[genome_name] = {
    #             'core': 0,
    #             'softcore': 0, 
    #             'dispensable': 0,
    #             'private': 0,
    #             'single_copy': [],    # 单拷贝基因 | Single copy genes
    #             'total': 0
    #         }
        
    #     # 统计每个类别中每个基因组的基因数量 | Count genes per genome in each category
    #     for category, gene_infos in classification_results.items():
    #         for og_id, gene_ids, genome_names in gene_infos:
    #             for gene_id, genome_name in zip(gene_ids, genome_names):
    #                 if genome_name in genome_stats:
    #                     genome_stats[genome_name][category] += 1
    #                     genome_stats[genome_name]['total'] += 1
        
    #     return genome_stats

    def get_genome_gene_counts(self, classification_results: Dict) -> Dict[str, Dict[str, int]]:
        """获取每个基因组的基因数量统计 | Get gene count statistics for each genome"""
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
        
        # 统计每个类别中每个基因组的基因数量 | Count genes per genome in each category
        for category, gene_infos in classification_results.items():
            for gene_info in gene_infos:
                try:
                    # 确保gene_info是3元组
                    if len(gene_info) != 3:
                        continue
                        
                    og_id, gene_ids, genome_names = gene_info
                    
                    # 确保gene_ids和genome_names是列表
                    if not isinstance(gene_ids, (list, tuple)):
                        continue
                    if not isinstance(genome_names, (list, tuple)):
                        continue
                        
                    # 确保两个列表长度相同
                    if len(gene_ids) != len(genome_names):
                        continue
                        
                    for gene_id, genome_name in zip(gene_ids, genome_names):
                        if genome_name in genome_stats:
                            genome_stats[genome_name][category] += 1
                            genome_stats[genome_name]['total'] += 1
                            
                except Exception as e:
                    # 跳过有问题的数据
                    continue
        
        return genome_stats
    
    def calculate_frequency_distribution(self, classification_results: Dict) -> Dict[str, Dict[int, int]]:
        """计算每个分类中同源群的出现频率分布 | Calculate frequency distribution of orthogroups in each category"""
        frequency_distributions = {}
        
        for category, gene_infos in classification_results.items():
            frequency_dist = defaultdict(int)
            
            for og_id, gene_ids, genome_names in gene_infos:
                # 计算该同源群在多少个基因组中出现 | Count in how many genomes this orthogroup appears
                unique_genomes = set(genome_names)
                frequency = len(unique_genomes)
                frequency_dist[frequency] += 1
            
            frequency_distributions[category] = dict(frequency_dist)
        
        return frequency_distributions
    
    def _is_single_copy_orthogroup(self, og_id: str) -> bool:
        """检查是否为单拷贝同源群 | Check if single copy orthogroup"""
        if og_id not in self.gene_count_data.index:
            return False
        
        row = self.gene_count_data.loc[og_id]
        # 检查所有基因组是否都有且只有1个基因 | Check if all genomes have exactly 1 gene
        return all(count == 1 for count in row)