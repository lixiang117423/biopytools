# # ===== FILE: orthofinder_pangenome/results.py =====
# """
# 结果处理和报告生成模块 | Results Processing and Report Generation Module
# """

# import pandas as pd
# from pathlib import Path
# from typing import Dict, List, Optional
# from datetime import datetime

# class ResultsProcessor:
#     """结果处理器 | Results Processor"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def generate_comprehensive_report(self, classification_results: Dict, genome_stats: Dict, 
#                                     frequency_distributions: Dict, rarefaction_results: Optional[Dict],
#                                     output_dir: Path) -> Path:
#         """生成综合分析报告 | Generate comprehensive analysis report"""
#         self.logger.info("生成综合分析报告 | Generating comprehensive analysis report")
        
#         report_file = output_dir / "pangenome_analysis_report.txt"
        
#         with open(report_file, 'w', encoding='utf-8') as f:
#             self._write_report_header(f)
#             self._write_input_summary(f)
#             self._write_analysis_parameters(f)
#             self._write_pangenome_overview(f, classification_results)
#             self._write_frequency_distribution_summary(f, frequency_distributions)
#             if rarefaction_results:
#                 self._write_rarefaction_summary(f, rarefaction_results)
#             self._write_genome_statistics(f, genome_stats)
#             self._write_output_files_summary(f, output_dir)
#             self._write_interpretation_guide(f)
        
#         self.logger.info(f"综合报告已生成 | Comprehensive report generated: {report_file}")
#         return report_file
    
#     def _write_report_header(self, f):
#         """写入报告头部 | Write report header"""
#         f.write("OrthoFinder泛基因组分析报告 | OrthoFinder Pangenome Analysis Report\n")
#         f.write("=" * 80 + "\n")
#         f.write(f"分析时间 | Analysis Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
#         f.write(f"分析工具 | Analysis Tool: OrthoFinder Pangenome Analyzer v1.0\n")
#         f.write(f"分析项目 | Project Name: {self.config.project_name}\n")
#         f.write("\n")
    
#     def _write_input_summary(self, f):
#         """写入输入数据摘要 | Write input data summary"""
#         f.write("输入数据摘要 | Input Data Summary\n")
#         f.write("-" * 40 + "\n")
#         f.write(f"输入目录 | Input Directory: {self.config.input_dir}\n")
#         f.write(f"输出目录 | Output Directory: {self.config.output_dir}\n")
#         f.write(f"序列类型 | Sequence Type: {self.config.sequence_type}\n")
#         f.write("\n")
    
#     def _write_analysis_parameters(self, f):
#         """写入分析参数 | Write analysis parameters"""
#         f.write("分析参数 | Analysis Parameters\n")
#         f.write("-" * 40 + "\n")
#         f.write(f"线程数 | Threads: {self.config.threads}\n")
#         f.write(f"搜索程序 | Search Program: {self.config.search_program}\n")
#         f.write(f"MCL Inflation: {self.config.mcl_inflation}\n")
#         f.write(f"Softcore缺失阈值 | Softcore Missing Threshold: <={self.config.softcore_missing_threshold}\n")
#         f.write(f"Dispensable缺失阈值 | Dispensable Missing Threshold: >{self.config.dispensable_missing_threshold}\n")
#         f.write(f"基础分析模式 | Basic Analysis Only: {'是 | Yes' if self.config.basic_analysis_only else '否 | No'}\n")
#         f.write(f"稀释曲线分析 | Rarefaction Analysis: {'是 | Yes' if self.config.enable_rarefaction else '否 | No'}\n")
#         f.write("\n")
    
#     def _write_pangenome_overview(self, f, classification_results: Dict):
#         """写入泛基因组概览 | Write pangenome overview"""
#         f.write("泛基因组概览 | Pangenome Overview\n")
#         f.write("-" * 40 + "\n")
        
#         total_orthogroups = sum(len(results) for results in classification_results.values())
#         total_genes = sum(sum(len(gene_info[1]) for gene_info in results) 
#                          for results in classification_results.values())
        
#         f.write(f"同源基因群总数 | Total Orthogroups: {total_orthogroups:,}\n")
#         f.write(f"基因总数 | Total Genes: {total_genes:,}\n")
#         f.write("\n")
        
#         f.write("分类统计 | Classification Statistics:\n")
#         # for category in ['core', 'softcore', 'dispensable', 'private']:
#         for category in ['core', 'softcore', 'dispensable', 'private', 'single_copy']:
#             og_count = len(classification_results[category])
#             gene_count = sum(len(gene_info[1]) for gene_info in classification_results[category])
#             og_percentage = (og_count / total_orthogroups * 100) if total_orthogroups > 0 else 0
#             gene_percentage = (gene_count / total_genes * 100) if total_genes > 0 else 0
            
#             # category_name = {
#             #     'core': '核心基因 | Core Genes',
#             #     'softcore': '软核心基因 | Softcore Genes',
#             #     'dispensable': '非必需基因 | Dispensable Genes', 
#             #     'private': '私有基因 | Private Genes'
#             # }[category]
#             category_name = {
#                 'core': '核心基因 | Core genes',
#                 'softcore': '软核心基因 | Softcore genes',
#                 'dispensable': '非必需基因 | Dispensable genes',
#                 'private': '私有基因 | Private genes',
#                 'single_copy': '单拷贝基因 | Single copy genes'
#             }[category]
            
#             f.write(f"  {category_name}:\n")
#             f.write(f"    同源基因群: {og_count:,} ({og_percentage:.2f}%)\n")
#             f.write(f"    基因数量: {gene_count:,} ({gene_percentage:.2f}%)\n")
#         f.write("\n")
    
#     def _write_frequency_distribution_summary(self, f, frequency_distributions: Dict):
#         """写入频率分布摘要 | Write frequency distribution summary"""
#         f.write("频率分布摘要 | Frequency Distribution Summary\n")
#         f.write("-" * 40 + "\n")
        
#         for category in ['core', 'softcore', 'dispensable', 'private']:
#             freq_dist = frequency_distributions.get(category, {})
#             if freq_dist:
#                 category_name = {
#                     'core': '核心基因 | Core Genes',
#                     'softcore': '软核心基因 | Softcore Genes',
#                     'dispensable': '非必需基因 | Dispensable Genes',
#                     'private': '私有基因 | Private Genes'
#                 }[category]
                
#                 f.write(f"{category_name} 频率分布:\n")
#                 for frequency in sorted(freq_dist.keys()):
#                     count = freq_dist[frequency]
#                     f.write(f"  出现在{frequency}个基因组: {count}个同源群\n")
#                 f.write("\n")
    
#     # ===== 替换文件: orthofinder_pangenome/results.py =====

#     def _write_rarefaction_summary(self, f, rarefaction_results: Dict):
#         """写入稀释分析摘要 | Write rarefaction analysis summary"""
#         f.write("稀释曲线分析摘要 | Rarefaction Curve Analysis Summary\n")
#         f.write("-" * 40 + "\n")
        
#         # 从详细结果中提取数据
#         df_detailed = pd.DataFrame(rarefaction_results)
        
#         # 获取样本大小范围
#         min_sample_size = df_detailed['sample_sizes'].min()
#         max_sample_size = df_detailed['sample_sizes'].max()
        
#         # 计算最大Pan基因组和最小Core基因组的均值
#         max_pan_mean = df_detailed.groupby('sample_sizes')['pan_counts'].mean().max()
#         min_core_mean = df_detailed.groupby('sample_sizes')['core_counts'].mean().min()
        
#         # 获取全部样本时的Pan和Core大小 (只有一个迭代，所以直接取值)
#         full_dataset_pan = df_detailed[df_detailed['sample_sizes'] == max_sample_size]['pan_counts'].iloc[0]
#         full_dataset_core = df_detailed[df_detailed['sample_sizes'] == max_sample_size]['core_counts'].iloc[0]
        
#         f.write(f"样本大小范围 | Sample size range: {min_sample_size} - {max_sample_size}\n")
#         f.write(f"最大Pan基因组大小 (均值) | Maximum pan-genome size (mean): {max_pan_mean:.0f}\n")
#         f.write(f"最小Core基因组大小 (均值) | Minimum core-genome size (mean): {min_core_mean:.0f}\n")
#         f.write(f"全部样本Pan大小 | Full dataset pan size: {full_dataset_pan:.0f}\n")
#         f.write(f"全部样本Core大小 | Full dataset core size: {full_dataset_core:.0f}\n")
#         f.write("\n")
    
#     def _write_genome_statistics(self, f, genome_stats: Dict):
#         """写入基因组统计信息 | Write genome statistics"""
#         f.write("基因组统计信息 | Genome Statistics\n")
#         f.write("-" * 40 + "\n")
        
#         # 创建统计表格 | Create statistics table
#         f.write(f"{'基因组名称':<20} {'核心':<8} {'软核心':<8} {'非必需':<8} {'私有':<8} {'总计':<8}\n")
#         f.write(f"{'Genome Name':<20} {'Core':<8} {'Soft':<8} {'Disp.':<8} {'Priv.':<8} {'Total':<8}\n")
#         f.write("-" * 70 + "\n")
        
#         for genome_name, stats in genome_stats.items():
#             f.write(f"{genome_name:<20} {stats['core']:<8} {stats['softcore']:<8} "
#                    f"{stats['dispensable']:<8} {stats['private']:<8} {stats['total']:<8}\n")
#         f.write("\n")
    
#     def _write_output_files_summary(self, f, output_dir: Path):
#         """写入输出文件摘要 | Write output files summary"""
#         f.write("输出文件摘要 | Output Files Summary\n")
#         f.write("-" * 40 + "\n")
        
#         output_files = [
#             ("pangenome_gene_families.txt", "泛基因组基因家族详细信息 | Detailed pangenome gene families"),
#             ("gene_families_detailed.tsv", "基因家族详细表格 | Detailed gene families table"),
#             ("pangenome_classification_summary.txt", "泛基因组分类统计摘要 | Pangenome classification summary"),
#             ("rarefaction_curve_data.tsv", "稀释曲线数据 | Rarefaction curve data"),
#             ("pangenome_rarefaction_curve.png", "稀释曲线图 | Rarefaction curve plot"),
#             ("pangenome_classification_pie.png", "分类饼图 | Classification pie chart"),
#             ("frequency_distribution.png", "频率分布图 | Frequency distribution plot"),
#             ("pangenome_analysis_report.txt", "综合分析报告 | Comprehensive analysis report"),
#             ("pangenome_analysis.log", "分析日志文件 | Analysis log file")
#         ]
        
#         for filename, description in output_files:
#             file_path = output_dir / filename
#             status = "存在" if file_path.exists() else "缺失"
#             f.write(f"{status} {filename}: {description}\n")
#         f.write("\n")
    
#     def _write_interpretation_guide(self, f):
#         """写入结果解读指南 | Write interpretation guide"""
#         f.write("结果解读指南 | Interpretation Guide\n")
#         f.write("-" * 40 + "\n")
        
#         f.write("泛基因组分析解读要点:\n")
#         f.write("Key Points for Pangenome Analysis Interpretation:\n\n")
        
#         f.write("1. 核心基因组 (Core Genome):\n")
#         f.write("   - 所有样本共享的基因家族\n")
#         f.write("   - 代表基本生命功能\n")
#         f.write("   - Gene families shared by all samples\n")
#         f.write("   - Represent essential life functions\n\n")
        
#         f.write("2. 软核心基因组 (Softcore Genome):\n")
#         f.write(f"   - 最多缺失{self.config.softcore_missing_threshold}个基因组的基因家族\n")
#         f.write("   - 重要但非绝对必需的功能\n")
#         f.write(f"   - Gene families missing in ≤{self.config.softcore_missing_threshold} genomes\n")
#         f.write("   - Important but not absolutely essential functions\n\n")
        
#         f.write("3. 非必需基因组 (Dispensable Genome):\n")
#         f.write(f"   - 缺失超过{self.config.dispensable_missing_threshold}个基因组的基因家族\n")
#         f.write("   - 提供功能多样性和环境适应性\n")
#         f.write(f"   - Gene families missing in >{self.config.dispensable_missing_threshold} genomes\n")
#         f.write("   - Provide functional diversity and environmental adaptation\n\n")
        
#         f.write("4. 私有基因组 (Private Genome):\n")
#         f.write("   - 仅在单个基因组中存在的基因家族\n")
#         f.write("   - 可能代表独特的适应性特征\n")
#         f.write("   - Gene families present in only one genome\n")
#         f.write("   - May represent unique adaptive features\n\n")
        
#         f.write("稀释曲线解读:\n")
#         f.write("Rarefaction Curve Interpretation:\n")
#         f.write("- Pan曲线上升表示新基因家族的发现\n")
#         f.write("- Core曲线下降表示共有基因家族的减少\n")
#         f.write("- Pan curve rising indicates discovery of new gene families\n")
#         f.write("- Core curve declining indicates reduction of shared gene families\n\n")
        
#         f.write("建议后续分析:\n")
#         f.write("Recommended Follow-up Analyses:\n")
#         f.write("- 功能富集分析 | Functional enrichment analysis\n")
#         f.write("- 系统发育关联分析 | Phylogenetic association analysis\n") 
#         f.write("- 环境因子相关性分析 | Environmental factor correlation analysis\n")
#         f.write("- 基因家族进化分析 | Gene family evolution analysis\n")
    
#     # def save_gene_families_table(self, classification_results: Dict, output_dir: Path) -> Path:
#     #     """保存基因家族表格 | Save gene families table"""
#     #     self.logger.info("保存基因家族详细表格 | Saving detailed gene families table")
        
#     #     # 准备数据 | Prepare data
#     #     data_rows = []
#     #     for category in ['core', 'softcore', 'dispensable', 'private']:
#     #         for og_id, gene_ids, genome_names in classification_results[category]:
#     #             for gene_id, genome_name in zip(gene_ids, genome_names):
#     #                 data_rows.append({
#     #                     'Category': category,
#     #                     'Orthogroup_ID': og_id,
#     #                     'Gene_ID': gene_id,
#     #                     'Genome_Name': genome_name
#     #                 })
        
#     #     # 创建DataFrame并保存 | Create DataFrame and save
#     #     df = pd.DataFrame(data_rows)
#     #     output_file = output_dir / "gene_families_detailed.tsv"
#     #     df.to_csv(output_file, sep='\t', index=False)
        
#     #     self.logger.info(f"基因家族表格已保存 | Gene families table saved: {output_file}")
#     #     return output_file
    
#     def save_gene_families_table(self, classification_results: Dict, output_dir: Path) -> Path:
#         """保存基因家族表格 | Save gene families table"""
#         self.logger.info("保存基因家族详细表格 | Saving detailed gene families table")
        
#         # 准备数据 | Prepare data
#         data_rows = []
#         for category in ['core', 'softcore', 'dispensable', 'private', 'single_copy']:
#             if category not in classification_results:
#                 continue
                
#             for gene_info in classification_results[category]:
#                 try:
#                     # 安全解包元组 | Safe tuple unpacking
#                     if len(gene_info) == 3:
#                         og_id, gene_ids, genome_names = gene_info
#                     else:
#                         self.logger.warning(f"跳过格式异常的基因信息 | Skipping malformed gene info: {gene_info}")
#                         continue
                    
#                     # 确保gene_ids和genome_names是列表 | Ensure gene_ids and genome_names are lists
#                     if not isinstance(gene_ids, (list, tuple)):
#                         self.logger.warning(f"基因ID不是列表格式，跳过 | Gene IDs not in list format, skipping: {og_id}")
#                         continue
                        
#                     if not isinstance(genome_names, (list, tuple)):
#                         self.logger.warning(f"基因组名称不是列表格式，跳过 | Genome names not in list format, skipping: {og_id}")
#                         continue
                    
#                     # 检查长度一致性 | Check length consistency
#                     if len(gene_ids) != len(genome_names):
#                         self.logger.warning(f"基因ID和基因组名称数量不匹配，跳过 | Gene IDs and genome names count mismatch, skipping: {og_id}")
#                         continue
                    
#                     # 添加数据行 | Add data rows
#                     for gene_id, genome_name in zip(gene_ids, genome_names):
#                         data_rows.append({
#                             'Category': category,
#                             'Orthogroup_ID': og_id,
#                             'Gene_ID': str(gene_id),
#                             'Genome_Name': str(genome_name)
#                         })
                        
#                 except Exception as e:
#                     self.logger.error(f"处理基因信息时出错 | Error processing gene info: {gene_info}, 错误: {e}")
#                     continue
        
#         # 创建DataFrame并保存 | Create DataFrame and save
#         if data_rows:
#             df = pd.DataFrame(data_rows)
#             output_file = output_dir / "gene_families_detailed.tsv"
#             df.to_csv(output_file, sep='\t', index=False)
            
#             self.logger.info(f"基因家族表格已保存 | Gene families table saved: {output_file}")
#             self.logger.info(f"保存了 {len(data_rows)} 行数据 | Saved {len(data_rows)} rows of data")
#         else:
#             self.logger.warning("没有有效数据可保存 | No valid data to save")
#             output_file = output_dir / "gene_families_detailed.tsv"
#             # 创建空文件 | Create empty file
#             with open(output_file, 'w') as f:
#                 f.write("Category\tOrthogroup_ID\tGene_ID\tGenome_Name\n")
        
#         return output_file
    
#     def save_frequency_distribution_table(self, frequency_distributions: Dict, output_dir: Path) -> Path:
#         """保存频率分布表格 | Save frequency distribution table"""
#         self.logger.info("保存频率分布表格 | Saving frequency distribution table")
        
#         # 准备数据 | Prepare data
#         all_frequencies = set()
#         for freq_dist in frequency_distributions.values():
#             all_frequencies.update(freq_dist.keys())
        
#         all_frequencies = sorted(all_frequencies)
        
#         # 创建表格 | Create table
#         data = {'Frequency': all_frequencies}
#         for category in ['core', 'softcore', 'dispensable', 'private']:
#             freq_dist = frequency_distributions.get(category, {})
#             data[category.capitalize()] = [freq_dist.get(freq, 0) for freq in all_frequencies]
        
#         df = pd.DataFrame(data)
#         output_file = output_dir / "frequency_distribution_table.tsv"
#         df.to_csv(output_file, sep='\t', index=False)
        
#         self.logger.info(f"频率分布表格已保存 | Frequency distribution table saved: {output_file}")
#         return output_file

# # ===== END FILE =====

"""
📝 结果处理和报告生成模块 | Results Processing and Report Generation Module
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime

class ResultsProcessor:
    """📝 结果处理器 | Results Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_comprehensive_report(self, classification_results: Dict, genome_stats: Dict, 
                                    frequency_distributions: Dict, rarefaction_results: Optional[Dict],
                                    output_dir: Path) -> Path:
        """📄 生成综合分析报告 | Generate comprehensive analysis report"""
        self.logger.info("📄 生成综合分析报告 | Generating comprehensive analysis report")
        
        report_file = output_dir / "pangenome_analysis_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            self._write_report_header(f)
            self._write_input_summary(f)
            self._write_analysis_parameters(f)
            self._write_pangenome_overview(f, classification_results)
            self._write_frequency_distribution_summary(f, frequency_distributions)
            if rarefaction_results:
                self._write_rarefaction_summary(f, rarefaction_results)
            self._write_genome_statistics(f, genome_stats)
            self._write_output_files_summary(f, output_dir)
            self._write_interpretation_guide(f)
        
        self.logger.info(f"✅ 综合报告已生成 | Comprehensive report generated: {report_file}")
        return report_file
    
    def _write_report_header(self, f):
        """📋 写入报告头部 | Write report header"""
        f.write("📝 OrthoFinder泛基因组分析报告 | OrthoFinder Pangenome Analysis Report\n")
        f.write("=" * 80 + "\n")
        f.write(f"⏰ 分析时间 | Analysis Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"🔧 分析工具 | Analysis Tool: OrthoFinder Pangenome Analyzer v1.0\n")
        f.write(f"🏷️ 分析项目 | Project Name: {self.config.project_name}\n")
        f.write("\n")
    
    def _write_input_summary(self, f):
        """📁 写入输入数据摘要 | Write input data summary"""
        f.write("📁 输入数据摘要 | Input Data Summary\n")
        f.write("-" * 40 + "\n")
        f.write(f"📂 输入目录 | Input Directory: {self.config.input_dir}\n")
        f.write(f"📂 输出目录 | Output Directory: {self.config.output_dir}\n")
        f.write(f"🧬 序列类型 | Sequence Type: {self.config.sequence_type}\n")
        f.write("\n")
    
    def _write_analysis_parameters(self, f):
        """⚙️ 写入分析参数 | Write analysis parameters"""
        f.write("⚙️ 分析参数 | Analysis Parameters\n")
        f.write("-" * 40 + "\n")
        f.write(f"🧵 线程数 | Threads: {self.config.threads}\n")
        f.write(f"🔍 搜索程序 | Search Program: {self.config.search_program}\n")
        f.write(f"🔧 MCL Inflation: {self.config.mcl_inflation}\n")
        f.write(f"🟠 Softcore缺失阈值 | Softcore Missing Threshold: <={self.config.softcore_missing_threshold}\n")
        f.write(f"🟡 Dispensable缺失阈值 | Dispensable Missing Threshold: >{self.config.dispensable_missing_threshold}\n")
        f.write(f"🔧 基础分析模式 | Basic Analysis Only: {'是 | Yes' if self.config.basic_analysis_only else '否 | No'}\n")
        f.write(f"📈 稀释曲线分析 | Rarefaction Analysis: {'是 | Yes' if self.config.enable_rarefaction else '否 | No'}\n")
        f.write("\n")
    
    def _write_pangenome_overview(self, f, classification_results: Dict):
        """🧬 写入泛基因组概览 | Write pangenome overview"""
        f.write("🧬 泛基因组概览 | Pangenome Overview\n")
        f.write("-" * 40 + "\n")
        
        total_orthogroups = sum(len(results) for results in classification_results.values())
        total_genes = sum(sum(len(gene_info[1]) for gene_info in results) 
                         for results in classification_results.values())
        
        f.write(f"📦 同源基因群总数 | Total Orthogroups: {total_orthogroups:,}\n")
        f.write(f"🧬 基因总数 | Total Genes: {total_genes:,}\n")
        f.write("\n")
        
        f.write("📊 分类统计 | Classification Statistics:\n")
        # for category in ['core', 'softcore', 'dispensable', 'private']:
        for category in ['core', 'softcore', 'dispensable', 'private', 'single_copy']:
            og_count = len(classification_results[category])
            gene_count = sum(len(gene_info[1]) for gene_info in classification_results[category])
            og_percentage = (og_count / total_orthogroups * 100) if total_orthogroups > 0 else 0
            gene_percentage = (gene_count / total_genes * 100) if total_genes > 0 else 0
            
            # category_name = {
            #     'core': '🔴 核心基因 | Core Genes',
            #     'softcore': '🟠 软核心基因 | Softcore Genes',
            #     'dispensable': '🟡 非必需基因 | Dispensable Genes', 
            #     'private': '🟢 私有基因 | Private Genes'
            # }[category]
            category_name = {
                'core': '🔴 核心基因 | Core genes',
                'softcore': '🟠 软核心基因 | Softcore genes',
                'dispensable': '🟡 非必需基因 | Dispensable genes',
                'private': '🟢 私有基因 | Private genes',
                'single_copy': '🔵 单拷贝基因 | Single copy genes'
            }[category]
            
            f.write(f"  {category_name}:\n")
            f.write(f"    📦 同源基因群: {og_count:,} ({og_percentage:.2f}%)\n")
            f.write(f"    🧬 基因数量: {gene_count:,} ({gene_percentage:.2f}%)\n")
        f.write("\n")
    
    def _write_frequency_distribution_summary(self, f, frequency_distributions: Dict):
        """📈 写入频率分布摘要 | Write frequency distribution summary"""
        f.write("📈 频率分布摘要 | Frequency Distribution Summary\n")
        f.write("-" * 40 + "\n")
        
        for category in ['core', 'softcore', 'dispensable', 'private']:
            freq_dist = frequency_distributions.get(category, {})
            if freq_dist:
                category_name = {
                    'core': '🔴 核心基因 | Core Genes',
                    'softcore': '🟠 软核心基因 | Softcore Genes',
                    'dispensable': '🟡 非必需基因 | Dispensable Genes',
                    'private': '🟢 私有基因 | Private Genes'
                }[category]
                
                f.write(f"{category_name} 📊 频率分布:\n")
                for frequency in sorted(freq_dist.keys()):
                    count = freq_dist[frequency]
                    f.write(f"  📊 出现在{frequency}个基因组: {count}个同源群\n")
                f.write("\n")
    
    # ===== 替换文件: orthofinder_pangenome/results.py =====

    def _write_rarefaction_summary(self, f, rarefaction_results: Dict):
        """📈 写入稀释分析摘要 | Write rarefaction analysis summary"""
        f.write("📈 稀释曲线分析摘要 | Rarefaction Curve Analysis Summary\n")
        f.write("-" * 40 + "\n")
        
        # 从详细结果中提取数据
        df_detailed = pd.DataFrame(rarefaction_results)
        
        # 获取样本大小范围
        min_sample_size = df_detailed['sample_sizes'].min()
        max_sample_size = df_detailed['sample_sizes'].max()
        
        # 计算最大Pan基因组和最小Core基因组的均值
        max_pan_mean = df_detailed.groupby('sample_sizes')['pan_counts'].mean().max()
        min_core_mean = df_detailed.groupby('sample_sizes')['core_counts'].mean().min()
        
        # 获取全部样本时的Pan和Core大小 (只有一个迭代，所以直接取值)
        full_dataset_pan = df_detailed[df_detailed['sample_sizes'] == max_sample_size]['pan_counts'].iloc[0]
        full_dataset_core = df_detailed[df_detailed['sample_sizes'] == max_sample_size]['core_counts'].iloc[0]
        
        f.write(f"📊 样本大小范围 | Sample size range: {min_sample_size} - {max_sample_size}\n")
        f.write(f"📈 最大Pan基因组大小 (均值) | Maximum pan-genome size (mean): {max_pan_mean:.0f}\n")
        f.write(f"📉 最小Core基因组大小 (均值) | Minimum core-genome size (mean): {min_core_mean:.0f}\n")
        f.write(f"📊 全部样本Pan大小 | Full dataset pan size: {full_dataset_pan:.0f}\n")
        f.write(f"📊 全部样本Core大小 | Full dataset core size: {full_dataset_core:.0f}\n")
        f.write("\n")
    
    def _write_genome_statistics(self, f, genome_stats: Dict):
        """📊 写入基因组统计信息 | Write genome statistics"""
        f.write("📊 基因组统计信息 | Genome Statistics\n")
        f.write("-" * 40 + "\n")
        
        # 创建统计表格 | Create statistics table
        f.write(f"{'基因组名称':<20} {'🔴核心':<8} {'🟠软核心':<8} {'🟡非必需':<8} {'🟢私有':<8} {'📊总计':<8}\n")
        f.write(f"{'Genome Name':<20} {'Core':<8} {'Soft':<8} {'Disp.':<8} {'Priv.':<8} {'Total':<8}\n")
        f.write("-" * 70 + "\n")
        
        for genome_name, stats in genome_stats.items():
            f.write(f"{genome_name:<20} {stats['core']:<8} {stats['softcore']:<8} "
                   f"{stats['dispensable']:<8} {stats['private']:<8} {stats['total']:<8}\n")
        f.write("\n")
    
    def _write_output_files_summary(self, f, output_dir: Path):
        """📁 写入输出文件摘要 | Write output files summary"""
        f.write("📁 输出文件摘要 | Output Files Summary\n")
        f.write("-" * 40 + "\n")
        
        output_files = [
            ("pangenome_gene_families.txt", "🧬 泛基因组基因家族详细信息 | Detailed pangenome gene families"),
            ("gene_families_detailed.tsv", "📋 基因家族详细表格 | Detailed gene families table"),
            ("pangenome_classification_summary.txt", "📊 泛基因组分类统计摘要 | Pangenome classification summary"),
            ("rarefaction_curve_data.tsv", "📈 稀释曲线数据 | Rarefaction curve data"),
            ("pangenome_rarefaction_curve.png", "📈 稀释曲线图 | Rarefaction curve plot"),
            ("pangenome_classification_pie.png", "🥧 分类饼图 | Classification pie chart"),
            ("frequency_distribution.png", "📊 频率分布图 | Frequency distribution plot"),
            ("pangenome_analysis_report.txt", "📝 综合分析报告 | Comprehensive analysis report"),
            ("pangenome_analysis.log", "📋 分析日志文件 | Analysis log file")
        ]
        
        for filename, description in output_files:
            file_path = output_dir / filename
            status = "✅ 存在" if file_path.exists() else "❌ 缺失"
            f.write(f"{status} {filename}: {description}\n")
        f.write("\n")
    
    def _write_interpretation_guide(self, f):
        """📚 写入结果解读指南 | Write interpretation guide"""
        f.write("📚 结果解读指南 | Interpretation Guide\n")
        f.write("-" * 40 + "\n")
        
        f.write("🧬 泛基因组分析解读要点:\n")
        f.write("Key Points for Pangenome Analysis Interpretation:\n\n")
        
        f.write("1. 🔴 核心基因组 (Core Genome):\n")
        f.write("   - 🧬 所有样本共享的基因家族\n")
        f.write("   - ⚡ 代表基本生命功能\n")
        f.write("   - Gene families shared by all samples\n")
        f.write("   - Represent essential life functions\n\n")
        
        f.write("2. 🟠 软核心基因组 (Softcore Genome):\n")
        f.write(f"   - 📊 最多缺失{self.config.softcore_missing_threshold}个基因组的基因家族\n")
        f.write("   - 🔧 重要但非绝对必需的功能\n")
        f.write(f"   - Gene families missing in ≤{self.config.softcore_missing_threshold} genomes\n")
        f.write("   - Important but not absolutely essential functions\n\n")
        
        f.write("3. 🟡 非必需基因组 (Dispensable Genome):\n")
        f.write(f"   - 📊 缺失超过{self.config.dispensable_missing_threshold}个基因组的基因家族\n")
        f.write("   - 🌱 提供功能多样性和环境适应性\n")
        f.write(f"   - Gene families missing in >{self.config.dispensable_missing_threshold} genomes\n")
        f.write("   - Provide functional diversity and environmental adaptation\n\n")
        
        f.write("4. 🟢 私有基因组 (Private Genome):\n")
        f.write("   - 🎯 仅在单个基因组中存在的基因家族\n")
        f.write("   - ✨ 可能代表独特的适应性特征\n")
        f.write("   - Gene families present in only one genome\n")
        f.write("   - May represent unique adaptive features\n\n")
        
        f.write("📈 稀释曲线解读:\n")
        f.write("Rarefaction Curve Interpretation:\n")
        f.write("- 📈 Pan曲线上升表示新基因家族的发现\n")
        f.write("- 📉 Core曲线下降表示共有基因家族的减少\n")
        f.write("- Pan curve rising indicates discovery of new gene families\n")
        f.write("- Core curve declining indicates reduction of shared gene families\n\n")
        
        f.write("🔬 建议后续分析:\n")
        f.write("Recommended Follow-up Analyses:\n")
        f.write("- 🧪 功能富集分析 | Functional enrichment analysis\n")
        f.write("- 🌳 系统发育关联分析 | Phylogenetic association analysis\n") 
        f.write("- 🌍 环境因子相关性分析 | Environmental factor correlation analysis\n")
        f.write("- 🧬 基因家族进化分析 | Gene family evolution analysis\n")
    
    # def save_gene_families_table(self, classification_results: Dict, output_dir: Path) -> Path:
    #     """💾 保存基因家族表格 | Save gene families table"""
    #     self.logger.info("💾 保存基因家族详细表格 | Saving detailed gene families table")
        
    #     # 准备数据 | Prepare data
    #     data_rows = []
    #     for category in ['core', 'softcore', 'dispensable', 'private']:
    #         for og_id, gene_ids, genome_names in classification_results[category]:
    #             for gene_id, genome_name in zip(gene_ids, genome_names):
    #                 data_rows.append({
    #                     'Category': category,
    #                     'Orthogroup_ID': og_id,
    #                     'Gene_ID': gene_id,
    #                     'Genome_Name': genome_name
    #                 })
        
    #     # 创建DataFrame并保存 | Create DataFrame and save
    #     df = pd.DataFrame(data_rows)
    #     output_file = output_dir / "gene_families_detailed.tsv"
    #     df.to_csv(output_file, sep='\t', index=False)
        
    #     self.logger.info(f"✅ 基因家族表格已保存 | Gene families table saved: {output_file}")
    #     return output_file
    
    def save_gene_families_table(self, classification_results: Dict, output_dir: Path) -> Path:
        """💾 保存基因家族表格 | Save gene families table"""
        self.logger.info("💾 保存基因家族详细表格 | Saving detailed gene families table")
        
        # 📋 准备数据 | Prepare data
        data_rows = []
        for category in ['core', 'softcore', 'dispensable', 'private', 'single_copy']:
            if category not in classification_results:
                continue
                
            for gene_info in classification_results[category]:
                try:
                    # 🔍 安全解包元组 | Safe tuple unpacking
                    if len(gene_info) == 3:
                        og_id, gene_ids, genome_names = gene_info
                    else:
                        self.logger.warning(f"⚠️ 跳过格式异常的基因信息 | Skipping malformed gene info: {gene_info}")
                        continue
                    
                    # ✅ 确保gene_ids和genome_names是列表 | Ensure gene_ids and genome_names are lists
                    if not isinstance(gene_ids, (list, tuple)):
                        self.logger.warning(f"⚠️ 基因ID不是列表格式，跳过 | Gene IDs not in list format, skipping: {og_id}")
                        continue
                        
                    if not isinstance(genome_names, (list, tuple)):
                        self.logger.warning(f"⚠️ 基因组名称不是列表格式，跳过 | Genome names not in list format, skipping: {og_id}")
                        continue
                    
                    # 🔍 检查长度一致性 | Check length consistency
                    if len(gene_ids) != len(genome_names):
                        self.logger.warning(f"⚠️ 基因ID和基因组名称数量不匹配，跳过 | Gene IDs and genome names count mismatch, skipping: {og_id}")
                        continue
                    
                    # ➕ 添加数据行 | Add data rows
                    for gene_id, genome_name in zip(gene_ids, genome_names):
                        data_rows.append({
                            'Category': category,
                            'Orthogroup_ID': og_id,
                            'Gene_ID': str(gene_id),
                            'Genome_Name': str(genome_name)
                        })
                        
                except Exception as e:
                    self.logger.error(f"❌ 处理基因信息时出错 | Error processing gene info: {gene_info}, 错误: {e}")
                    continue
        
        # 📊 创建DataFrame并保存 | Create DataFrame and save
        if data_rows:
            df = pd.DataFrame(data_rows)
            output_file = output_dir / "gene_families_detailed.tsv"
            df.to_csv(output_file, sep='\t', index=False)
            
            self.logger.info(f"✅ 基因家族表格已保存 | Gene families table saved: {output_file}")
            self.logger.info(f"📊 保存了 {len(data_rows)} 行数据 | Saved {len(data_rows)} rows of data")
        else:
            self.logger.warning("⚠️ 没有有效数据可保存 | No valid data to save")
            output_file = output_dir / "gene_families_detailed.tsv"
            # 📄 创建空文件 | Create empty file
            with open(output_file, 'w') as f:
                f.write("Category\tOrthogroup_ID\tGene_ID\tGenome_Name\n")
        
        return output_file
    
    def save_frequency_distribution_table(self, frequency_distributions: Dict, output_dir: Path) -> Path:
        """📊 保存频率分布表格 | Save frequency distribution table"""
        self.logger.info("📊 保存频率分布表格 | Saving frequency distribution table")
        
        # 📋 准备数据 | Prepare data
        all_frequencies = set()
        for freq_dist in frequency_distributions.values():
            all_frequencies.update(freq_dist.keys())
        
        all_frequencies = sorted(all_frequencies)
        
        # 📊 创建表格 | Create table
        data = {'Frequency': all_frequencies}
        for category in ['core', 'softcore', 'dispensable', 'private']:
            freq_dist = frequency_distributions.get(category, {})
            data[category.capitalize()] = [freq_dist.get(freq, 0) for freq in all_frequencies]
        
        df = pd.DataFrame(data)
        output_file = output_dir / "frequency_distribution_table.tsv"
        df.to_csv(output_file, sep='\t', index=False)
        
        self.logger.info(f"✅ 频率分布表格已保存 | Frequency distribution table saved: {output_file}")
        return output_file

# ===== END FILE =====