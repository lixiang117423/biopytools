"""
📊 VCF基因型统计结果处理模块 | VCF Genotype Statistics Results Processing Module
"""

import os
import pandas as pd
from collections import defaultdict

class ResultsExporter:
    """💾 结果导出器 | Results Exporter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def export_summary_statistics(self, stats_results: dict, detailed_stats: dict):
        """📊 导出汇总统计结果 | Export summary statistics"""
        if not self.config.output_summary:
            return
        
        self.logger.info("💾 导出汇总统计结果 | Exporting summary statistics")
        
        # 📋 创建汇总DataFrame | Create summary DataFrame
        summary_data = []
        for sample_name, stats in stats_results.items():
            summary_data.append({
                'Sample': sample_name,
                'Total_Sites': stats['total_sites'],
                'Valid_Calls': stats['valid_calls'],
                'Missing_Calls': stats['missing_calls'],
                'Call_Rate': f"{stats['call_rate']:.4f}",
                'Missing_Rate': f"{stats['missing_rate']:.4f}",
                'Hom_Ref_Count': stats['hom_ref_count'],
                'Het_Count': stats['het_count'],
                'Hom_Alt_Count': stats['hom_alt_count'],
                'Hom_Ref_Rate': f"{stats['hom_ref_rate']:.4f}",
                'Het_Rate': f"{stats['het_rate']:.4f}",
                'Hom_Alt_Rate': f"{stats['hom_alt_rate']:.4f}",
                'Heterozygosity_Rate': f"{stats['heterozygosity_rate']:.4f}",
                'Homozygosity_Rate': f"{stats['homozygosity_rate']:.4f}",
            })
        
        summary_df = pd.DataFrame(summary_data)
        
        # 💾 保存汇总统计 | Save summary statistics
        summary_file = os.path.join(self.config.output_dir, "genotype_summary_statistics.txt")
        summary_df.to_csv(summary_file, sep='\t', index=False)
        self.logger.info(f"✅ 汇总统计已保存 | Summary statistics saved: {summary_file}")
        
        # 💾 保存简化版本（仅包含主要比率） | Save simplified version (main rates only)
        simple_data = []
        for sample_name, stats in stats_results.items():
            simple_data.append({
                'Sample': sample_name,
                'Total_SNPs': stats['total_sites'],
                'Hom_Ref_0/0': f"{stats['hom_ref_count']} ({stats['hom_ref_rate']:.4f})",
                'Het_0/1': f"{stats['het_count']} ({stats['het_rate']:.4f})",
                'Hom_Alt_1/1': f"{stats['hom_alt_count']} ({stats['hom_alt_rate']:.4f})",
                'Missing_./.': f"{stats['missing_calls']} ({stats['missing_rate']:.4f})",
                'Heterozygosity_Rate': f"{stats['heterozygosity_rate']:.4f}",
                'Homozygosity_Rate': f"{stats['homozygosity_rate']:.4f}",
            })
        
        simple_df = pd.DataFrame(simple_data)
        simple_file = os.path.join(self.config.output_dir, "genotype_rates_simple.txt")
        simple_df.to_csv(simple_file, sep='\t', index=False)
        self.logger.info(f"✅ 简化统计已保存 | Simplified statistics saved: {simple_file}")
    
    def export_detailed_statistics(self, detailed_stats: dict):
        """📋 导出详细统计结果 | Export detailed statistics"""
        if not self.config.output_detailed:
            return
        
        self.logger.info("📋 导出详细统计结果 | Exporting detailed statistics")
        
        # 🔍 收集所有基因型类别 | Collect all genotype categories
        all_genotypes = set()
        for sample_stats in detailed_stats.values():
            all_genotypes.update(sample_stats.keys())
        
        all_genotypes = sorted(list(all_genotypes))
        
        # 📊 创建详细统计DataFrame | Create detailed statistics DataFrame
        detailed_data = []
        for sample_name, stats in detailed_stats.items():
            total_calls = sum(stats.values())
            row = {'Sample': sample_name, 'Total_Calls': total_calls}
            
            for gt in all_genotypes:
                count = stats.get(gt, 0)
                rate = count / total_calls if total_calls > 0 else 0
                row[f'{gt}_Count'] = count
                row[f'{gt}_Rate'] = f"{rate:.4f}"
            
            detailed_data.append(row)
        
        detailed_df = pd.DataFrame(detailed_data)
        detailed_file = os.path.join(self.config.output_dir, "genotype_detailed_statistics.txt")
        detailed_df.to_csv(detailed_file, sep='\t', index=False)
        self.logger.info(f"✅ 详细统计已保存 | Detailed statistics saved: {detailed_file}")
    
    def export_per_sample_files(self, stats_results: dict, detailed_stats: dict):
        """👥 为每个样本导出单独的统计文件 | Export separate statistics files for each sample"""
        self.logger.info("👥 为每个样本导出单独统计文件 | Exporting individual sample statistics")
        
        sample_dir = os.path.join(self.config.output_dir, "per_sample_stats")
        os.makedirs(sample_dir, exist_ok=True)
        
        for sample_name in stats_results.keys():
            sample_file = os.path.join(sample_dir, f"{sample_name}_genotype_stats.txt")
            
            with open(sample_file, 'w', encoding='utf-8') as f:
                f.write(f"📊 基因型统计报告 - 样本: {sample_name}\n")
                f.write(f"📊 Genotype Statistics Report - Sample: {sample_name}\n")
                f.write("=" * 60 + "\n\n")
                
                # 📋 基本统计 | Basic statistics
                stats = stats_results[sample_name]
                f.write("📋 基本统计 | Basic Statistics:\n")
                f.write(f"  📊 总位点数 | Total Sites: {stats['total_sites']}\n")
                f.write(f"  ✅ 有效调用 | Valid Calls: {stats['valid_calls']}\n")
                f.write(f"  ❓ 缺失调用 | Missing Calls: {stats['missing_calls']}\n")
                f.write(f"  📈 调用率 | Call Rate: {stats['call_rate']:.4f}\n")
                f.write(f"  📉 缺失率 | Missing Rate: {stats['missing_rate']:.4f}\n\n")
                
                # 🔢 基因型计数 | Genotype counts
                f.write("🔢 基因型计数 | Genotype Counts:\n")
                f.write(f"  🟢 参考纯合子 (0/0) | Hom Ref: {stats['hom_ref_count']}\n")
                f.write(f"  🟡 杂合子 (0/1) | Heterozygous: {stats['het_count']}\n")
                f.write(f"  🔴 变异纯合子 (1/1) | Hom Alt: {stats['hom_alt_count']}\n\n")
                
                # 📊 基因型比率 | Genotype rates
                f.write("📊 基因型比率 | Genotype Rates:\n")
                f.write(f"  🟢 参考纯合子率 | Hom Ref Rate: {stats['hom_ref_rate']:.4f}\n")
                f.write(f"  🟡 杂合子率 | Het Rate: {stats['het_rate']:.4f}\n")
                f.write(f"  🔴 变异纯合子率 | Hom Alt Rate: {stats['hom_alt_rate']:.4f}\n")
                f.write(f"  🌟 杂合率 | Heterozygosity Rate: {stats['heterozygosity_rate']:.4f}\n")
                f.write(f"  🔹 纯合率 | Homozygosity Rate: {stats['homozygosity_rate']:.4f}\n\n")
                
                # 📋 详细基因型分布 | Detailed genotype distribution
                if sample_name in detailed_stats:
                    f.write("📋 详细基因型分布 | Detailed Genotype Distribution:\n")
                    sample_detailed = detailed_stats[sample_name]
                    total = sum(sample_detailed.values())
                    
                    for gt in sorted(sample_detailed.keys()):
                        count = sample_detailed[gt]
                        rate = count / total if total > 0 else 0
                        f.write(f"  🧬 {gt}: {count} ({rate:.4f})\n")

class SummaryGenerator:
    """📝 总结生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_analysis_summary(self, stats_results: dict):
        """📄 生成分析总结报告 | Generate analysis summary report"""
        summary_file = os.path.join(self.config.output_dir, "analysis_summary.txt")
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("📊 VCF基因型统计分析总结报告\n")
            f.write("📊 VCF Genotype Statistics Analysis Summary Report\n")
            f.write("=" * 60 + "\n\n")
            
            # 📁 输入文件信息 | Input file information
            f.write("📁 输入文件 | Input Files:\n")
            f.write(f"  📄 VCF文件 | VCF file: {self.config.vcf_file}\n\n")
            
            # ⚙️ 分析参数 | Analysis parameters
            f.write("⚙️ 分析参数 | Analysis Parameters:\n")
            f.write(f"  📏 最小深度过滤 | Min depth filter: {self.config.min_depth}\n")
            f.write(f"  🎯 最小质量过滤 | Min quality filter: {self.config.min_qual}\n")
            f.write(f"  ❌ 排除缺失基因型 | Exclude missing genotypes: {'是 | Yes' if self.config.exclude_missing else '否 | No'}\n\n")
            
            # 📊 数据总览 | Data overview
            f.write("📊 数据总览 | Data Overview:\n")
            f.write(f"  👥 样本数量 | Number of samples: {len(self.config.sample_names) if self.config.sample_names else 0}\n")
            f.write(f"  🧬 总变异位点数 | Total variant sites: {self.config.total_snps}\n\n")
            
            # 📈 样本统计概览 | Sample statistics overview
            if stats_results:
                f.write("📈 样本统计概览 | Sample Statistics Overview:\n")
                
                # 📊 计算全体样本的平均统计 | Calculate average statistics across all samples
                total_samples = len(stats_results)
                avg_call_rate = sum(s['call_rate'] for s in stats_results.values()) / total_samples
                avg_het_rate = sum(s['heterozygosity_rate'] for s in stats_results.values()) / total_samples
                avg_hom_rate = sum(s['homozygosity_rate'] for s in stats_results.values()) / total_samples
                
                f.write(f"  📈 平均调用率 | Average call rate: {avg_call_rate:.4f}\n")
                f.write(f"  🌟 平均杂合率 | Average heterozygosity rate: {avg_het_rate:.4f}\n")
                f.write(f"  🔹 平均纯合率 | Average homozygosity rate: {avg_hom_rate:.4f}\n\n")
                
                # ⚠️ 找出异常样本 | Identify outlier samples
                f.write("⚠️ 潜在异常样本 | Potential Outlier Samples:\n")
                low_call_samples = [name for name, stats in stats_results.items() if stats['call_rate'] < 0.9]
                high_het_samples = [name for name, stats in stats_results.items() if stats['heterozygosity_rate'] > 0.7]
                
                if low_call_samples:
                    f.write(f"  📉 低调用率样本 (<0.9) | Low call rate samples: {', '.join(low_call_samples[:5])}{'...' if len(low_call_samples) > 5 else ''}\n")
                if high_het_samples:
                    f.write(f"  ⬆️ 高杂合率样本 (>0.7) | High heterozygosity samples: {', '.join(high_het_samples[:5])}{'...' if len(high_het_samples) > 5 else ''}\n")
                
                if not low_call_samples and not high_het_samples:
                    f.write("  ✅ 未发现明显异常样本 | No obvious outlier samples detected\n")
            
            f.write(f"\n📁 输出文件 | Output Files:\n")
            f.write(f"  📊 - genotype_summary_statistics.txt: 汇总统计表 | Summary statistics table\n")
            f.write(f"  📋 - genotype_rates_simple.txt: 简化比率表 | Simplified rates table\n")
            if self.config.output_detailed:
                f.write(f"  📄 - genotype_detailed_statistics.txt: 详细统计表 | Detailed statistics table\n")
                f.write(f"  📂 - per_sample_stats/: 每个样本的详细报告 | Individual sample reports\n")
        
        self.logger.info(f"✅ 分析总结报告已生成 | Analysis summary report generated: {summary_file}")