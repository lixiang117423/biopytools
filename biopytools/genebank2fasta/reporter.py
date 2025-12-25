"""
统计报告生成模块 📊 | Statistical Report Generation Module
"""

import os
from collections import defaultdict

class StatisticsReporter:
    """统计报告生成器 📈 | Statistics Reporter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_comprehensive_report(self, all_sample_data, output_dir):
        """生成详细的提取统计报告 📋 | Generate detailed extraction statistics report"""
        
        # 重组数据用于统计 | Reorganize data for statistics
        gene_sequences = defaultdict(list)
        for sample_data in all_sample_data:
            for seq_info in sample_data['cds_sequences']:
                gene_name = seq_info['gene']
                gene_sequences[gene_name].append(seq_info)
        
        report_file = f"{output_dir}/extraction_summary.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("🧬 GenBank序列提取报告 | GenBank Sequence Extraction Report\n")
            f.write("=" * 80 + "\n\n")
            
            # 总体统计 | Overall statistics
            f.write("📊 总体统计 | Overall Statistics:\n")
            f.write("-" * 40 + "\n")
            f.write(f"处理样本总数 | Total samples processed: {len(all_sample_data)}\n")
            f.write(f"提取基因种类数 | Gene types extracted: {len(gene_sequences)}\n")
            f.write(f"总序列条数 | Total sequences: {sum(len(seqs) for seqs in gene_sequences.values())}\n")
            f.write(f"使用线程数 | Threads used: {self.config.threads}\n\n")
            
            # 输出目录信息 | Output directory information
            f.write("📁 输出文件组织 | Output File Organization:\n")
            f.write("-" * 40 + "\n")
            if self.config.separate_by_sample:
                f.write(f"按样品分类 | By sample: {self.config.cds_dir}/by_sample/ 和 {self.config.pep_dir}/by_sample/\n")
            if self.config.separate_by_gene:
                f.write(f"按基因分类 | By gene: {self.config.cds_dir}/by_gene/ 和 {self.config.pep_dir}/by_gene/\n")
            f.write("\n")
            
            # 基因统计 | Gene statistics
            self.write_gene_statistics(f, gene_sequences)
            
            # 样本统计 | Sample statistics
            self.write_sample_statistics(f, all_sample_data)
            
            # 核心基因识别 | Core gene identification
            self.write_core_genes(f, gene_sequences, len(all_sample_data), output_dir)
        
        self.logger.info(f"📋 详细报告已生成: {report_file}")
    
    def write_gene_statistics(self, f, gene_sequences):
        """写入基因统计信息 🧬 | Write gene statistics"""
        f.write("🧬 基因统计 | Gene Statistics (按样本数排序 | sorted by sample count):\n")
        f.write("-" * 50 + "\n")
        f.write(f"{'基因名':<20} {'样本数':<10} {'平均CDS长度':<15} {'平均蛋白长度':<15}\n")
        f.write(f"{'Gene':<20} {'Samples':<10} {'Avg CDS len':<15} {'Avg Prot len':<15}\n")
        f.write("-" * 65 + "\n")
        
        gene_stats = []
        for gene_name, seq_list in gene_sequences.items():
            sample_count = len(seq_list)
            avg_cds_len = sum(seq['cds_length'] for seq in seq_list) // sample_count
            avg_prot_len = sum(seq['protein_length'] for seq in seq_list) // sample_count
            gene_stats.append((gene_name, sample_count, avg_cds_len, avg_prot_len))
        
        # 按样本数排序 | Sort by sample count
        gene_stats.sort(key=lambda x: x[1], reverse=True)
        
        for gene_name, count, avg_cds, avg_prot in gene_stats:
            f.write(f"{gene_name:<20} {count:<10} {avg_cds:<15} {avg_prot:<15}\n")
        
        f.write("\n")
    
    def write_sample_statistics(self, f, all_sample_data):
        """写入样本统计信息 🔬 | Write sample statistics"""
        f.write("🔬 样本统计 | Sample Statistics (前20个样本 | Top 20 samples):\n")
        f.write("-" * 50 + "\n")
        f.write(f"{'样本名':<35} {'提取基因数':<15}\n")
        f.write(f"{'Sample':<35} {'Genes':<15}\n")
        f.write("-" * 55 + "\n")
        
        # 按基因数排序 | Sort by gene count
        sorted_samples = sorted(all_sample_data, key=lambda x: x['genes_extracted'], reverse=True)
        
        for i, sample_data in enumerate(sorted_samples[:20]):
            f.write(f"{sample_data['sample']:<35} {sample_data['genes_extracted']:<15}\n")
        
        if len(sorted_samples) > 20:
            f.write(f"... 还有 {len(sorted_samples)-20} 个样本 | ... {len(sorted_samples)-20} more samples\n")
        
        f.write("\n")
    
    def write_core_genes(self, f, gene_sequences, total_samples, output_dir):
        """写入核心基因信息 💎 | Write core gene information"""
        core_threshold = int(total_samples * 0.9)
        
        f.write(f"💎 核心基因 | Core Genes (出现在>{core_threshold}个样本中 | present in >{core_threshold} samples):\n")
        f.write("-" * 50 + "\n")
        
        core_genes = []
        for gene_name, seq_list in gene_sequences.items():
            if len(seq_list) >= core_threshold:
                core_genes.append((gene_name, len(seq_list)))
        
        core_genes.sort(key=lambda x: x[1], reverse=True)
        for gene, count in core_genes:
            f.write(f"{gene}: {count}/{total_samples} 样本 | samples\n")
        
        f.write(f"\n核心基因总数 | Total core genes: {len(core_genes)}\n")
        
        # 生成核心基因列表文件 | Generate core genes list file
        core_genes_file = f"{output_dir}/core_genes_list.txt"
        with open(core_genes_file, 'w') as cf:
            for gene, count in core_genes:
                cf.write(f"{gene}\n")
        
        self.logger.info(f"💎 核心基因列表: {core_genes_file}")
