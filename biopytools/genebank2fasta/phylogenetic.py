"""
系统发育分析矩阵构建模块 🌲 | Phylogenetic Analysis Matrix Construction Module
"""

import os
from Bio import SeqIO
from collections import defaultdict

class PhylogeneticMatrixBuilder:
    """系统发育矩阵构建器 🧬 | Phylogenetic Matrix Builder"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def create_phylogenetic_matrix(self):
        """创建系统发育分析矩阵 🌲 | Create phylogenetic analysis matrix"""
        self.logger.info("🌲 创建系统发育分析矩阵...")
        
        # 读取核心基因列表 | Read core genes list
        core_genes_file = f"{self.config.output_dir}/core_genes_list.txt"
        if not os.path.exists(core_genes_file):
            self.logger.error("❌ 未找到核心基因列表，请先运行序列提取")
            return False
        
        with open(core_genes_file, 'r') as f:
            core_genes = [line.strip() for line in f if line.strip()]
        
        self.logger.info(f"📊 使用 {len(core_genes)} 个核心基因构建矩阵")
        
        # 创建输出目录 | Create output directory
        matrix_dir = f"{self.config.output_dir}/phylogenetic_matrix"
        os.makedirs(matrix_dir, exist_ok=True)
        
        # 读取所有样本的序列 | Read sequences from all samples
        sample_sequences = self.collect_sequences(core_genes)
        
        if not sample_sequences:
            self.logger.error("❌ 未找到样本序列数据")
            return False
        
        # 写出超基因组矩阵 | Write supermatrix
        self.write_supermatrix(sample_sequences, core_genes, matrix_dir)
        
        return True
    
    def collect_sequences(self, core_genes):
        """收集序列数据 📚 | Collect sequence data"""
        sample_sequences = defaultdict(dict)
        cds_by_gene_dir = f"{self.config.cds_dir}/by_gene"
        
        for gene in core_genes:
            gene_file = f"{cds_by_gene_dir}/{gene}.fasta"
            if os.path.exists(gene_file):
                for record in SeqIO.parse(gene_file, "fasta"):
                    sample_name = record.id.split()[0]  # 取第一个空格前的部分作为样本名
                    sample_sequences[sample_name][gene] = str(record.seq)
        
        return sample_sequences
    
    def write_supermatrix(self, sample_sequences, core_genes, matrix_dir):
        """写出超基因组矩阵 🧬 | Write supermatrix"""
        supermatrix_file = f"{matrix_dir}/supermatrix.fasta"
        partition_file = f"{matrix_dir}/partitions.txt"
        
        with open(supermatrix_file, 'w') as fasta_out, open(partition_file, 'w') as part_out:
            
            samples = sorted(sample_sequences.keys())
            
            # 写序列 | Write sequences
            for sample in samples:
                concatenated_seq = ""
                
                for gene in core_genes:
                    if gene in sample_sequences[sample]:
                        concatenated_seq += sample_sequences[sample][gene]
                    else:
                        # 用N填充缺失的基因 | Fill missing genes with N
                        avg_length = 1500  # 假设平均基因长度 | Assumed average gene length
                        concatenated_seq += "N" * avg_length
                
                fasta_out.write(f">{sample}\n")
                for i in range(0, len(concatenated_seq), 80):
                    fasta_out.write(f"{concatenated_seq[i:i+80]}\n")
            
            # 写分区文件 | Write partition file
            current_pos = 1
            for gene in core_genes:
                if samples and gene in sample_sequences[samples[0]]:
                    gene_length = len(sample_sequences[samples[0]][gene])
                else:
                    gene_length = 1500  # 默认长度 | Default length
                
                part_out.write(f"DNA, {gene} = {current_pos}-{current_pos + gene_length - 1}\n")
                current_pos += gene_length
        
        self.logger.info(f"🧬 超基因组矩阵已创建: {supermatrix_file}")
        self.logger.info(f"📋 分区文件已创建: {partition_file}")
