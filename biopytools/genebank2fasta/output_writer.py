"""
输出文件写入模块 📝 | Output File Writer Module
"""

import os
from collections import defaultdict

class SequenceWriter:
    """序列写入器 ✍️ | Sequence Writer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def write_sequences(self, all_sample_data):
        """写入所有序列文件 💾 | Write all sequence files"""
        # 重组数据结构 | Reorganize data structure
        sample_cds_sequences = {}
        sample_protein_sequences = {}
        gene_cds_sequences = defaultdict(list)
        gene_protein_sequences = defaultdict(list)
        
        for sample_data in all_sample_data:
            sample_name = sample_data['sample']
            
            # 按样品组织 | Organize by sample
            sample_cds_sequences[sample_name] = sample_data['cds_sequences']
            sample_protein_sequences[sample_name] = sample_data['protein_sequences']
            
            # 按基因组织 | Organize by gene
            for seq_info in sample_data['cds_sequences']:
                gene_name = seq_info['gene']
                gene_cds_sequences[gene_name].append(seq_info)
                gene_protein_sequences[gene_name].append(seq_info)
        
        # 按样品输出 | Output by sample
        if self.config.separate_by_sample:
            self.write_sample_sequences(sample_cds_sequences, self.config.cds_dir, "cds")
            self.write_sample_sequences(sample_protein_sequences, self.config.pep_dir, "pep")
        
        # 按基因输出 | Output by gene
        if self.config.separate_by_gene:
            self.write_gene_sequences(gene_cds_sequences, self.config.cds_dir, "cds")
            self.write_gene_sequences(gene_protein_sequences, self.config.pep_dir, "pep")
    
    def write_sample_sequences(self, sample_sequences_dict, output_base_dir, seq_type):
        """按样品写出序列文件 🔬 | Write sequence files by sample"""
        by_sample_dir = f"{output_base_dir}/by_sample"
        
        for sample_name, seq_list in sample_sequences_dict.items():
            sample_file = f"{by_sample_dir}/{sample_name}_{seq_type}.fasta"
            
            with open(sample_file, 'w') as f:
                for seq_info in seq_list:
                    self.write_fasta_entry(f, seq_info, seq_type)
        
        self.logger.info(f"✅ 已为 {len(sample_sequences_dict)} 个样品创建独立的{seq_type}序列文件")
    
    def write_gene_sequences(self, sequences_dict, output_base_dir, seq_type):
        """按基因写出序列文件 🧬 | Write sequence files by gene"""
        by_gene_dir = f"{output_base_dir}/by_gene"
        
        for gene_name, seq_list in sequences_dict.items():
            gene_file = f"{by_gene_dir}/{gene_name}.fasta"
            
            with open(gene_file, 'w') as f:
                for seq_info in seq_list:
                    self.write_fasta_entry(f, seq_info, seq_type, include_sample=True)
        
        self.logger.info(f"✅ 输出 {len(sequences_dict)} 个基因的{seq_type}序列")
    
    def write_fasta_entry(self, file_handle, seq_info, seq_type, include_sample=False):
        """写入单个FASTA条目 📄 | Write single FASTA entry"""
        gene = seq_info['gene']
        organism = seq_info['organism']
        product = seq_info['product']
        
        if seq_type == "cds":
            sequence = seq_info['cds_seq']
            seq_len = seq_info['cds_length']
            unit = "bp"
        else:  # pep
            sequence = seq_info['protein_seq']
            seq_len = seq_info['protein_length']
            unit = "aa"
        
        if include_sample:
            sample = seq_info['sample']
            header = f">{sample} {organism} {gene} | {product} | {seq_len}{unit}"
        else:
            header = f">{gene} {organism} | {product} | {seq_len}{unit}"
        
        file_handle.write(f"{header}\n")
        
        # 每行80个字符 | 80 characters per line
        for i in range(0, len(sequence), 80):
            file_handle.write(f"{sequence[i:i+80]}\n")
