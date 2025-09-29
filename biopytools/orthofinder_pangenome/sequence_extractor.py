# # ===== FILE: orthofinder_pangenome/sequence_extractor.py =====
# """
# 序列提取模块 | Sequence Extraction Module
# """

# import os
# from pathlib import Path
# from typing import Dict, List, Tuple
# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq

# class SequenceExtractor:
#     """序列提取器 | Sequence Extractor"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
#         self.genome_sequences = {}  # 基因组序列字典
        
#     def load_genome_sequences(self):
#         """加载所有基因组序列 | Load all genome sequences"""
#         self.logger.info("加载基因组序列文件 | Loading genome sequence files")
        
#         input_path = Path(self.config.input_dir)
#         fasta_extensions = ['*.fa', '*.faa', '*.fas', '*.fasta', '*.pep', '*.protein']
        
#         fasta_files = []
#         for ext in fasta_extensions:
#             fasta_files.extend(input_path.glob(ext))
        
#         for fasta_file in fasta_files:
#             genome_name = fasta_file.stem
#             self.logger.info(f"  加载 {genome_name} 序列")
            
#             sequences = {}
#             try:
#                 for record in SeqIO.parse(fasta_file, "fasta"):
#                     sequences[record.id] = record
#                 self.genome_sequences[genome_name] = sequences
#             except Exception as e:
#                 self.logger.error(f"加载 {fasta_file} 失败: {e}")
        
#         total_genomes = len(self.genome_sequences)
#         total_sequences = sum(len(seqs) for seqs in self.genome_sequences.values())
#         self.logger.info(f"序列加载完成 | Sequence loading completed: {total_genomes} genomes, {total_sequences} sequences")
    
#     def extract_single_copy_sequences(self, single_copy_genes: List[Tuple], output_dir: Path):
#         """提取单拷贝基因序列 | Extract single copy gene sequences"""
#         if not self.config.extract_sequences:
#             return
        
#         self.logger.info("提取单拷贝基因序列 | Extracting single copy gene sequences")
        
#         # 按同源群组织数据
#         orthogroup_genes = {}
#         for og_id, gene_ids, genome_names in single_copy_genes:
#             if og_id not in orthogroup_genes:
#                 orthogroup_genes[og_id] = []
#             for gene_id, genome_name in zip(gene_ids, genome_names):
#                 orthogroup_genes[og_id].append((gene_id, genome_name))
        
#         # 创建输出目录
#         if self.config.single_copy_output_format in ['by_orthogroup', 'both']:
#             og_dir = output_dir / "single_copy_by_orthogroup"
#             og_dir.mkdir(exist_ok=True)
#             self._extract_by_orthogroup(orthogroup_genes, og_dir)
        
#         if self.config.single_copy_output_format in ['by_genome', 'both']:
#             genome_dir = output_dir / "single_copy_by_genome"
#             genome_dir.mkdir(exist_ok=True)
#             self._extract_by_genome(orthogroup_genes, genome_dir)
        
#         # 生成汇总文件
#         self._generate_summary_file(orthogroup_genes, output_dir)
        
#     def _extract_by_orthogroup(self, orthogroup_genes: Dict, output_dir: Path):
#         """按同源群提取序列 | Extract sequences by orthogroup"""
#         self.logger.info("按同源群提取序列 | Extracting sequences by orthogroup")
        
#         for og_id, gene_info_list in orthogroup_genes.items():
#             sequences = []
            
#             for gene_id, genome_name in gene_info_list:
#                 if genome_name in self.genome_sequences:
#                     if gene_id in self.genome_sequences[genome_name]:
#                         seq_record = self.genome_sequences[genome_name][gene_id]
#                         # 重命名序列ID包含基因组信息
#                         new_record = SeqRecord(
#                             seq_record.seq,
#                             id=f"{genome_name}_{gene_id}",
#                             description=f"{og_id} | {genome_name} | {seq_record.description}"
#                         )
#                         sequences.append(new_record)
            
#             if sequences:
#                 output_file = output_dir / f"{og_id}.fasta"
#                 SeqIO.write(sequences, output_file, "fasta")
        
#         self.logger.info(f"按同源群序列提取完成 | Orthogroup-wise extraction completed: {len(orthogroup_genes)} files")
    
#     def _extract_by_genome(self, orthogroup_genes: Dict, output_dir: Path):
#         """按基因组提取序列 | Extract sequences by genome"""
#         self.logger.info("按基因组提取序列 | Extracting sequences by genome")
        
#         genome_sequences = {}
        
#         for og_id, gene_info_list in orthogroup_genes.items():
#             for gene_id, genome_name in gene_info_list:
#                 if genome_name not in genome_sequences:
#                     genome_sequences[genome_name] = []
                
#                 if genome_name in self.genome_sequences:
#                     if gene_id in self.genome_sequences[genome_name]:
#                         seq_record = self.genome_sequences[genome_name][gene_id]
#                         new_record = SeqRecord(
#                             seq_record.seq,
#                             id=f"{og_id}_{gene_id}",
#                             description=f"{og_id} | {seq_record.description}"
#                         )
#                         genome_sequences[genome_name].append(new_record)
        
#         for genome_name, sequences in genome_sequences.items():
#             if sequences:
#                 output_file = output_dir / f"{genome_name}_single_copy.fasta"
#                 SeqIO.write(sequences, output_file, "fasta")
        
#         self.logger.info(f"按基因组序列提取完成 | Genome-wise extraction completed: {len(genome_sequences)} files")
    
#     def _generate_summary_file(self, orthogroup_genes: Dict, output_dir: Path):
#         """生成汇总文件 | Generate summary file"""
#         summary_file = output_dir / "single_copy_genes_summary.txt"
        
#         with open(summary_file, 'w') as f:
#             f.write("单拷贝基因汇总 | Single Copy Genes Summary\n")
#             f.write("=" * 50 + "\n\n")
#             f.write(f"同源群数量 | Number of orthogroups: {len(orthogroup_genes)}\n")
#             f.write(f"总基因数量 | Total genes: {sum(len(genes) for genes in orthogroup_genes.values())}\n\n")
            
#             f.write("Orthogroup_ID\tGene_Count\tGenomes\n")
#             for og_id, gene_info_list in orthogroup_genes.items():
#                 genomes = [genome for _, genome in gene_info_list]
#                 f.write(f"{og_id}\t{len(gene_info_list)}\t{','.join(genomes)}\n")
        
#         self.logger.info(f"单拷贝基因汇总文件已生成 | Single copy genes summary generated: {summary_file}")

# # ===== END FILE =====

"""
🧬 序列提取模块 | Sequence Extraction Module
"""

import os
from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class SequenceExtractor:
    """🔍 序列提取器 | Sequence Extractor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.genome_sequences = {}  # 🧬 基因组序列字典
        
    def load_genome_sequences(self):
        """📥 加载所有基因组序列 | Load all genome sequences"""
        self.logger.info("📥 加载基因组序列文件 | Loading genome sequence files")
        
        input_path = Path(self.config.input_dir)
        fasta_extensions = ['*.fa', '*.faa', '*.fas', '*.fasta', '*.pep', '*.protein']
        
        fasta_files = []
        for ext in fasta_extensions:
            fasta_files.extend(input_path.glob(ext))
        
        for fasta_file in fasta_files:
            genome_name = fasta_file.stem
            self.logger.info(f"  📄 加载 {genome_name} 序列")
            
            sequences = {}
            try:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    sequences[record.id] = record
                self.genome_sequences[genome_name] = sequences
            except Exception as e:
                self.logger.error(f"❌ 加载 {fasta_file} 失败: {e}")
        
        total_genomes = len(self.genome_sequences)
        total_sequences = sum(len(seqs) for seqs in self.genome_sequences.values())
        self.logger.info(f"✅ 序列加载完成 | Sequence loading completed: {total_genomes} genomes, {total_sequences} sequences")
    
    def extract_single_copy_sequences(self, single_copy_genes: List[Tuple], output_dir: Path):
        """🔵 提取单拷贝基因序列 | Extract single copy gene sequences"""
        if not self.config.extract_sequences:
            return
        
        self.logger.info("🔵 提取单拷贝基因序列 | Extracting single copy gene sequences")
        
        # 📋 按同源群组织数据
        orthogroup_genes = {}
        for og_id, gene_ids, genome_names in single_copy_genes:
            if og_id not in orthogroup_genes:
                orthogroup_genes[og_id] = []
            for gene_id, genome_name in zip(gene_ids, genome_names):
                orthogroup_genes[og_id].append((gene_id, genome_name))
        
        # 📁 创建输出目录
        if self.config.single_copy_output_format in ['by_orthogroup', 'both']:
            og_dir = output_dir / "single_copy_by_orthogroup"
            og_dir.mkdir(exist_ok=True)
            self._extract_by_orthogroup(orthogroup_genes, og_dir)
        
        if self.config.single_copy_output_format in ['by_genome', 'both']:
            genome_dir = output_dir / "single_copy_by_genome"
            genome_dir.mkdir(exist_ok=True)
            self._extract_by_genome(orthogroup_genes, genome_dir)
        
        # 📝 生成汇总文件
        self._generate_summary_file(orthogroup_genes, output_dir)
        
    def _extract_by_orthogroup(self, orthogroup_genes: Dict, output_dir: Path):
        """📦 按同源群提取序列 | Extract sequences by orthogroup"""
        self.logger.info("📦 按同源群提取序列 | Extracting sequences by orthogroup")
        
        for og_id, gene_info_list in orthogroup_genes.items():
            sequences = []
            
            for gene_id, genome_name in gene_info_list:
                if genome_name in self.genome_sequences:
                    if gene_id in self.genome_sequences[genome_name]:
                        seq_record = self.genome_sequences[genome_name][gene_id]
                        # 🏷️ 重命名序列ID包含基因组信息
                        new_record = SeqRecord(
                            seq_record.seq,
                            id=f"{genome_name}_{gene_id}",
                            description=f"{og_id} | {genome_name} | {seq_record.description}"
                        )
                        sequences.append(new_record)
            
            if sequences:
                output_file = output_dir / f"{og_id}.fasta"
                SeqIO.write(sequences, output_file, "fasta")
        
        self.logger.info(f"✅ 按同源群序列提取完成 | Orthogroup-wise extraction completed: {len(orthogroup_genes)} files")
    
    def _extract_by_genome(self, orthogroup_genes: Dict, output_dir: Path):
        """🧬 按基因组提取序列 | Extract sequences by genome"""
        self.logger.info("🧬 按基因组提取序列 | Extracting sequences by genome")
        
        genome_sequences = {}
        
        for og_id, gene_info_list in orthogroup_genes.items():
            for gene_id, genome_name in gene_info_list:
                if genome_name not in genome_sequences:
                    genome_sequences[genome_name] = []
                
                if genome_name in self.genome_sequences:
                    if gene_id in self.genome_sequences[genome_name]:
                        seq_record = self.genome_sequences[genome_name][gene_id]
                        new_record = SeqRecord(
                            seq_record.seq,
                            id=f"{og_id}_{gene_id}",
                            description=f"{og_id} | {seq_record.description}"
                        )
                        genome_sequences[genome_name].append(new_record)
        
        for genome_name, sequences in genome_sequences.items():
            if sequences:
                output_file = output_dir / f"{genome_name}_single_copy.fasta"
                SeqIO.write(sequences, output_file, "fasta")
        
        self.logger.info(f"✅ 按基因组序列提取完成 | Genome-wise extraction completed: {len(genome_sequences)} files")
    
    def _generate_summary_file(self, orthogroup_genes: Dict, output_dir: Path):
        """📋 生成汇总文件 | Generate summary file"""
        summary_file = output_dir / "single_copy_genes_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("🔵 单拷贝基因汇总 | Single Copy Genes Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"📦 同源群数量 | Number of orthogroups: {len(orthogroup_genes)}\n")
            f.write(f"🧬 总基因数量 | Total genes: {sum(len(genes) for genes in orthogroup_genes.values())}\n\n")
            
            f.write("Orthogroup_ID\tGene_Count\tGenomes\n")
            for og_id, gene_info_list in orthogroup_genes.items():
                genomes = [genome for _, genome in gene_info_list]
                f.write(f"{og_id}\t{len(gene_info_list)}\t{','.join(genomes)}\n")
        
        self.logger.info(f"✅ 单拷贝基因汇总文件已生成 | Single copy genes summary generated: {summary_file}")

# ===== END FILE =====