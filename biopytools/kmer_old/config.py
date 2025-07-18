"""
K-mer分析配置管理模块 | K-mer Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class KmerConfig:
    """K-mer分析配置类 | K-mer Analysis Configuration Class"""
    
    # 必需文件 | Required files
    gene_fasta: str
    fastq_dir: str
    output_dir: str
    
    # K-mer参数 | K-mer parameters
    kmer_size: int = 51
    hard_min: int = 2
    threads: int = 32
    
    # 分析参数 | Analysis parameters
    n_components: int = 5
    project_name: Optional[str] = None
    
    # 流程控制 | Process control
    skip_build: bool = False
    run_haplotype: bool = False
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.gene_fasta = os.path.normpath(os.path.abspath(self.gene_fasta))
        self.fastq_dir = os.path.normpath(os.path.abspath(self.fastq_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        
        # 设置项目名称 | Set project name
        if not self.project_name:
            self.project_name = self.output_path.name or "kmer_db_analysis"
        
        # 设置各种路径 | Set various paths
        self.kmtricks_run_dir = self.output_path / f"{self.project_name}.k{self.kmer_size}"
        self.kmer_matrix_file = self.output_path / f"{self.project_name}.kmer.matrix.txt.gz"
        self.header_file = self.output_path / "header.txt"
        self.rocksdb_dir = self.output_path / f"{self.project_name}.kmer.matrix.rocksdb"
        
        # 输入输出文件 | Input/output files
        self.fof_file = self.output_path / "samples.fof"
        self.gene_kmer_file = self.output_path / "gene_kmers.txt"
        self.gene_kmer_pos_file = self.output_path / "gene_kmers_positions.txt"
        self.query_result_file = self.output_path / "query_results.txt"
        self.kmer_matrix_final = self.output_path / f"{self.project_name}_kmer_matrix.tsv"
        
        # 分析结果 | Analysis results
        self.haplotype_results = self.output_path / "haplotype_results.tsv"
        self.heatmap_file = self.output_path / "haplotype_heatmap.png"
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查必需文件 | Check required files
        if not os.path.exists(self.gene_fasta):
            errors.append(f"基因FASTA文件不存在 | Gene FASTA file does not exist: {self.gene_fasta}")
        
        if not os.path.exists(self.fastq_dir):
            errors.append(f"FASTQ目录不存在 | FASTQ directory does not exist: {self.fastq_dir}")
        
        # 检查参数范围 | Check parameter ranges
        if self.kmer_size <= 0:
            errors.append(f"k-mer大小必须为正整数 | k-mer size must be positive integer: {self.kmer_size}")
        
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Thread count must be positive integer: {self.threads}")
        
        if self.hard_min < 0:
            errors.append(f"最小频次不能为负数 | Hard min cannot be negative: {self.hard_min}")
        
        if self.n_components <= 0:
            errors.append(f"聚类数必须为正整数 | Number of components must be positive integer: {self.n_components}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True
