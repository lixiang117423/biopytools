"""
K-mer数据处理模块 | K-mer Data Processing Module
"""

import os
import glob
from pathlib import Path
from typing import List, Tuple, Dict
from collections import defaultdict

try:
    import pyfastx
except ImportError:
    pyfastx = None

class FastqFileDetector:
    """FASTQ文件检测器 | FASTQ File Detector"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def extract_sample_name(self, file_name: str) -> str:
        """从文件名提取样本名 | Extract sample name from filename"""
        # 常见的read标识符后缀
        read_suffixes = ['_R1', '_R2', '_1', '_2', '.R1', '.R2', '.1', '.2', '_f1', '_f2', '_r1''_r2']
        
        # 常见的文件扩展名
        file_extensions = [
            '.fastq.gz', '.fq.gz', '.fastq', '.fq',
            '.FASTQ.gz', '.FQ.gz', '.FASTQ', '.FQ',
            '.gz', '.bz2'
        ]
        
        # 首先移除文件扩展名
        base_name = file_name
        for ext in file_extensions:
            if base_name.lower().endswith(ext.lower()):
                base_name = base_name[:-len(ext)]
                break
        
        # 然后移除read标识符
        for suffix in read_suffixes:
            if base_name.upper().endswith(suffix.upper()):
                base_name = base_name[:-len(suffix)]
                break
        
        return base_name if base_name else file_name
    
    def auto_detect_fastq_files(self, fastq_dir: str) -> List[Tuple[str, List[str]]]:
        """自动检测FASTQ文件并生成样本信息 | Auto detect FASTQ files and generate sample info"""
        self.logger.info("自动检测FASTQ文件 | Auto detecting FASTQ files...")
        
        fastq_dir = Path(fastq_dir)
        if not fastq_dir.exists():
            raise FileNotFoundError(f"FASTQ目录不存在 | FASTQ directory does not exist: {fastq_dir}")
        
        # 支持的文件扩展名
        patterns = [
            "*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq",
            "*.FASTQ.gz", "*.FQ.gz", "*.FASTQ", "*.FQ"
        ]
        
        # 查找所有FASTQ文件
        all_files = []
        for pattern in patterns:
            all_files.extend(fastq_dir.glob(pattern))
        
        if not all_files:
            self.logger.error(f"在{fastq_dir}中未找到FASTQ文件 | No FASTQ files found in {fastq_dir}")
            raise FileNotFoundError(f"在{fastq_dir}中未找到FASTQ文件 | No FASTQ files found in {fastq_dir}")
        
        self.logger.info(f"找到{len(all_files)}个FASTQ文件 | Found {len(all_files)} FASTQ files")
        
        # 分组样本文件
        samples = defaultdict(list)
        
        for file_path in all_files:
            file_name = file_path.name
            sample_name = self.extract_sample_name(file_name)
            
            if sample_name:
                samples[sample_name].append(str(file_path))
        
        # 整理样本信息
        sample_info = []
        for sample_name, files in samples.items():
            files.sort()  # 确保R1在R2前面
            sample_info.append((sample_name, files))
        
        sample_info.sort(key=lambda x: x[0])
        
        self.logger.info(f"检测到{len(sample_info)}个样本 | Detected {len(sample_info)} samples")
        if len(sample_info) > 100:
            self.logger.info(f"样本示例 | Sample examples: {[x[0] for x in sample_info[:5]]}... (显示前5个 | showing first 5)")
        
        return sample_info

class FOFFileGenerator:
    """FOF文件生成器 | FOF File Generator"""
    
    def __init__(self, config, logger, fastq_detector: FastqFileDetector):
        self.config = config
        self.logger = logger
        self.fastq_detector = fastq_detector
    
    def create_fof_file(self) -> str:
        """创建FOF文件 | Create FOF file"""
        self.logger.info("创建FOF文件 | Creating FOF file...")
        
        sample_info = self.fastq_detector.auto_detect_fastq_files(self.config.fastq_dir)
        
        with open(self.config.fof_file, 'w') as fof:
            for sample_name, files in sample_info:
                if len(files) == 1:
                    fof.write(f"{sample_name}:{files[0]}\n")
                else:
                    file_str = ";".join(files)
                    fof.write(f"{sample_name}:{file_str}\n")
        
        self.logger.info(f"FOF文件已创建 | FOF file created: {self.config.fof_file}")
        self.logger.info(f"包含 {len(sample_info)} 个样本 | Contains {len(sample_info)} samples")
        return str(self.config.fof_file)

class GeneKmerExtractor:
    """基因k-mer提取器 | Gene K-mer Extractor"""
    
    def __init__(self, config, logger, seq_utils):
        self.config = config
        self.logger = logger
        self.seq_utils = seq_utils
    
    def extract_gene_kmers(self) -> Dict[str, List[Tuple[str, int]]]:
        """从基因FASTA文件提取k-mer | Extract k-mers from gene FASTA file"""
        if not pyfastx:
            raise ImportError("pyfastx未安装，请安装: pip install pyfastx | pyfastx not installed, please install: pip install pyfastx")
        
        self.logger.info("从基因序列提取k-mer | Extracting k-mers from gene sequences...")
        
        gene_kmers = {}
        all_kmers = []
        kmer_positions = []
        kmer_id = 1
        
        try:
            for gene_name, seq in pyfastx.Fasta(self.config.gene_fasta, build_index=False):
                seq = seq.upper()
                gene_kmer_list = []
                seen_kmers = set()
                
                self.logger.info(f"处理基因 | Processing gene: {gene_name} (长度 | length: {len(seq)} bp)")
                
                if len(seq) < self.config.kmer_size:
                    self.logger.warning(f"基因 {gene_name} 长度不足，跳过 | Gene {gene_name} too short, skipping")
                    continue
                
                for i in range(len(seq) - self.config.kmer_size + 1):
                    kmer = seq[i:i+self.config.kmer_size]
                    if 'N' not in kmer:
                        canonical_kmer = self.seq_utils.get_canonical_kmer(kmer)
                        
                        if canonical_kmer not in seen_kmers:
                            seen_kmers.add(canonical_kmer)
                            gene_kmer_list.append((canonical_kmer, i + 1))
                            all_kmers.append(canonical_kmer)
                            kmer_positions.append({
                                'kmer_id': kmer_id,
                                'gene': gene_name,
                                'kmer': canonical_kmer,
                                'position': i + 1
                            })
                            kmer_id += 1
                
                gene_kmers[gene_name] = gene_kmer_list
                self.logger.info(f"基因 {gene_name}: 提取到 {len(gene_kmer_list)} 个唯一k-mer | Gene {gene_name}: extracted {len(gene_kmer_list)} unique k-mers")
        
        except Exception as e:
            self.logger.error(f"提取基因k-mer失败 | Failed to extract gene k-mers: {e}")
            raise
        
        # 保存k-mer列表
        with open(self.config.gene_kmer_file, 'w') as f:
            f.write("kmer\n")
            for kmer in all_kmers:
                f.write(f"{kmer}\n")
        
        # 保存位置信息
        with open(self.config.gene_kmer_pos_file, 'w') as f:
            f.write("kmer_id\tgene\tkmer\tposition\n")
            for pos_info in kmer_positions:
                f.write(f"{pos_info['kmer_id']}\t{pos_info['gene']}\t{pos_info['kmer']}\t{pos_info['position']}\n")
        
        self.logger.info(f"总计提取到 {len(all_kmers)} 个唯一k-mer | Total extracted {len(all_kmers)} unique k-mers")
        return gene_kmers
