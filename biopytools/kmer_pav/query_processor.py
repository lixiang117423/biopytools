"""
查询样本处理模块 | Query Sample Processing Module
"""

import os
from pathlib import Path
from typing import List, Dict
from .utils import KMCRunner, FileDetector, FastaSequenceSplitter, find_sequence_files, cleanup_files

class QueryProcessor:
    """查询样本处理器 | Query Sample Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.kmc_runner = KMCRunner(config, logger)
        self.file_detector = FileDetector(logger)
        self.fasta_splitter = FastaSequenceSplitter(logger, self.config.output_path)
    
    def get_query_files(self) -> List[Dict[str, str]]:
        """获取查询文件列表 | Get query file list"""
        files = find_sequence_files(self.config.query_input, self.config.query_pattern)
        
        if not files:
            raise ValueError(f"未找到查询文件 | No query files found in: {self.config.query_input}")
        
        query_files = []
        for file_path in files:
            file_format = self.file_detector.detect_format(file_path)
            
            query_files.append({
                'name': Path(file_path).stem,
                'file': file_path,
                'format': file_format
            })
        
        self.logger.info(f"找到 {len(query_files)} 个查询文件 | Found {len(query_files)} query files")
        return query_files
    
    def prepare_samples(self, query_files: List[Dict[str, str]]) -> List[Dict[str, str]]:
        """准备样本列表 | Prepare sample list"""
        self.logger.info("准备样本列表 | Preparing sample list")
        
        samples = []
        
        for query_file in query_files:
            file_format = query_file['format']
            file_path = query_file['file']
            file_name = query_file['name']
            
            if file_format == 'fastq':
                # FASTQ文件：整个文件作为一个样本 | FASTQ file: entire file as one sample
                samples.append({
                    'sample_name': file_name,
                    'file_path': file_path,
                    'file_format': file_format,
                    'sample_type': 'fastq_file'
                })
                self.logger.info(f"FASTQ样本: {file_name}")
                
            elif file_format == 'fasta':
                # FASTA文件：每条序列作为一个样本 | FASTA file: each sequence as one sample
                sequence_files = self.fasta_splitter.split_fasta_file(file_path)
                
                for seq_info in sequence_files:
                    samples.append({
                        'sample_name': seq_info['sample_name'],
                        'file_path': seq_info['file_path'],
                        'file_format': 'fasta',
                        'sample_type': 'fasta_sequence',
                        'original_file': seq_info['original_file'],
                        'sequence_id': seq_info['sequence_id'],
                        'sequence_length': seq_info['sequence_length']
                    })
                
                self.logger.info(f"FASTA文件 {file_name}: 分割为 {len(sequence_files)} 个序列样本")
        
        self.logger.info(f"总共准备了 {len(samples)} 个样本 | Total {len(samples)} samples prepared")
        return samples
    
    def process_sample(self, sample_info: Dict[str, str], database_name: str) -> str:
        """处理单个样本 | Process single sample"""
        sample_name = sample_info['sample_name']
        file_path = sample_info['file_path']
        file_format = sample_info['file_format']
        
        self.logger.info(f"处理样本: {sample_name}")
        
        # 1. 构建样本k-mer数据库 | Build sample k-mer database
        sample_db_name = f"sample_{sample_name.replace('.', '_').replace('-', '_')}"
        format_flag = "-fm" if file_format == "fasta" else "-fq"
        
        # 使用绝对路径 | Use absolute paths
        sample_db_path = str(self.config.output_path / sample_db_name)
        
        cmd = [
            self.config.kmc_path,
            f"-k{self.config.kmer_size}",
            f"-ci{self.config.min_count}",
            f"-cx{self.config.max_count}",
            f"-t{self.config.threads}",
            f"-m{self.config.kmc_memory_gb}",
            format_flag,  # 添加格式标志 | Add format flag
            file_path,
            sample_db_path,  # 使用绝对路径 | Use absolute path
            str(self.config.tmp_path)
        ]
        
        if not self.kmc_runner.run_command(cmd, f"构建样本k-mer数据库: {sample_name}"):
            return None
        
        # 2. 与数据库k-mer取交集 | Intersect with database k-mers
        intersect_db_name = f"intersect_{sample_name.replace('.', '_').replace('-', '_')}"
        
        # 使用绝对路径 | Use absolute paths
        database_path = str(self.config.output_path / database_name)
        intersect_db_path = str(self.config.output_path / intersect_db_name)
        
        cmd = [
            self.config.kmc_tools_path,
            "simple",
            database_path,      # 使用绝对路径 | Use absolute path
            sample_db_path,     # 使用绝对路径 | Use absolute path
            "intersect",
            intersect_db_path   # 使用绝对路径 | Use absolute path
        ]
        
        if not self.kmc_runner.run_command(cmd, f"计算k-mer交集: {sample_name}"):
            return None
        
        # 3. 导出k-mer计数 | Export k-mer counts
        count_file = f"{sample_name.replace('.', '_').replace('-', '_')}_counts.txt"
        count_file_path = str(self.config.output_path / count_file)
        
        cmd = [
            self.config.kmc_tools_path,
            "transform",
            intersect_db_path,  # 使用绝对路径 | Use absolute path
            "dump",
            "-s",
            count_file_path     # 使用绝对路径 | Use absolute path
        ]
        
        if self.kmc_runner.run_command(cmd, f"导出k-mer计数: {sample_name}"):
            # 清理中间文件 | Clean up intermediate files
            if not self.config.keep_intermediate:
                cleanup_files([
                    f"{sample_db_name}.kmc*",
                    f"{intersect_db_name}.kmc*"
                ], self.config.output_path, self.logger)
            
            return count_file
        else:
            return None
    
    def process_all_samples(self, samples: List[Dict[str, str]], database_name: str) -> List[str]:
        """处理所有样本 | Process all samples"""
        self.logger.info("=" * 60)
        self.logger.info("阶段2: 处理查询样本 | Phase 2: Processing query samples")
        self.logger.info("=" * 60)
        
        count_files = []
        sample_names = []
        
        for i, sample_info in enumerate(samples, 1):
            self.logger.info(f"进度: {i}/{len(samples)}")
            
            count_file = self.process_sample(sample_info, database_name)
            if count_file:
                count_files.append(count_file)
                sample_names.append(sample_info['sample_name'])
            else:
                self.logger.error(f"样本处理失败: {sample_info['sample_name']}")
                count_files.append(None)
                sample_names.append(sample_info['sample_name'])
        
        self.logger.info(f"样本处理完成，成功处理 {len([f for f in count_files if f])} 个样本")
        return count_files, sample_names
