"""
Jellyfish处理模块|Jellyfish Processing Module
"""

import os
import subprocess
from pathlib import Path
from typing import List
from modelscope.hub.errors import FileIntegrityError

class JellyfishProcessor:
    """Jellyfish处理器|Jellyfish Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    #     return jf_file
    def count_kmers(self, sample_name: str, fastq_files: List[str]) -> Path:
        """使用Jellyfish计数k-mer|Count k-mers using Jellyfish（支持断点续传）"""
        self.logger.info(f"开始k-mer计数|Starting k-mer counting for sample: {sample_name}")

        # 输出文件：使用持久化cache目录而非temp_dir|Output file: use persistent cache directory
        cache_dir = self.config.output_dir / "jellyfish_cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        jf_file = cache_dir / f"{sample_name}.jf"

        # 检查是否已存在（断点续传）|Check if already exists for checkpoint resumption
        if jf_file.exists() and jf_file.stat().st_size > 0:
            self.logger.info(f"发现已有jellyfish计数文件，跳过|Found existing jellyfish file, skipping: {jf_file}")
            return jf_file

        # 检查是否有压缩文件，使用generator方式|Check for compressed files, use generator method
        input_files = [f for f in fastq_files if f is not None]
        has_compressed = any(f.endswith('.gz') for f in input_files)

        # generator文件也放在cache目录|Generator file also in cache directory
        generator_file = cache_dir / f"{sample_name}_generator.txt"

        if has_compressed:
            # 创建generator文件|Create generator file
            with open(generator_file, 'w') as f:
                for file_path in input_files:
                    if file_path.endswith('.gz'):
                        f.write(f"zcat {file_path}\n")
                    else:
                        f.write(f"cat {file_path}\n")
            
            self.logger.info(f" 创建generator文件|Created generator file: {generator_file}")
            
            # 构建命令|Build command
            cmd_parts = [
                self.config.jellyfish_path,
                "count",
                f"-m {self.config.kmer_size}",
                f"-s {self.config.hash_size}",
                f"-t {self.config.threads}",
                f"-g {generator_file}",  # 使用generator文件
            ]
            
            if self.config.canonical:
                cmd_parts.append("-C")
            
            cmd_parts.extend(["-o", str(jf_file)])
            
        else:
            # 普通方式处理未压缩文件|Normal method for uncompressed files
            cmd_parts = [
                self.config.jellyfish_path,
                "count",
                f"-m {self.config.kmer_size}",
                f"-s {self.config.hash_size}",
                f"-t {self.config.threads}",
            ]
            
            if self.config.canonical:
                cmd_parts.append("-C")
            
            cmd_parts.extend(input_files)
            cmd_parts.extend(["-o", str(jf_file)])
        
        cmd = " ".join(cmd_parts)
        
        # 执行命令|Execute command
        self.cmd_runner.run(cmd, f" K-mer计数|K-mer counting for {sample_name}")
        
        if not jf_file.exists():
            raise FileNotFoundError(f" Jellyfish输出文件未生成|Jellyfish output file not generated: {jf_file}")
        
        self.logger.info(f" K-mer计数完成|K-mer counting completed: {jf_file}")
        return jf_file
    
    def query_kmers(self, sample_name: str, jf_file: Path, kmer_lib: Path) -> Path:
        """ 查询k-mer丰度|Query k-mer abundance（支持断点续传）"""
        self.logger.info(f" 查询k-mer丰度|Querying k-mer abundance for sample: {sample_name}")

        # 输出文件：使用持久化cache目录|Output file: use persistent cache directory
        cache_dir = self.config.output_dir / "jellyfish_cache"
        count_file = cache_dir / f"{sample_name}.kmers.count"

        # 检查是否已存在（断点续传）|Check if already exists for checkpoint resumption
        if count_file.exists() and count_file.stat().st_size > 0:
            self.logger.info(f"发现已有k-mer计数文件，跳过|Found existing count file, skipping: {count_file}")
            return count_file

        # 构建命令|Build command
        cmd = f"{self.config.jellyfish_path} query -s {kmer_lib} {jf_file} -o {count_file}"

        # 执行命令|Execute command
        self.cmd_runner.run(cmd, f" 查询k-mer丰度|Querying k-mer abundance for {sample_name}")

        if not count_file.exists():
            raise FileNotFoundError(f" 查询输出文件未生成|Query output file not generated: {count_file}")
        
        self.logger.info(f" K-mer丰度查询完成|K-mer abundance query completed: {count_file}")
        return count_file
    
