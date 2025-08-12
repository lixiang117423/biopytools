"""
📁 文件处理模块 | File Processing Module
"""

import os
import re
import shutil
import glob
import gzip
from pathlib import Path
from typing import List

class FileProcessor:
    """📂 文件处理器 | File Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def find_fastq_files(self) -> List[tuple]:
        """🔍 查找FASTQ文件 | Find FASTQ files"""
        self.logger.info(f"🔍 查找FASTQ文件 | Finding FASTQ files")
        self.logger.info(f"📁 输入目录: {self.config.input_dir}")
        self.logger.info(f"🔍 文件模式: {self.config.pattern}")
        
        # 构建完整搜索路径
        input_dir = Path(self.config.input_dir).resolve()
        search_pattern = str(input_dir / self.config.pattern)
        
        self.logger.info(f"🔍 完整搜索路径: {search_pattern}")
        self.logger.info(f"📂 输入目录是否存在: {input_dir.exists()}")
        
        # 列出目录内容用于调试
        if input_dir.exists():
            files_in_dir = list(input_dir.glob("*"))
            self.logger.info(f"📄 目录中的文件: {[f.name for f in files_in_dir[:10]]}")  # 只显示前10个
        
        # 1. 先找到所有匹配模式的R1文件
        r1_files = glob.glob(search_pattern)
        self.logger.info(f"🔍 找到{len(r1_files)}个R1文件")
        
        if not r1_files:
            raise FileNotFoundError(f"📂 未找到匹配的文件 | No files found matching pattern: {search_pattern}")
        
        samples = []
        for r1_file in r1_files:
            # 2. 从R1文件提取样本名（*的部分）
            r1_filename = os.path.basename(r1_file)
            
            # 提取样本名：去掉模式中*后面的部分
            pattern_suffix = self.config.pattern.split('*')[1]  # 例如："_1.clean.fq.gz"
            sample_name = r1_filename.replace(pattern_suffix, "")  # 例如："OV8-105"
            
            self.logger.info(f"📝 提取样本名: {sample_name}")
            
            # 3. 查找对应的R2文件
            r2_pattern_suffix = pattern_suffix.replace('_1.', '_2.').replace('_R1.', '_R2.')
            r2_filename = sample_name + r2_pattern_suffix
            r2_file = input_dir / r2_filename
            
            if r2_file.exists():
                samples.append((sample_name, r1_file, str(r2_file)))
                self.logger.info(f"👥 找到双端样本: {sample_name}")
                self.logger.info(f"  📄 R1: {r1_file}")
                self.logger.info(f"  📄 R2: {r2_file}")
            else:
                samples.append((sample_name, r1_file, None))
                self.logger.info(f"📄 找到单端样本: {sample_name}")
                self.logger.info(f"  📄 R1: {r1_file}")
                self.logger.info(f"  ❌ R2不存在: {r2_file}")
        
        if not samples:
            raise FileNotFoundError("❌ 未找到有效的FASTQ文件 | No valid FASTQ files found")
        
        self.config.samples = samples
        return samples
    
    def decompress_files(self, sample_files: List[str]) -> List[str]:
        """📦 解压缩文件到临时目录 | Decompress files to temporary directory"""
        decompressed_files = []
        
        for file_path in sample_files:
            if file_path is None:
                decompressed_files.append(None)
                continue
                
            if file_path.endswith('.gz'):
                # 解压缩文件 | Decompress file
                base_name = os.path.basename(file_path).replace('.gz', '')
                temp_file = self.config.temp_dir / base_name
                
                self.logger.info(f"📦 解压缩文件 | Decompressing file: {file_path}")
                with gzip.open(file_path, 'rt') as f_in:
                    with open(temp_file, 'w') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                
                decompressed_files.append(str(temp_file))
            else:
                # 复制文件到临时目录 | Copy file to temporary directory
                base_name = os.path.basename(file_path)
                temp_file = self.config.temp_dir / base_name
                shutil.copy2(file_path, temp_file)
                decompressed_files.append(str(temp_file))
        
        return decompressed_files
    
    def prepare_kmer_library(self) -> Path:
        """🧬 准备k-mer库文件 | Prepare k-mer library file"""
        self.logger.info("🧬 准备k-mer库文件 | Preparing k-mer library file")
        
        # 读取原始文件并转换为大写 | Read original file and convert to uppercase
        upper_kmer_file = self.config.temp_dir / "kmers_upper.fasta"
        
        with open(self.config.kmer_lib, 'r') as f_in:
            with open(upper_kmer_file, 'w') as f_out:
                for line in f_in:
                    if line.startswith('>'):
                        f_out.write(line)
                    else:
                        f_out.write(line.upper())
        
        self.logger.info(f"✅ K-mer库文件已准备完成 | K-mer library file prepared: {upper_kmer_file}")
        return upper_kmer_file
