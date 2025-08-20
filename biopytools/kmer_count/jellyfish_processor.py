"""
🐙 Jellyfish处理模块 | Jellyfish Processing Module
"""

import os
import subprocess
from pathlib import Path
from typing import List

class JellyfishProcessor:
    """🐙 Jellyfish处理器 | Jellyfish Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    #     return jf_file
    def count_kmers(self, sample_name: str, fastq_files: List[str]) -> Path:
        """🧮 使用Jellyfish计数k-mer | Count k-mers using Jellyfish"""
        self.logger.info(f"🧮 开始k-mer计数 | Starting k-mer counting for sample: {sample_name}")
        
        # 输出文件 | Output file
        jf_file = self.config.temp_dir / f"{sample_name}.jf"
        
        # 🔥 检查是否有压缩文件，使用generator方式 | Check for compressed files, use generator method
        input_files = [f for f in fastq_files if f is not None]
        has_compressed = any(f.endswith('.gz') for f in input_files)
        
        if has_compressed:
            # 创建generator文件 | Create generator file
            generator_file = self.config.temp_dir / f"{sample_name}_generator.txt"
            
            with open(generator_file, 'w') as f:
                for file_path in input_files:
                    if file_path.endswith('.gz'):
                        f.write(f"zcat {file_path}\n")
                    else:
                        f.write(f"cat {file_path}\n")
            
            self.logger.info(f"📝 创建generator文件 | Created generator file: {generator_file}")
            
            # 构建命令 | Build command
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
            # 普通方式处理未压缩文件 | Normal method for uncompressed files
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
        
        # 执行命令 | Execute command
        self.cmd_runner.run(cmd, f"🧮 K-mer计数 | K-mer counting for {sample_name}")
        
        if not jf_file.exists():
            raise FileNotFoundError(f"❌ Jellyfish输出文件未生成 | Jellyfish output file not generated: {jf_file}")
        
        self.logger.info(f"✅ K-mer计数完成 | K-mer counting completed: {jf_file}")
        return jf_file
    
    def query_kmers(self, sample_name: str, jf_file: Path, kmer_lib: Path) -> Path:
        """🔍 查询k-mer丰度 | Query k-mer abundance"""
        self.logger.info(f"🔍 查询k-mer丰度 | Querying k-mer abundance for sample: {sample_name}")
        
        # 输出文件 | Output file
        count_file = self.config.temp_dir / f"{sample_name}.kmers.count"
        
        # 构建命令 | Build command
        cmd = f"{self.config.jellyfish_path} query -s {kmer_lib} {jf_file} -o {count_file}"
        
        # 执行命令 | Execute command
        self.cmd_runner.run(cmd, f"🔍 查询k-mer丰度 | Querying k-mer abundance for {sample_name}")
        
        if not count_file.exists():
            raise FileNotFoundError(f"❌ 查询输出文件未生成 | Query output file not generated: {count_file}")
        
        self.logger.info(f"✅ K-mer丰度查询完成 | K-mer abundance query completed: {count_file}")
        return count_file
    
    def count_kmers(self, sample_name: str, fastq_files: List[str]) -> Path:
        """🧮 使用Jellyfish计数k-mer（带文件完整性检查）| Count k-mers using Jellyfish (with integrity check)"""
        self.logger.info(f"🧮 开始k-mer计数 | Starting k-mer counting for sample: {sample_name}")
        
        jf_file = self.config.temp_dir / f"{sample_name}.jf"
        input_files = [f for f in fastq_files if f is not None]
        has_compressed = any(f.endswith('.gz') for f in input_files)
        
        if has_compressed:
            generator_file = self.config.temp_dir / f"{sample_name}_generator.txt"
            
            with open(generator_file, 'w') as f:
                for file_path in input_files:
                    if file_path.endswith('.gz'):
                        f.write(f"zcat {file_path}\n")
                    else:
                        f.write(f"cat {file_path}\n")
            
            self.logger.info(f"📝 创建generator文件 | Created generator file: {generator_file}")
            
            cmd_parts = [
                self.config.jellyfish_path,
                "count",
                f"-m {self.config.kmer_size}",
                f"-s {self.config.hash_size}",
                f"-t {self.config.threads}",
                f"-g {generator_file}",
            ]
            
            if self.config.canonical:
                cmd_parts.append("-C")
            
            cmd_parts.extend(["-o", str(jf_file)])
            
        else:
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
        
        try:
            # 尝试执行命令
            self.cmd_runner.run(cmd, f"🧮 K-mer计数 | K-mer counting for {sample_name}")
            
        except subprocess.CalledProcessError as e:
            # 检查是否是SIGPIPE错误
            if "killed by signal 13" in str(e.stderr) or "Some generator commands failed" in str(e.stderr):
                self.logger.warning(f"⚠️ 检测到SIGPIPE错误，可能是文件完整性问题: {sample_name}")
                
                # 检查当前样本的文件完整性
                corrupted_files = []
                for file_path in input_files:
                    # 需要从主程序传入file_processor引用
                    if hasattr(self, 'file_processor'):
                        self.logger.info(f"🔍 检查文件完整性: {os.path.basename(file_path)}")
                        if not self.file_processor.check_file_integrity(file_path):
                            corrupted_files.append(file_path)
                
                if corrupted_files:
                    error_msg = f"💥 样本 {sample_name} 包含损坏的文件，无法处理:"
                    for corrupted_file in corrupted_files:
                        error_msg += f"\n  💔 {os.path.basename(corrupted_file)}"
                    self.logger.error(error_msg)
                    
                    # 抛出特定的异常，让主程序知道这是文件完整性问题
                    raise FileIntegrityError(f"🚫 文件完整性检查失败: {sample_name}")
                else:
                    # 如果文件完整性检查通过，可能是其他问题
                    self.logger.warning(f"🤔 SIGPIPE错误但文件完整性检查通过，可能是其他问题")
                    raise
            else:
                # 其他类型的错误，直接抛出
                raise
        
        if not jf_file.exists():
            raise FileNotFoundError(f"❌ Jellyfish输出文件未生成: {jf_file}")
        
        self.logger.info(f"✅ K-mer计数完成: {jf_file}")
        return jf_file
