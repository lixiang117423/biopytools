"""
🧬 K-mer提取核心模块 | K-mer Extraction Core Module
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple
from .utils import CommandRunner

class UnikmrExtractor:
    """🧬 基于Unikmer的K-mer提取器 | Unikmer-based K-mer Extractor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def extract_kmers_from_fasta(self, fasta_files: List[str]) -> str:
        """从FASTA文件提取k-mer | Extract k-mers from FASTA files"""
        self.logger.info(f"🧬 从 {len(fasta_files)} 个FASTA文件提取k-mer | Extracting k-mers from {len(fasta_files)} FASTA files")
        
        # 准备输出文件 | Prepare output files
        output_prefix = self.config.output_path / self.config.base_name
        unik_file = f"{output_prefix}.unik"
        
        # 构建unikmer命令 | Build unikmer command
        cmd_parts = [
            self.config.unikmer_path,
            "count",
            f"-k {self.config.kmer_length}",
            f"-j {self.config.threads}",
            "-o", str(output_prefix)
        ]
        
        if self.config.canonical:
            cmd_parts.append("--canonical")
        
        if not self.config.compress_output:
            cmd_parts.append("--no-compress")
        
        # 添加输入文件 | Add input files
        cmd_parts.extend(fasta_files)
        
        cmd = " ".join(cmd_parts)
        
        # 执行提取 | Execute extraction
        if not self.cmd_runner.run(cmd, f"🧬 提取k-mer | Extracting k-mers (k={self.config.kmer_length})"):
            raise RuntimeError("❌ K-mer提取失败 | K-mer extraction failed")
        
        return unik_file
    
    def get_kmer_statistics(self, unik_file: str) -> Dict[str, int]:
        """获取k-mer统计信息 | Get k-mer statistics"""
        info_cmd = f"{self.config.unikmer_path} info {unik_file}"
        
        try:
            result = subprocess.run(info_cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                stats = {'total_kmers': 0}
                
                for line in lines:
                    if 'k-mers' in line.lower():
                        import re
                        numbers = re.findall(r'\d+', line)
                        if numbers:
                            stats['total_kmers'] = int(numbers[0])
                        break
                
                self.logger.info(f"📊 K-mer统计 | K-mer statistics: {stats['total_kmers']} k-mers")
                return stats
            else:
                self.logger.warning(f"⚠️ 无法获取k-mer统计 | Cannot get k-mer statistics: {result.stderr}")
                return {'total_kmers': 0}
        except Exception as e:
            self.logger.warning(f"⚠️ 获取统计信息时出错 | Error getting statistics: {e}")
            return {'total_kmers': 0}

class JellyfishExtractor:
    """🐟 基于Jellyfish的K-mer提取器 | Jellyfish-based K-mer Extractor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def extract_kmers_from_fastq(self, merged_r1: str, merged_r2: str = None) -> str:
        """从合并的FASTQ文件提取k-mer | Extract k-mers from merged FASTQ files"""
        self.logger.info(f"🐟 使用Jellyfish从FASTQ文件提取k-mer | Extracting k-mers from FASTQ files using Jellyfish")
        
        # 准备输出文件 | Prepare output files
        jf_file = str(self.config.output_path / f"{self.config.base_name}.jf")
        
        # 构建jellyfish count命令 | Build jellyfish count command
        cmd_parts = [
            self.config.jellyfish_path,
            "count",
            f"-m {self.config.kmer_length}",
            f"-s {self.config.jellyfish_hash_size}",
            f"-t {self.config.threads}",
            "-o", jf_file
        ]
        
        if self.config.canonical:
            cmd_parts.append("-C")
        
        # 添加输入文件 | Add input files
        cmd_parts.append(merged_r1)
        if merged_r2:
            cmd_parts.append(merged_r2)
            self.logger.info(f"📄 处理双端测序数据 | Processing paired-end data")
        else:
            self.logger.info(f"📄 处理单端测序数据 | Processing single-end data")
        
        cmd = " ".join(cmd_parts)
        
        # 执行提取 | Execute extraction
        if not self.cmd_runner.run(cmd, f"🐟 Jellyfish提取k-mer | Jellyfish k-mer extraction (k={self.config.kmer_length})"):
            raise RuntimeError("❌ Jellyfish K-mer提取失败 | Jellyfish K-mer extraction failed")
        
        return jf_file
    
    def get_kmer_statistics(self, jf_file: str) -> Dict[str, int]:
        """获取k-mer统计信息 | Get k-mer statistics"""
        stats_cmd = f"{self.config.jellyfish_path} stats {jf_file}"
        
        try:
            result = subprocess.run(stats_cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                stats = {'total_kmers': 0, 'unique_kmers': 0}
                
                for line in lines:
                    if 'Unique:' in line:
                        import re
                        numbers = re.findall(r'\d+', line)
                        if numbers:
                            stats['unique_kmers'] = int(numbers[0])
                    elif 'Distinct:' in line:
                        import re
                        numbers = re.findall(r'\d+', line)
                        if numbers:
                            stats['total_kmers'] = int(numbers[0])
                
                self.logger.info(f"📊 Jellyfish K-mer统计 | Jellyfish K-mer statistics: {stats['total_kmers']} distinct k-mers, {stats['unique_kmers']} unique k-mers")
                return stats
            else:
                self.logger.warning(f"⚠️ 无法获取Jellyfish统计 | Cannot get Jellyfish statistics: {result.stderr}")
                return {'total_kmers': 0, 'unique_kmers': 0}
        except Exception as e:
            self.logger.warning(f"⚠️ 获取Jellyfish统计信息时出错 | Error getting Jellyfish statistics: {e}")
            return {'total_kmers': 0, 'unique_kmers': 0}
