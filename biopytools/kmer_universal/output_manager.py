"""
输出管理模块 | Output Management Module 🖨️📄
"""

import os
import csv
import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
import pandas as pd
import numpy as np

from .position_tracker import KmerInfo, KmerPosition

class OutputManager:
    """输出管理器 🖨️"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # 确保输出目录存在 📂
        Path(self.config.output_dir).mkdir(parents=True, exist_ok=True)
    
    def write_kmer_library(self, kmers: Dict[str, KmerInfo], output_file: str):
        """输出k-mer库文件（FASTA格式） 🧬📄"""
        self.logger.info(f"✍️ Writing k-mer library to {output_file}")
        
        with open(output_file, 'w') as f:
            for kmer_seq, kmer_info in kmers.items():
                # 生成FASTA标题 🏷️
                if kmer_info.source_type == 'fasta' and kmer_info.positions:
                    # FASTA来源：使用位置信息
                    pos = kmer_info.positions[0]  # 使用第一个位置
                    header = f">{pos.seq_name}_{pos.start_pos}_{pos.end_pos}"
                else:
                    # FASTQ来源：使用样品名和序号
                    header = f">{kmer_info.sample_name}_kmer_{abs(hash(kmer_seq)) % 1000000:06d}"
                
                f.write(f"{header}\n{kmer_seq}\n")
        
        self.logger.info(f"✅ Written {len(kmers)} k-mers to library file")
    
    def write_abundance_matrix(self, kmer_library: Dict[str, KmerInfo], 
                             target_abundances: Dict[str, Dict[str, int]], 
                             output_file: str):
        """输出丰度矩阵 📊"""
        self.logger.info(f"✍️ Writing abundance matrix to {output_file}")
        
        # 获取所有目标样品名
        target_samples = list(target_abundances.keys())
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # 写入表头 🏷️
            if any(info.source_type == 'fasta' for info in kmer_library.values()):
                # 包含FASTA来源的k-mer，需要位置信息
                header = ['seq_name', 'start_pos', 'end_pos', 'kmer_seq'] + target_samples
            else:
                # 仅FASTQ来源
                header = ['sample_name', 'kmer_id', 'kmer_seq'] + target_samples
            
            writer.writerow(header)
            
            # 写入数据行 📈
            for kmer_seq, kmer_info in kmer_library.items():
                if kmer_info.source_type == 'fasta' and kmer_info.positions:
                    # FASTA来源
                    pos = kmer_info.positions[0]
                    row = [pos.seq_name, pos.start_pos, pos.end_pos, kmer_seq]
                else:
                    # FASTQ来源
                    pos = kmer_info.positions[0] if kmer_info.positions else None
                    sample_name = kmer_info.sample_name
                    kmer_id = pos.seq_name if pos else f"kmer_{abs(hash(kmer_seq)) % 1000000:06d}"
                    row = [sample_name, kmer_id, kmer_seq]
                
                # 添加各样品的丰度
                for sample in target_samples:
                    abundance = target_abundances[sample].get(kmer_seq, 0)
                    row.append(abundance)
                
                writer.writerow(row)
        
        self.logger.info("✅ Abundance matrix written successfully")
    
    def write_presence_matrix(self, abundance_file: str, output_file: str):
        """基于丰度矩阵生成存在/缺失矩阵 ✅❌"""
        self.logger.info(f"⚙️ Generating presence matrix from {abundance_file}")
        
        # 读取丰度矩阵 🐼
        df = pd.read_csv(abundance_file)
        
        # 确定样品列（数值列） 🎯
        if 'seq_name' in df.columns:
            sample_cols = df.columns[4:]  # FASTA格式
        else:
            sample_cols = df.columns[3:]  # FASTQ格式
        
        # 转换为0/1矩阵 🔄
        presence_df = df.copy()
        for col in sample_cols:
            presence_df[col] = (df[col] > 0).astype(int)
        
        # 写入文件 ✍️
        presence_df.to_csv(output_file, index=False)
        self.logger.info("✅ Presence matrix written successfully")
    
    def write_sliding_window_analysis(self, kmer_library: Dict[str, KmerInfo], 
                                    target_abundances: Dict[str, Dict[str, int]], 
                                    output_file: str):
        """滑窗分析输出 🪟"""
        self.logger.info(f"🏃‍♀️ Performing sliding window analysis, output to {output_file}")
        
        # 只对FASTA来源的k-mer进行滑窗分析 🧬
        fasta_kmers = {k: v for k, v in kmer_library.items() if v.source_type == 'fasta'}
        
        if not fasta_kmers:
            self.logger.warning("⚠️ No FASTA-derived k-mers found for sliding window analysis")
            return
        
        # 按序列和位置组织k-mer 🗺️
        seq_kmers = {}
        for kmer_seq, kmer_info in fasta_kmers.items():
            for pos in kmer_info.positions:
                seq_name = pos.seq_name
                if seq_name not in seq_kmers:
                    seq_kmers[seq_name] = []
                seq_kmers[seq_name].append((pos.start_pos, kmer_seq))
        
        # 为每个序列执行滑窗分析 ⚙️
        results = []
        for seq_name, kmer_positions in seq_kmers.items():
            # 按位置排序
            kmer_positions.sort()
            
            if not kmer_positions:
                continue
            
            # 确定序列长度
            max_pos = max(pos + self.config.kmer_size - 1 for pos, _ in kmer_positions)
            
            # 滑窗分析
            for window_size in self.config.window_sizes:
                window_results = self._sliding_window_for_sequence(
                    seq_name, kmer_positions, target_abundances, 
                    max_pos, window_size
                )
                results.extend(window_results)
        
        # 写入结果 ✍️
        self._write_sliding_window_results(results, output_file)
    
    def _sliding_window_for_sequence(self, seq_name: str, kmer_positions: List[Tuple[int, str]], 
                                   target_abundances: Dict[str, Dict[str, int]], 
                                   seq_length: int, window_size: int) -> List[Dict]:
        """对单个序列执行滑窗分析 🧬⚙️"""
        results = []
        
        # 滑窗步长（可配置，默认为窗口大小的一半） 📏
        step_size = window_size // 2
        
        for start in range(1, seq_length, step_size):
            end = min(start + window_size - 1, seq_length)
            
            # 找到窗口内的k-mer 🔍
            window_kmers = [
                kmer for pos, kmer in kmer_positions 
                if start <= pos <= end - self.config.kmer_size + 1
            ]
            
            if not window_kmers:
                continue
            
            # 计算每个样品的统计 🧮
            for sample_name, abundances in target_abundances.items():
                present_kmers = sum(1 for kmer in window_kmers if abundances.get(kmer, 0) > 0)
                total_kmers = len(window_kmers)
                ratio = present_kmers / total_kmers if total_kmers > 0 else 0.0
                
                results.append({
                    'seq_name': seq_name,
                    'window_start': start,
                    'window_end': end,
                    'window_size': window_size,
                    'sample_name': sample_name,
                    'total_kmers': total_kmers,
                    'present_kmers': present_kmers,
                    'presence_ratio': ratio
                })
        
        return results
    
    def _write_sliding_window_results(self, results: List[Dict], output_file: str):
        """写入滑窗分析结果 ✍️📄"""
        if not results:
            self.logger.warning("⚠️ No sliding window results to write")
            return
        
        df = pd.DataFrame(results)
        df.to_csv(output_file, index=False)
        self.logger.info(f"✅ Sliding window analysis results written to {output_file}")
    
    def write_summary_report(self, kmer_library: Dict[str, KmerInfo], 
                           target_abundances: Dict[str, Dict[str, int]], 
                           output_file: str):
        """生成分析摘要报告 📝"""
        self.logger.info(f"✍️ Generating summary report: {output_file}")
        
        with open(output_file, 'w') as f:
            f.write("=== K-mer Analysis Summary Report 📋 ===\n\n")
            
            # 基本信息 ℹ️
            f.write("Basic Information:\n")
            f.write(f"  K-mer size: {self.config.kmer_size}\n")
            f.write(f"  Total unique k-mers: {len(kmer_library)}\n")
            f.write(f"  Target samples: {len(target_abundances)}\n\n")
            
            # K-mer来源统计 📚
            fasta_kmers = sum(1 for info in kmer_library.values() if info.source_type == 'fasta')
            fastq_kmers = sum(1 for info in kmer_library.values() if info.source_type == 'fastq')
            
            f.write("K-mer Sources:\n")
            f.write(f"  From FASTA files: {fasta_kmers}\n")
            f.write(f"  From FASTQ files: {fastq_kmers}\n\n")
            
            # 样品统计 🎯
            f.write("Sample Statistics:\n")
            for sample_name, abundances in target_abundances.items():
                total_abundance = sum(abundances.values())
                present_kmers = sum(1 for count in abundances.values() if count > 0)
                
                f.write(f"  {sample_name}:\n")
                f.write(f"    Present k-mers: {present_kmers}/{len(kmer_library)} ({present_kmers/len(kmer_library)*100:.1f}%)\n")
                f.write(f"    Total abundance: {total_abundance}\n")
            
            f.write("\n")
            
            # 配置信息 ⚙️
            f.write("Configuration:\n")
            f.write(f"  Threads: {self.config.threads}\n")
            f.write(f"  Memory: {self.config.memory_gb}GB\n")
            f.write(f"  Window sizes: {self.config.window_sizes}\n")
            f.write(f"  Output directory: {self.config.output_dir}\n")
        
        self.logger.info("✅ Summary report generated successfully")