"""
K-mer分析结果处理模块 | K-mer Analysis Results Processing Module
"""

import pandas as pd
from pathlib import Path
from typing import Dict

class ResultsWriter:
    """结果写入器 | Results Writer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def save_matrices(self, count_matrix: pd.DataFrame, pa_matrix: pd.DataFrame):
        """保存矩阵文件 | Save matrix files"""
        output_dir = Path(self.config.output_dir)
        prefix = self.config.output_prefix
        
        # 保存计数矩阵 | Save count matrix
        count_file = output_dir / f"{prefix}_counts.csv"
        count_matrix.to_csv(count_file)
        self.logger.info(f"计数矩阵已保存 | Count matrix saved: {count_file}")
        
        # 保存PA矩阵 | Save PA matrix
        pa_file = output_dir / f"{prefix}_presence_absence.csv"
        pa_matrix.to_csv(pa_file)
        self.logger.info(f"存在/缺失矩阵已保存 | Presence/absence matrix saved: {pa_file}")
        
        # 保存转置矩阵 (样本为行) | Save transposed matrices (samples as rows)
        count_t_file = output_dir / f"{prefix}_counts_transposed.csv"
        count_matrix.T.to_csv(count_t_file)
        
        pa_t_file = output_dir / f"{prefix}_presence_absence_transposed.csv"
        pa_matrix.T.to_csv(pa_t_file)
        
        self.logger.info(f"转置矩阵已保存 | Transposed matrices saved")
        
        return {
            'count_matrix': str(count_file),
            'pa_matrix': str(pa_file),
            'count_matrix_t': str(count_t_file),
            'pa_matrix_t': str(pa_t_file)
        }
    
    def save_statistics(self, stats: Dict):
        """保存统计信息 | Save statistics"""
        output_dir = Path(self.config.output_dir)
        prefix = self.config.output_prefix
        
        # 保存基本统计 | Save basic statistics
        stats_file = output_dir / f"{prefix}_statistics.txt"
        with open(stats_file, 'w') as f:
            f.write("K-mer PAV分析统计报告 (双阶段设计) | K-mer PAV Analysis Statistics Report (Two-Stage Design)\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"总k-mer数 | Total k-mers: {stats['total_kmers']:,}\n")
            f.write(f"总样本数 | Total samples: {stats['total_samples']}\n")
            f.write(f"核心k-mer数 | Core k-mers: {stats['core_kmers']:,}\n")
            f.write(f"可变k-mer数 | Variable k-mers: {stats['variable_kmers']:,}\n")
            f.write(f"单例k-mer数 | Singleton k-mers: {stats['singleton_kmers']:,}\n\n")
            
            f.write("每个样本的k-mer数量 | K-mers per sample:\n")
            for sample, count in stats['kmers_per_sample'].items():
                f.write(f"  {sample}: {count:,}\n")
            
            f.write("\n每个样本的总计数 | Total counts per sample:\n")
            for sample, total in stats['total_counts_per_sample'].items():
                f.write(f"  {sample}: {total:,}\n")
        
        self.logger.info(f"统计信息已保存 | Statistics saved: {stats_file}")
        
        # 保存相似性矩阵 | Save similarity matrix
        similarity_file = output_dir / f"{prefix}_sample_similarity.csv"
        stats['sample_similarity'].to_csv(similarity_file)
        self.logger.info(f"样本相似性矩阵已保存 | Sample similarity matrix saved: {similarity_file}")
        
        return str(stats_file)

class SummaryGenerator:
    """总结生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self, output_files: Dict, stats: Dict):
        """生成总结报告 | Generate summary report"""
        output_dir = Path(self.config.output_dir)
        summary_file = output_dir / "analysis_summary.txt"
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("K-mer PAV分析总结报告 (双阶段设计) | K-mer PAV Analysis Summary Report (Two-Stage Design)\n")
            f.write("=" * 90 + "\n\n")
            
            # 分析设计说明 | Analysis design description
            f.write("分析设计 | Analysis Design:\n")
            f.write("  阶段1 | Phase 1: 从数据库文件构建统一k-mer数据库 | Build unified k-mer database from database files\n")
            f.write("  阶段2 | Phase 2: 查询文件与k-mer数据库比较分析 | Query files comparison against k-mer database\n")
            f.write("  样本处理 | Sample processing:\n")
            f.write("    - FASTQ文件: 整个文件作为一个样本 | FASTQ files: entire file as one sample\n")
            f.write("    - FASTA文件: 每条序列作为一个样本 | FASTA files: each sequence as one sample\n\n")
            
            # 输入参数 | Input parameters
            f.write("分析参数 | Analysis Parameters:\n")
            f.write(f"  数据库输入 | Database input: {self.config.database_input}\n")
            f.write(f"  查询输入 | Query input: {self.config.query_input}\n")
            f.write(f"  K-mer大小 | K-mer size: {self.config.kmer_size}\n")
            f.write(f"  最小计数 | Minimum count: {self.config.min_count}\n")
            f.write(f"  最大计数 | Maximum count: {self.config.max_count}\n")
            f.write(f"  反向互补 | Reverse complement: {'是 | Yes' if self.config.reverse_complement else '否 | No'}\n")
            f.write(f"  线程数 | Threads: {self.config.threads}\n")
            f.write(f"  KMC内存限制 | KMC memory limit: {self.config.kmc_memory_gb}GB\n")
            if self.config.database_pattern:
                f.write(f"  数据库文件模式 | Database file pattern: {self.config.database_pattern}\n")
            if self.config.query_pattern:
                f.write(f"  查询文件模式 | Query file pattern: {self.config.query_pattern}\n")
            f.write("\n")
            
            # 分析结果 | Analysis results
            f.write("分析结果 | Analysis Results:\n")
            f.write(f"  总k-mer数 | Total k-mers: {stats['total_kmers']:,}\n")
            f.write(f"  总样本数 | Total samples: {stats['total_samples']}\n")
            f.write(f"  核心k-mer数 | Core k-mers: {stats['core_kmers']:,} ({stats['core_kmers']/stats['total_kmers']*100:.1f}%)\n")
            f.write(f"  可变k-mer数 | Variable k-mers: {stats['variable_kmers']:,} ({stats['variable_kmers']/stats['total_kmers']*100:.1f}%)\n")
            f.write(f"  单例k-mer数 | Singleton k-mers: {stats['singleton_kmers']:,} ({stats['singleton_kmers']/stats['total_kmers']*100:.1f}%)\n")
            f.write("\n")
            
            # 输出文件 | Output files
            f.write("输出文件 | Output Files:\n")
            for key, filepath in output_files.items():
                f.write(f"  {key}: {filepath}\n")
            
            f.write("\n")
            f.write("文件说明 | File Descriptions:\n")
            f.write("  *_counts.csv: K-mer计数矩阵 (k-mer为行，样本为列) | K-mer count matrix (k-mers as rows, samples as columns)\n")
            f.write("  *_presence_absence.csv: 存在/缺失矩阵 (0/1二进制) | Presence/absence matrix (0/1 binary)\n")
            f.write("  *_transposed.csv: 转置矩阵 (样本为行，k-mer为列) | Transposed matrices (samples as rows, k-mers as columns)\n")
            f.write("  *_sample_similarity.csv: 样本间Jaccard相似性矩阵 | Sample Jaccard similarity matrix\n")
            f.write("  *_statistics.txt: 详细统计信息 | Detailed statistics\n")
            
        self.logger.info(f"总结报告已生成 | Summary report generated: {summary_file}")
