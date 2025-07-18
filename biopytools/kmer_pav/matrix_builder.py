"""
K-mer矩阵构建模块 | K-mer Matrix Builder Module
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Set, List, Tuple
from .utils import cleanup_files

class KmerMatrixBuilder:
    """K-mer矩阵构建器 | K-mer Matrix Builder"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def build_matrices(self, count_files: List[str], sample_names: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """构建k-mer矩阵 | Build k-mer matrices"""
        self.logger.info("=" * 60)
        self.logger.info("阶段3: 构建k-mer矩阵 | Phase 3: Building k-mer matrices")
        self.logger.info("=" * 60)
        
        # 1. 收集所有k-mer | Collect all k-mers
        all_kmers = self._collect_all_kmers(count_files)
        
        # 2. 读取样本数据 | Read sample data
        sample_data = self._read_sample_data(count_files, sample_names)
        
        # 3. 构建DataFrame | Build DataFrames
        count_matrix, pa_matrix = self._build_dataframes(sample_data, all_kmers)
        
        # 4. 清理计数文件 | Clean up count files
        if not self.config.keep_intermediate:
            self._cleanup_count_files(count_files)
        
        return count_matrix, pa_matrix
    
    def _collect_all_kmers(self, count_files: List[str]) -> Set[str]:
        """收集所有k-mer | Collect all k-mers"""
        self.logger.info("收集所有k-mer | Collecting all k-mers")
        all_kmers = set()
        
        for count_file in count_files:
            if count_file is None:
                continue
                
            file_path = self.config.output_path / count_file
            if file_path.exists():
                try:
                    with open(file_path, 'r') as f:
                        for line in f:
                            if line.strip():
                                kmer = line.split('\t')[0]
                                all_kmers.add(kmer)
                except Exception as e:
                    self.logger.error(f"读取文件失败 | Failed to read file {count_file}: {e}")
        
        self.logger.info(f"收集到 {len(all_kmers)} 个唯一k-mer | Collected {len(all_kmers)} unique k-mers")
        return all_kmers
    
    def _read_sample_data(self, count_files: List[str], sample_names: List[str]) -> Dict[str, Dict[str, int]]:
        """读取样本数据 | Read sample data"""
        self.logger.info("读取样本数据 | Reading sample data")
        sample_data = {}
        
        for i, count_file in enumerate(count_files):
            if count_file is None:
                sample_data[sample_names[i]] = {}
                continue
                
            sample_name = sample_names[i]
            file_path = self.config.output_path / count_file
            sample_kmers = {}
            
            if file_path.exists():
                try:
                    with open(file_path, 'r') as f:
                        for line in f:
                            if line.strip():
                                parts = line.strip().split('\t')
                                if len(parts) >= 2:
                                    kmer, count = parts[0], int(parts[1])
                                    sample_kmers[kmer] = count
                    
                    sample_data[sample_name] = sample_kmers
                    self.logger.info(f"样本 {sample_name}: {len(sample_kmers)} 个k-mer | Sample {sample_name}: {len(sample_kmers)} k-mers")
                    
                except Exception as e:
                    self.logger.error(f"处理样本数据失败 | Failed to process sample data {sample_name}: {e}")
                    sample_data[sample_name] = {}
            else:
                self.logger.warning(f"样本计数文件不存在 | Sample count file not found: {count_file}")
                sample_data[sample_name] = {}
        
        return sample_data
    
    def _build_dataframes(self, sample_data: Dict[str, Dict[str, int]], 
                         all_kmers: Set[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """构建pandas DataFrame | Build pandas DataFrames"""
        self.logger.info("构建数据矩阵 | Building data matrices")
        
        # 排序以确保一致性 | Sort for consistency
        sample_names = sorted(sample_data.keys())
        kmer_list = sorted(all_kmers)
        
        # 初始化矩阵 | Initialize matrices
        count_data = np.zeros((len(kmer_list), len(sample_names)), dtype=int)
        
        # 填充数据 | Fill data
        for j, sample_name in enumerate(sample_names):
            sample_kmers = sample_data[sample_name]
            for i, kmer in enumerate(kmer_list):
                count_data[i, j] = sample_kmers.get(kmer, 0)
        
        # 创建DataFrame | Create DataFrames
        count_matrix = pd.DataFrame(
            count_data,
            index=kmer_list,
            columns=sample_names
        )
        count_matrix.index.name = 'kmer'
        
        # 创建PA矩阵 (0/1) | Create PA matrix (0/1)
        pa_matrix = (count_matrix > 0).astype(int)
        
        self.logger.info(f"数据矩阵构建完成 | Data matrices built: {len(kmer_list)} k-mers × {len(sample_names)} samples")
        
        return count_matrix, pa_matrix
    
    def _cleanup_count_files(self, count_files: List[str]):
        """清理计数文件 | Clean up count files"""
        self.logger.info("清理计数文件 | Cleaning up count files")
        for count_file in count_files:
            if count_file:
                file_path = self.config.output_path / count_file
                if file_path.exists():
                    try:
                        file_path.unlink()
                        self.logger.debug(f"删除计数文件 | Removed count file: {count_file}")
                    except Exception as e:
                        self.logger.warning(f"删除计数文件失败 | Failed to remove count file {count_file}: {e}")

class StatisticsCalculator:
    """统计计算器 | Statistics Calculator"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def calculate_statistics(self, count_matrix: pd.DataFrame, 
                           pa_matrix: pd.DataFrame) -> Dict:
        """计算统计信息 | Calculate statistics"""
        self.logger.info("计算统计信息 | Calculating statistics")
        
        stats = {}
        
        # 基本统计 | Basic statistics
        stats['total_kmers'] = len(count_matrix)
        stats['total_samples'] = len(count_matrix.columns)
        
        # K-mer统计 | K-mer statistics
        stats['kmers_per_sample'] = pa_matrix.sum(axis=0).to_dict()
        stats['samples_per_kmer'] = pa_matrix.sum(axis=1).to_dict()
        
        # 核心/可变k-mer | Core/variable k-mers
        kmer_presence = pa_matrix.sum(axis=1)
        stats['core_kmers'] = len(kmer_presence[kmer_presence == stats['total_samples']])
        stats['variable_kmers'] = len(kmer_presence[kmer_presence < stats['total_samples']])
        stats['singleton_kmers'] = len(kmer_presence[kmer_presence == 1])
        
        # 平均计数统计 | Average count statistics
        stats['avg_count_per_kmer'] = count_matrix.mean(axis=1).to_dict()
        stats['total_counts_per_sample'] = count_matrix.sum(axis=0).to_dict()
        
        # 样本相似性 | Sample similarity
        sample_jaccard = self._calculate_jaccard_matrix(pa_matrix)
        stats['sample_similarity'] = sample_jaccard
        
        self.logger.info(f"统计完成 | Statistics completed:")
        self.logger.info(f"  - 总k-mer数 | Total k-mers: {stats['total_kmers']:,}")
        self.logger.info(f"  - 核心k-mer数 | Core k-mers: {stats['core_kmers']:,}")
        self.logger.info(f"  - 可变k-mer数 | Variable k-mers: {stats['variable_kmers']:,}")
        self.logger.info(f"  - 单例k-mer数 | Singleton k-mers: {stats['singleton_kmers']:,}")
        
        return stats
    
    def _calculate_jaccard_matrix(self, pa_matrix: pd.DataFrame) -> pd.DataFrame:
        """计算Jaccard相似性矩阵 | Calculate Jaccard similarity matrix"""
        samples = pa_matrix.columns
        jaccard_matrix = pd.DataFrame(index=samples, columns=samples, dtype=float)
        
        for i, sample1 in enumerate(samples):
            for j, sample2 in enumerate(samples):
                if i <= j:
                    set1 = set(pa_matrix.index[pa_matrix[sample1] == 1])
                    set2 = set(pa_matrix.index[pa_matrix[sample2] == 1])
                    
                    intersection = len(set1 & set2)
                    union = len(set1 | set2)
                    
                    jaccard = intersection / union if union > 0 else 0
                    jaccard_matrix.loc[sample1, sample2] = jaccard
                    jaccard_matrix.loc[sample2, sample1] = jaccard
        
        return jaccard_matrix
