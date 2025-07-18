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
    
    def build_matrices(self, count_files: List[str], sample_names: List[str], 
                      kmer_sources: Dict[str, str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """构建k-mer矩阵 | Build k-mer matrices"""
        self.logger.info("=" * 60)
        self.logger.info("阶段3: 构建k-mer矩阵 | Phase 3: Building k-mer matrices")
        self.logger.info("=" * 60)
        
        # 1. 收集所有k-mer | Collect all k-mers
        all_kmers = self._collect_all_kmers(count_files)
        
        # 2. 读取样本数据 | Read sample data
        sample_data = self._read_sample_data(count_files, sample_names)
        
        # 3. 构建DataFrame | Build DataFrames
        count_matrix, pa_matrix = self._build_dataframes(sample_data, all_kmers, kmer_sources)
        
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
                         all_kmers: Set[str], kmer_sources: Dict[str, str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
        
        # 创建基础DataFrame | Create basic DataFrames
        count_matrix = pd.DataFrame(
            count_data,
            index=kmer_list,
            columns=sample_names
        )
        count_matrix.index.name = 'kmer'
        
        # 创建PA矩阵 (0/1) | Create PA matrix (0/1)
        pa_matrix = (count_matrix > 0).astype(int)
        
        # 分析k-mer分布并添加特征列 | Analyze k-mer distribution and add feature column
        kmer_features = self._analyze_kmer_features(pa_matrix, sample_names)
        
        # 添加from列（k-mer来源） | Add from column (k-mer source)
        kmer_from_sources = self._get_kmer_from_sources(kmer_list, kmer_sources)
        
        # 重置索引以便操作 | Reset index for manipulation
        count_matrix_with_features = count_matrix.reset_index()
        pa_matrix_with_features = pa_matrix.reset_index()
        
        # 在第二列位置插入from列，在第三列位置插入feature列 | Insert from column at second position, feature column at third position
        count_matrix_with_features.insert(1, 'from', kmer_from_sources)
        count_matrix_with_features.insert(2, 'feature', kmer_features)
        
        pa_matrix_with_features.insert(1, 'from', kmer_from_sources)
        pa_matrix_with_features.insert(2, 'feature', kmer_features)
        
        # 重新设置索引 | Reset index
        count_matrix_with_features.set_index('kmer', inplace=True)
        pa_matrix_with_features.set_index('kmer', inplace=True)
        
        self.logger.info(f"数据矩阵构建完成 | Data matrices built: {len(kmer_list)} k-mers × {len(sample_names)} samples")
        
        # 统计from列 | Count from column
        from_counts = pd.Series(kmer_from_sources).value_counts()
        self.logger.info(f"k-mer来源统计 | K-mer source statistics:")
        for source, count in from_counts.items():
            self.logger.info(f"  - {source}: {count:,} k-mers")
        
        # 统计feature列 | Count feature column
        feature_counts = pd.Series(kmer_features).value_counts()
        self.logger.info(f"k-mer特征统计 | K-mer feature statistics:")
        for feature, count in feature_counts.items():
            self.logger.info(f"  - {feature}: {count:,} k-mers")
        
        return count_matrix_with_features, pa_matrix_with_features
    
    def _analyze_kmer_features(self, pa_matrix: pd.DataFrame, sample_names: List[str]) -> List[str]:
        """分析k-mer特征 | Analyze k-mer features"""
        self.logger.info("分析k-mer分布特征 | Analyzing k-mer distribution features")
        
        features = []
        total_samples = len(sample_names)
        
        for kmer in pa_matrix.index:
            # 统计该k-mer在多少个样本中出现 | Count how many samples contain this k-mer
            present_samples = pa_matrix.loc[kmer]
            present_count = present_samples.sum()
            
            if present_count == 1:
                # 单例k-mer，找到唯一包含它的样本 | Singleton k-mer, find the unique sample containing it
                unique_sample = present_samples[present_samples == 1].index[0]
                features.append(unique_sample)
            elif present_count == total_samples:
                # 核心k-mer，所有样本都包含 | Core k-mer, present in all samples
                features.append("common")
            else:
                # 可变k-mer，在部分样本中出现 | Variable k-mer, present in some samples
                # 如果样本数不多，可以列出所有包含该k-mer的样本 | If not too many samples, list all containing samples
                if present_count <= 3:  # 如果包含该k-mer的样本数 <= 3，列出所有样本名
                    present_sample_names = present_samples[present_samples == 1].index.tolist()
                    features.append("|".join(sorted(present_sample_names)))
                else:
                    # 否则标记为variable | Otherwise mark as variable
                    features.append(f"variable({present_count}/{total_samples})")
        
        return features
    
    def _get_kmer_from_sources(self, kmer_list: List[str], kmer_sources: Dict[str, str]) -> List[str]:
        """获取k-mer来源信息 | Get k-mer source information"""
        self.logger.info("添加k-mer来源信息 | Adding k-mer source information")
        
        from_sources = []
        for kmer in kmer_list:
            # 从kmer_sources字典获取来源信息，如果没有则标记为unknown
            source = kmer_sources.get(kmer, "unknown")
            from_sources.append(source)
        
        return from_sources
    
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
        
        # 获取实际的样本列（排除kmer索引、from列和feature列） | Get actual sample columns (excluding kmer index, from and feature columns)
        sample_columns = [col for col in count_matrix.columns if col not in ['from', 'feature']]
        
        # 基本统计 | Basic statistics
        stats['total_kmers'] = len(count_matrix)
        stats['total_samples'] = len(sample_columns)
        
        # K-mer统计（只使用样本列） | K-mer statistics (only use sample columns)
        pa_sample_data = pa_matrix[sample_columns]
        stats['kmers_per_sample'] = pa_sample_data.sum(axis=0).to_dict()
        stats['samples_per_kmer'] = pa_sample_data.sum(axis=1).to_dict()
        
        # 核心/可变k-mer（基于feature列） | Core/variable k-mers (based on feature column)
        feature_counts = count_matrix['feature'].value_counts()
        stats['core_kmers'] = feature_counts.get('common', 0)
        
        # 计算单例k-mer（feature列中不是common也不是variable的） | Calculate singleton k-mers
        singleton_count = 0
        variable_count = 0
        
        for feature, count in feature_counts.items():
            if feature == 'common':
                continue
            elif feature.startswith('variable('):
                variable_count += count
            else:
                # 单例或少数样本特有的k-mer | Singleton or k-mers specific to few samples
                singleton_count += count
        
        stats['singleton_kmers'] = singleton_count
        stats['variable_kmers'] = variable_count
        
        # 特征分布统计 | Feature distribution statistics
        stats['feature_distribution'] = feature_counts.to_dict()
        
        # from列分布统计 | From column distribution statistics
        from_counts = count_matrix['from'].value_counts()
        stats['from_distribution'] = from_counts.to_dict()
        
        # 平均计数统计（只使用样本列） | Average count statistics (only use sample columns)
        count_sample_data = count_matrix[sample_columns]
        stats['avg_count_per_kmer'] = count_sample_data.mean(axis=1).to_dict()
        stats['total_counts_per_sample'] = count_sample_data.sum(axis=0).to_dict()
        
        # 样本相似性（只使用样本列） | Sample similarity (only use sample columns)
        sample_jaccard = self._calculate_jaccard_matrix(pa_sample_data)
        stats['sample_similarity'] = sample_jaccard
        
        self.logger.info(f"统计完成 | Statistics completed:")
        self.logger.info(f"  - 总k-mer数 | Total k-mers: {stats['total_kmers']:,}")
        self.logger.info(f"  - 核心k-mer数 | Core k-mers: {stats['core_kmers']:,}")
        self.logger.info(f"  - 可变k-mer数 | Variable k-mers: {stats['variable_kmers']:,}")
        self.logger.info(f"  - 特有k-mer数 | Specific k-mers: {stats['singleton_kmers']:,}")
        
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
