"""
KMC k-mer查询模块|KMC K-mer Query Module
"""

import os
import h5py
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Union
from concurrent.futures import ThreadPoolExecutor, as_completed

from .config import KMCConfig
from .utils import KMCLogger, format_number


# 模块级多线程查询工作函数|Module-level multithreading query worker function
def _query_kmer_id_worker(args):
    """多线程查询kmer_id工作函数|Multithreading kmer_id query worker function

    Args:
        args: (items, db) 元组，items是(kmer_name, kmer_seq)列表，db是LevelDB对象

    Returns:
        {kmer_name: kmer_id} 字典
    """
    items, db = args

    result = {}
    for kmer_name, kmer_seq in items:
        if db is not None:
            kid_bytes = db.get(kmer_seq.encode())
            if kid_bytes is not None:
                result[kmer_name] = int(kid_bytes)
        # 内存字典模式在主线程处理（已经很快）

    return result


class KMCQuery:
    """KMC k-mer查询器|KMC K-mer Query"""

    def __init__(self, config: KMCConfig, logger: Optional[KMCLogger] = None):
        """初始化查询器|Initialize query"""
        self.config = config
        self.config.validate()

        if logger is None:
            self.logger_manager = KMCLogger(self.config.output_path)
            self.logger = self.logger_manager.get_logger()
        else:
            self.logger_manager = logger
            self.logger = logger.get_logger()

        self.matrix_file = self.config.output_path / 'abundance_matrix.h5'
        self.kmer_dict_file = self.config.output_path / 'kmer_dictionary.h5'

        # 加载索引|Load indices
        self._kmer_to_id = None
        self._sample_names = None
        self._db = None  # 数据库索引（LevelDB 或 RocksDB）|Database index (LevelDB or RocksDB)
        self._db_type = None  # 'leveldb' 或 'rocksdb'

    def _load_indices(self):
        """加载索引|Load indices（优先使用数据库索引，回退到内存字典）|Load indices (prefer database index, fallback to memory dict)"""
        # 尝试加载数据库索引（LevelDB 或 RocksDB）|Try to load database index (LevelDB or RocksDB)
        if self._db is None:
            self._db, self._db_type = self._load_db_index()

            if self._db is not None:
                db_name = 'RocksDB' if self._db_type == 'rocksdb' else 'LevelDB'
                self.logger.info(f"使用 {db_name} 索引（低内存占用）|Using {db_name} index (low memory)")
            else:
                # 数据库不可用，使用内存字典|Database not available, use memory dict
                self.logger.info("数据库索引不可用，使用内存字典（高内存占用）|Database index not available, using memory dict (high memory)")
                if self._kmer_to_id is None:
                    self._kmer_to_id = {}
                    with h5py.File(self.kmer_dict_file, 'r') as f:
                        kmer_array = f['kmer'][:]
                        for i, kmer in enumerate(kmer_array):
                            self._kmer_to_id[kmer.decode('utf-8')] = i
                    self.logger.warning(f"已加载 {len(self._kmer_to_id):,} 个 k-mer 到内存|Loaded {len(self._kmer_to_id):,} k-mers to memory")

        # 加载样本名称|Load sample names
        if self._sample_names is None:
            with h5py.File(self.matrix_file, 'r') as f:
                self._sample_names = [s.decode('utf-8') for s in f['sample_names'][:]]

    def _load_db_index(self):
        """加载数据库索引（LevelDB 或 RocksDB）|Load database index (LevelDB or RocksDB)

        Returns:
            (数据库对象, 数据库类型) 或 (None, None)|(DB object, DB type) or (None, None)
        """
        import json

        # 检查元数据文件|Check metadata file
        metadata_file = self.config.output_path / 'kmer_metadata.json'
        if not metadata_file.exists():
            self.logger.info("未找到 kmer_metadata.json，数据库索引不可用|kmer_metadata.json not found, database index unavailable")
            return None, None

        try:
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)

            # 检查索引类型|Check index type
            index_type = metadata.get('index_type')
            if index_type not in ['leveldb', 'rocksdb']:
                self.logger.info(f"不支持的索引类型: {index_type}|Unsupported index type: {index_type}")
                return None, None

            # 检查索引路径|Check index path
            index_path = metadata.get('index_path')
            if not index_path or not Path(index_path).exists():
                self.logger.info(f"数据库索引不存在: {index_path}|Database index not exists: {index_path}")
                return None, None

            # 根据类型加载数据库|Load database based on type
            if index_type == 'rocksdb':
                # 加载 RocksDB|Load RocksDB
                from .matrix_builder import DB
                db = DB(str(index_path), create_if_missing=False,
                       max_memory_gb=self.config.max_memory)
                self.logger.info(f"RocksDB 索引加载成功|RocksDB index loaded: {index_path}")
                return db, 'rocksdb'
            else:
                # 加载 LevelDB|Load LevelDB
                import plyvel
                db = plyvel.DB(index_path, create_if_missing=False, error_if_exists=False)
                self.logger.info(f"LevelDB 索引加载成功|LevelDB index loaded: {index_path}")
                return db, 'leveldb'

        except ImportError as e:
            if index_type == 'rocksdb':
                self.logger.warning(f"RocksDB 加载失败（可能python-rocksdb未安装）|RocksDB load failed (python-rocksdb not installed): {e}")
            else:
                self.logger.warning(f"plyvel 未安装，无法使用 LevelDB|plyvel not installed, cannot use LevelDB: {e}")
            return None, None
        except Exception as e:
            self.logger.warning(f"加载数据库索引失败|Failed to load database index: {e}")
            return None, None

    def query_kmer(self, kmer_seq: str) -> Dict[str, int]:
        """查询单个k-mer的丰度|Query single k-mer abundance

        Args:
            kmer_seq: k-mer序列|K-mer sequence

        Returns:
            样本名→丰度的映射|Sample name to abundance mapping
        """
        self._load_indices()

        self.logger.info(f"查询k-mer|Querying k-mer: {kmer_seq}")

        # 根据索引类型查找 kmer_id|Lookup kmer_id based on index type
        kmer_id = None
        if self._db is not None:
            # 使用数据库查找|Use database lookup
            kid_bytes = self._db.get(kmer_seq.encode())
            if kid_bytes is None:
                self.logger.warning(f"k-mer不存在|K-mer does not exist: {kmer_seq}")
                return {}
            kmer_id = int(kid_bytes)
        else:
            # 使用内存字典查找|Use memory dict lookup
            if kmer_seq not in self._kmer_to_id:
                self.logger.warning(f"k-mer不存在|K-mer does not exist: {kmer_seq}")
                return {}
            kmer_id = self._kmer_to_id[kmer_seq]

        abundances = {}

        with h5py.File(self.matrix_file, 'r') as f:
            storage_type = f.attrs['storage_type'].decode('utf-8')

            if storage_type == 'sparse':
                # 稀疏存储查询|Sparse storage query
                # 找到所有匹配的kmer_id|Find all matching kmer_id
                mask = f['kmer_id'][:] == kmer_id

                if mask.any():
                    sample_ids = f['sample_id'][:][mask]
                    abunds = f['abundance'][:][mask]

                    for sid, ab in zip(sample_ids, abunds):
                        if ab > 0:
                            sample_name = self._sample_names[sid]
                            abundances[sample_name] = int(ab)
            else:
                # 密集矩阵查询|Dense matrix query
                abunds = f['abundance'][kmer_id, :]

                for i, ab in enumerate(abunds):
                    if ab > 0:
                        sample_name = self._sample_names[i]
                        abundances[sample_name] = int(ab)

        self.logger.info(f"查询完成|Query completed: 在 {len(abundances)} 个样本中找到|found in {len(abundances)} samples")

        return abundances

    def query_kmers(self, kmer_seqs: List[str]) -> Dict[str, Dict[str, int]]:
        """查询多个k-mer的丰度|Query multiple k-mers abundances

        注意：此方法已弃用，请使用 query_kmers_from_fasta 进行批量查询
        Note: This method is deprecated, use query_kmers_from_fasta for batch queries
        """
        self.logger.info(f"查询 {len(kmer_seqs)} 个k-mer|Querying {len(kmer_seqs)} k-mers")

        results = {}

        for i, kmer_seq in enumerate(kmer_seqs, 1):
            self.logger.info(f"进度|Progress: {i}/{len(kmer_seqs)}")
            results[kmer_seq] = self.query_kmer(kmer_seq)

        return results

    def query_kmers_batch(self, kmer_dict: Dict[str, str]) -> Dict[str, Dict]:
        """批量查询k-mer丰度（优化版，支持百万级）|Batch query k-mer abundances (optimized, supports millions)

        Args:
            kmer_dict: kmer_name → kmer_sequence 的映射|kmer_name to kmer_sequence mapping
                      例如: {"seq1": "ACGT...", "seq2": "GCTA..."}

        Returns:
            kmer_name → {sequence: "...", sample_name: abundance, ...} 的映射
            kmer_name to {sequence: "...", sample_name: abundance, ...} mapping
        """
        self.logger.info(f"开始批量查询|Starting batch query: {len(kmer_dict)} k-mers")

        # 加载索引|Load indices
        self._load_indices()

        if not kmer_dict:
            return {}

        # 步骤1: 批量查找 kmer_id|Step 1: Batch lookup kmer_id
        self.logger.info("步骤1/3: 批量查找 kmer_id|Step 1/3: Batch lookup kmer_id")
        kmer_name_to_id = {}
        missing_kmers = []

        if self._db is not None:
            # 使用数据库多线程批量查询|Use database multithreading batch query
            # 转换为列表方便分块|Convert to list for chunking
            items_list = list(kmer_dict.items())

            # 根据 k-mer 数量决定是否使用多线程|Use multithreading based on k-mer count
            if len(items_list) > 10000:  # 超过1万条使用多线程
                n_threads = min(self.config.threads, 12)  # 最多12线程
                self.logger.info(f"使用 {n_threads} 线程并行查询|Using {n_threads} threads for parallel query")

                # 数据分块|Split data into chunks
                chunk_size = len(items_list) // n_threads + 1
                chunks = []
                for i in range(0, len(items_list), chunk_size):
                    chunks.append(items_list[i:i+chunk_size])

                self.logger.info(f"数据分块|Data chunks: {len(chunks)}, 每块约|~{chunk_size:,} 条")

                # 多线程并行查询|Multithreaded parallel query
                chunk_args = [(chunk, self._db) for chunk in chunks]

                with ThreadPoolExecutor(max_workers=n_threads) as executor:
                    future_to_chunk = {executor.submit(_query_kmer_id_worker, args): args[0] for args in chunk_args}

                    for future in as_completed(future_to_chunk):
                        chunk = future_to_chunk[future]
                        try:
                            result = future.result()
                            kmer_name_to_id.update(result)
                        except Exception as e:
                            self.logger.error(f"查询块出错|Error querying chunk: {e}")
                            raise
            else:
                # 少量数据直接查询|Direct query for small data
                for kmer_name, kmer_seq in kmer_dict.items():
                    kid_bytes = self._db.get(kmer_seq.encode())
                    if kid_bytes is not None:
                        kmer_name_to_id[kmer_name] = int(kid_bytes)
                    else:
                        missing_kmers.append(kmer_name)
        else:
            # 使用内存字典（已经很快，不需要多线程）|Use memory dict (fast enough, no need for multithreading)
            for kmer_name, kmer_seq in kmer_dict.items():
                if kmer_seq in self._kmer_to_id:
                    kmer_name_to_id[kmer_name] = self._kmer_to_id[kmer_seq]
                else:
                    missing_kmers.append(kmer_name)

        # 统计未找到的k-mer|Count missing k-mers
        if self._db is not None and len(kmer_dict) > 10000:
            # 多线程模式需要统计missing
            found_count = len(kmer_name_to_id)
            missing_count = len(kmer_dict) - found_count
            if missing_count > 0:
                # 找出哪些k-mer没找到（仅在少量时）|Find which k-mers are missing (only for small amounts)
                if missing_count <= 100:
                    missing_kmers = [name for name, seq in kmer_dict.items()
                                   if name not in kmer_name_to_id]
        else:
            # 单线程模式已经收集了missing_kmers|Single-threaded mode already collected missing_kmers
            pass

        found_count = len(kmer_name_to_id)
        self.logger.info(f"找到 kmer_id|Found kmer_id: {found_count}/{len(kmer_dict)} ({found_count/len(kmer_dict)*100:.1f}%)")

        if missing_kmers:
            self.logger.warning(f"未找到的 k-mer|Missing k-mers: {len(missing_kmers)}")
            if len(missing_kmers) <= 10:
                self.logger.warning(f"  {', '.join(missing_kmers)}")

        # 步骤2: 从 HDF5 批量查询丰度|Step 2: Batch query abundances from HDF5
        self.logger.info("步骤2/3: 从 HDF5 批量查询丰度|Step 2/3: Batch query abundances from HDF5")

        results = {}

        # 初始化所有 k-mer 的结果（包括未找到的），包含序列信息|Initialize results with sequence info
        for kmer_name, kmer_seq in kmer_dict.items():
            results[kmer_name] = {'sequence': kmer_seq}

        if not kmer_name_to_id:
            self.logger.warning("没有找到任何 k-mer，跳过 HDF5 查询|No k-mers found, skipping HDF5 query")
            return results

        with h5py.File(self.matrix_file, 'r') as f:
            storage_type = f.attrs['storage_type']
            # 兼容字符串和bytes|Compatible with both str and bytes
            if isinstance(storage_type, bytes):
                storage_type = storage_type.decode('utf-8')

            if storage_type == 'sparse':
                # 稀疏存储：批量筛选|Sparse storage: batch filter
                kmer_id_array = f['kmer_id'][:]
                sample_id_array = f['sample_id'][:]
                abundance_array = f['abundance'][:]

                # 构建目标 kmer_id 集合|Build target kmer_id set
                target_ids = set(kmer_name_to_id.values())
                self.logger.info(f"目标 kmer_id 数量|Target kmer_id count: {len(target_ids)}")

                # 批量筛选|Batch filter
                mask = np.isin(kmer_id_array, list(target_ids))

                self.logger.info(f"匹配的数据行数|Matching rows: {mask.sum()}")

                if mask.any():
                    matched_kmer_ids = kmer_id_array[mask]
                    matched_sample_ids = sample_id_array[mask]
                    matched_abundances = abundance_array[mask]

                    # 构建 kmer_id → kmer_name 反向映射|Build kmer_id → kmer_name reverse mapping
                    id_to_name = {v: k for k, v in kmer_name_to_id.items()}

                    # 填充结果|Fill results
                    for kid, sid, ab in zip(matched_kmer_ids, matched_sample_ids, matched_abundances):
                        if ab > 0:
                            kmer_name = id_to_name[kid]
                            sample_name = self._sample_names[sid]
                            results[kmer_name][sample_name] = int(ab)
            else:
                # 密集矩阵：批量索引|Dense matrix: batch indexing
                abundance_matrix = f['abundance'][:]

                # 构建反向映射|Build reverse mapping
                id_to_name = {v: k for k, v in kmer_name_to_id.items()}

                for kmer_name, kmer_id in kmer_name_to_id.items():
                    abunds = abundance_matrix[kmer_id, :]
                    for i, ab in enumerate(abunds):
                        if ab > 0:
                            sample_name = self._sample_names[i]
                            results[kmer_name][sample_name] = int(ab)

        # 步骤3: 统计结果|Step 3: Statistics
        self.logger.info("步骤3/3: 统计结果|Step 3/3: Statistics")
        total_abundances = sum(len(v) for v in results.values())
        self.logger.info(f"批量查询完成|Batch query completed: {total_abundances} 个丰度值")

        return results

    def query_kmers_from_fasta(self, fasta_file: str) -> Dict[str, Dict[str, int]]:
        """从FASTA文件批量查询k-mer|Batch query k-mers from FASTA file

        Args:
            fasta_file: FASTA文件路径|FASTA file path

        Returns:
            kmer_name (from FASTA header) → {sample_name: abundance} 的映射
        """
        self.logger.info(f"从FASTA文件查询|Querying from FASTA: {fasta_file}")

        from Bio import SeqIO

        # 读取FASTA文件|Read FASTA file
        kmer_dict = {}
        for record in SeqIO.parse(fasta_file, 'fasta'):
            kmer_seq = str(record.seq).upper()  # 转大写
            kmer_name = record.id  # 使用FASTA header作为名称
            kmer_dict[kmer_name] = kmer_seq

        self.logger.info(f"从FASTA读取|Read from FASTA: {len(kmer_dict)} sequences")

        # 调用批量查询|Call batch query
        return self.query_kmers_batch(kmer_dict)

    def get_sample_abundances(self, sample_name: str) -> Dict[str, int]:
        """获取单个样本的所有k-mer丰度|Get all k-mer abundances for a single sample"""
        self._load_indices()

        self.logger.info(f"获取样本丰度|Getting sample abundances: {sample_name}")

        if sample_name not in self._sample_names:
            self.logger.error(f"样本不存在|Sample does not exist: {sample_name}")
            return {}

        sample_id = self._sample_names.index(sample_name)
        abundances = {}

        with h5py.File(self.matrix_file, 'r') as f:
            storage_type = f.attrs['storage_type'].decode('utf-8')

            if storage_type == 'sparse':
                # 稀疏存储|Sparse storage
                mask = f['sample_id'][:] == sample_id

                if mask.any():
                    kmer_ids = f['kmer_id'][:][mask]
                    abunds = f['abundance'][:][mask]

                    for kid, ab in zip(kmer_ids, abunds):
                        if ab > 0:
                            # 需要从k-mer ID反向查找k-mer序列|Need to reverse lookup k-mer from ID
                            kmer_seq = self._id_to_kmer(kid)
                            if kmer_seq:
                                abundances[kmer_seq] = int(ab)
            else:
                # 密集矩阵|Dense matrix
                with h5py.File(self.kmer_dict_file, 'r') as kf:
                    kmer_array = kf['kmer'][:]

                abunds = f['abundance'][:, sample_id]

                for i, ab in enumerate(abunds):
                    if ab > 0:
                        kmer_seq = kmer_array[i].decode('utf-8')
                        abundances[kmer_seq] = int(ab)

        self.logger.info(f"样本包含|Sample contains: {len(abundances)} 个k-mer")

        return abundances

    def _id_to_kmer(self, kmer_id: int) -> Optional[str]:
        """k-mer ID转序列|k-mer ID to sequence"""
        with h5py.File(self.kmer_dict_file, 'r') as f:
            if kmer_id < f.attrs['n_kmers']:
                return f['kmer'][kmer_id].decode('utf-8')
        return None

    def get_samples_with_kmer(self, kmer_seq: str, min_abundance: int = 1) -> List[str]:
        """获取包含特定k-mer的样本列表|Get list of samples containing specific k-mer"""
        abundances = self.query_kmer(kmer_seq)

        # 过滤最小丰度|Filter by minimum abundance
        samples = [sample for sample, ab in abundances.items() if ab >= min_abundance]

        return samples

    def get_kmer_count_per_sample(self) -> Dict[str, int]:
        """获取每个样本的k-mer数量|Get k-mer count per sample"""
        self._load_indices()

        counts = {}

        with h5py.File(self.matrix_file, 'r') as f:
            storage_type = f.attrs['storage_type'].decode('utf-8')

            if storage_type == 'sparse':
                # 稀疏存储|Sparse storage
                for i, sample_name in enumerate(self._sample_names):
                    mask = f['sample_id'][:] == i
                    counts[sample_name] = int(np.sum(mask))
            else:
                # 密集矩阵|Dense matrix
                for i, sample_name in enumerate(self._sample_names):
                    abunds = f['abundance'][:, i]
                    counts[sample_name] = int(np.sum(abunds > 0))

        return counts

    def export_sample_abundances_to_tsv(self, sample_name: str, output_file: str):
        """导出单个样本的丰度到TSV|Export single sample abundances to TSV"""
        abundances = self.get_sample_abundances(sample_name)

        with open(output_file, 'w') as f:
            f.write("kmer\tabundance\n")
            for kmer, abundance in sorted(abundances.items()):
                f.write(f"{kmer}\t{abundance}\n")

        self.logger.info(f"已导出|Exported: {output_file}")

    def export_matrix_to_tsv(self, output_file: str, max_kmers: Optional[int] = None):
        """导出完整矩阵到TSV|Export full matrix to TSV"""
        self._load_indices()

        self.logger.info("导出矩阵到TSV|Exporting matrix to TSV")

        with h5py.File(self.matrix_file, 'r') as f:
            storage_type = f.attrs['storage_type'].decode('utf-8')

            with open(output_file, 'w') as out:
                # 写入表头|Write header
                out.write("kmer\t" + "\t".join(self._sample_names) + "\n")

                if storage_type == 'sparse':
                    # 稀疏矩阵|Sparse matrix
                    n_kmers = f.attrs['n_kmers']

                    if max_kmers:
                        n_kmers = min(n_kmers, max_kmers)

                    # 为每个k-mer创建行|Create row for each k-mer
                    for kmer_id in range(n_kmers):
                        kmer_seq = self._id_to_kmer(kmer_id)
                        if not kmer_seq:
                            continue

                        # 查询此k-mer在所有样本中的丰度|Query this k-mer in all samples
                        mask = f['kmer_id'][:] == kmer_id

                        row_values = [kmer_seq]

                        if mask.any():
                            sample_ids = f['sample_id'][:][mask]
                            abunds = f['abundance'][:][mask]

                            # 构建完整的丰度行|Build full abundance row
                            sample_abunds = dict(zip(sample_ids, abunds))
                            for i in range(len(self._sample_names)):
                                row_values.append(str(sample_abunds.get(i, 0)))
                        else:
                            # 全0|All zeros
                            row_values.extend(['0'] * len(self._sample_names))

                        out.write("\t".join(row_values) + "\n")

                else:
                    # 密集矩阵|Dense matrix
                    with h5py.File(self.kmer_dict_file, 'r') as kf:
                        kmer_array = kf['kmer'][:]

                    n_kmers = f.attrs['n_kmers']
                    if max_kmers:
                        n_kmers = min(n_kmers, max_kmers)

                    for i in range(n_kmers):
                        kmer_seq = kmer_array[i].decode('utf-8')
                        abunds = f['abundance'][i, :]

                        row = [kmer_seq] + [str(ab) for ab in abunds]
                        out.write("\t".join(row) + "\n")

        self.logger.info(f"矩阵已导出|Matrix exported: {output_file}")

    def search_kmers_by_pattern(self, pattern: str) -> List[str]:
        """根据模式搜索k-mer|Search k-mers by pattern"""
        self._load_indices()

        import re

        matching_kmers = []

        with h5py.File(self.kmer_dict_file, 'r') as f:
            kmer_array = f['kmer'][:]

            for kmer_bytes in kmer_array:
                kmer = kmer_bytes.decode('utf-8')
                if re.search(pattern, kmer):
                    matching_kmers.append(kmer)

        self.logger.info(f"找到 {len(matching_kmers)} 个匹配的k-mer|Found {len(matching_kmers)} matching k-mers")

        return matching_kmers

    def close(self):
        """关闭数据库连接|Close database connection"""
        if self._db is not None:
            try:
                self._db.close()
                self._db = None
                db_name = 'RocksDB' if self._db_type == 'rocksdb' else 'LevelDB'
                self.logger.info(f"{db_name} 连接已关闭|{db_name} connection closed")
            except Exception as e:
                self.logger.warning(f"关闭数据库失败|Failed to close database: {e}")

    def __enter__(self):
        """上下文管理器入口|Context manager entry"""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """上下文管理器退出|Context manager exit"""
        self.close()
        return False

    def __del__(self):
        """析构函数|Destructor"""
        self.close()
