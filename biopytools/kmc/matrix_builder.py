"""
KMC丰度矩阵构建模块|KMC Abundance Matrix Builder Module
"""

import os
import shutil
import h5py
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple

from .config import KMCConfig
from .utils import KMCLogger, format_number


# ============================================================================
# RocksDB包装类（兼容plyvel API）|RocksDB Wrapper Class (plyvel API compatible)
# ============================================================================

def _get_system_memory_gb():
    """获取系统内存（GB）|Get system memory in GB"""
    try:
        import psutil
        return psutil.virtual_memory().total / (1024**3)
    except ImportError:
        # 如果没有psutil，使用默认值|If no psutil, use default
        return 800  # 假设800GB

def _calculate_rocksdb_config(max_memory_gb: float) -> dict:
    """根据最大内存限制计算RocksDB配置|Calculate RocksDB config based on max memory limit

    Args:
        max_memory_gb: 最大内存使用量（GB）|Maximum memory usage in GB

    Returns:
        配置字典|Config dict
    """
    # 直接使用指定的最大内存（已经是由用户或系统内存决定的值）
    # Use the specified max memory directly (already determined by user or system memory)
    available_memory_gb = max_memory_gb

    # 配置分配（按比例）|Configuration allocation (proportional)
    config = {
        'block_cache_size': int(available_memory_gb * 0.8 * 1024**3),  # 80%给块缓存|80% for block cache
        'write_buffer_size': int(available_memory_gb * 0.08 * 1024**3),  # 8%给MemTable|8% for MemTable
        'max_write_buffer_number': 4,  # MemTable数量（4个 = 32%总内存）|MemTable count (4 = 32% total)
        'block_cache_compressed_size': int(available_memory_gb * 0.1 * 1024**3),  # 10%给压缩缓存|10% for compressed cache
        'max_background_compactions': min(16, max(4, int(max_memory_gb / 50))),  # 每50GB内存1个线程，最多16个|1 thread per 50GB, max 16
        'max_background_flushes': 4,
        'max_bytes_for_level_base': int(available_memory_gb * 0.08 * 1024**3),  # 8%给Level 0|8% for Level 0
    }

    return config

class RocksDBAdapter:
    """RocksDB适配器，提供与plyvel兼容的API|RocksDB adapter with plyvel-compatible API"""

    def __init__(self, db_path: str, create_if_missing: bool = True,
                 error_if_exists: bool = False, block_cache_size: int = None,
                 max_memory_gb: float = None):
        """
        初始化RocksDB数据库|Initialize RocksDB database

        Args:
            db_path: 数据库路径|Database path
            create_if_missing: 如果不存在则创建|Create if not exists
            error_if_exists: 如果已存在则报错|Error if exists
            block_cache_size: 块缓存大小（字节），None表示自动配置|Block cache size in bytes, None means auto-configure
            max_memory_gb: 最大内存使用量（GB），None表示自动检测系统内存|Max memory usage in GB, None means auto-detect system memory
        """
        import rocksdb

        self.db_path = db_path
        self._rocksdb_module = rocksdb

        # 如果没有指定max_memory_gb，使用系统内存
        # If max_memory_gb not specified, use system memory
        if max_memory_gb is None:
            max_memory_gb = _get_system_memory_gb()

        # 如果没有指定block_cache_size，使用自动配置
        # If block_cache_size not specified, use auto-configuration
        if block_cache_size is None:
            config = _calculate_rocksdb_config(max_memory_gb)
            block_cache_size = config['block_cache_size']
            write_buffer_size = config['write_buffer_size']
            max_write_buffer_number = config['max_write_buffer_number']
            block_cache_compressed_size = config['block_cache_compressed_size']
            max_background_compactions = config['max_background_compactions']
            max_background_flushes = config['max_background_flushes']
            max_bytes_for_level_base = config['max_bytes_for_level_base']
            self.config_info = {
                'system_memory_gb': max_memory_gb,
                'block_cache_size_gb': block_cache_size / (1024**3),
                'write_buffer_size_gb': write_buffer_size / (1024**3),
                'max_write_buffer_number': max_write_buffer_number,
                'block_cache_compressed_size_gb': block_cache_compressed_size / (1024**3),
                'max_background_compactions': max_background_compactions,
                'max_background_flushes': max_background_flushes,
            }
        else:
            # 使用指定的block_cache_size，其他参数使用默认值
            # Use specified block_cache_size, other parameters use defaults
            write_buffer_size = 64 * 1024 * 1024 * 1024  # 64GB
            max_write_buffer_number = 4
            block_cache_compressed_size = 100 * 1024 * 1024 * 1024  # 100GB
            max_background_compactions = 16
            max_background_flushes = 4
            max_bytes_for_level_base = 64 * 1024 * 1024 * 1024  # 64GB
            self.config_info = {
                'system_memory_gb': max_memory_gb if max_memory_gb else _get_system_memory_gb(),
                'block_cache_size_gb': block_cache_size / (1024**3),
                'write_buffer_size_gb': write_buffer_size / (1024**3),
                'max_write_buffer_number': max_write_buffer_number,
                'block_cache_compressed_size_gb': block_cache_compressed_size / (1024**3),
                'max_background_compactions': max_background_compactions,
                'max_background_flushes': max_background_flushes,
            }

        # 配置RocksDB优化选项|Configure RocksDB options for performance
        opts = rocksdb.Options()
        opts.create_if_missing = create_if_missing

        # 配置表工厂（缓存、布隆过滤器）|Configure table factory (cache, bloom filter)
        opts.table_factory = rocksdb.BlockBasedTableFactory(
            block_cache=rocksdb.LRUCache(block_cache_size),
            filter_policy=rocksdb.BloomFilterPolicy(10),  # 10位布隆过滤器|10-bit bloom filter
            block_size=16 * 1024,  # 16KB块大小|16KB block size
            block_cache_compressed=rocksdb.LRUCache(block_cache_compressed_size)  # 压缩缓存|Compressed cache
        )

        # 配置MemTable（加速写入）|Configure MemTable (speed up writes)
        opts.write_buffer_size = write_buffer_size
        opts.max_write_buffer_number = max_write_buffer_number
        opts.min_write_buffer_number_to_merge = 1

        # 配置压缩（节省空间）|Configure compression (save space)
        opts.compression = rocksdb.CompressionType.lz4_compression  # 快速压缩|Fast compression

        # 配置并发优化|Configure parallelism
        opts.max_background_compactions = max_background_compactions  # 后台压缩线程|Background compaction threads
        opts.max_background_flushes = max_background_flushes  # 后台刷盘线程|Background flush threads

        # 配置层级大小|Configure level sizes
        opts.num_levels = 7
        opts.level0_file_num_compaction_trigger = 4
        opts.max_bytes_for_level_base = max_bytes_for_level_base

        # 检查数据库是否已存在|Check if database already exists
        if error_if_exists and os.path.exists(db_path) and os.listdir(db_path):
            raise rocksdb.Error(f"Database already exists: {db_path}")

        # 如果是新建且目录存在，先删除|If creating new and dir exists, delete first
        if create_if_missing and not error_if_exists:
            if os.path.exists(db_path):
                import shutil as _shutil
                _shutil.rmtree(db_path)

        # 创建数据库|Create database
        os.makedirs(db_path, exist_ok=True)
        self.db = rocksdb.DB(db_path, opts)

        # 保存配置信息（用于日志）|Save config info (for logging)
        self.config_info = {
            'system_memory_gb': max_memory_gb,
            'block_cache_size_gb': block_cache_size / (1024**3),
            'write_buffer_size_gb': write_buffer_size / (1024**3),
            'max_write_buffer_number': max_write_buffer_number,
            'block_cache_compressed_size_gb': block_cache_compressed_size / (1024**3),
            'max_background_compactions': max_background_compactions,
            'max_background_flushes': max_background_flushes,
            'max_bytes_for_level_base_gb': max_bytes_for_level_base / (1024**3),
        }

    def get(self, key: bytes, default=None):
        """获取值|Get value"""
        value = self.db.get(key)
        if value is None:
            return default
        return value

    def put(self, key: bytes, value: bytes):
        """写入键值对|Put key-value pair"""
        self.db.put(key, value)

    def delete(self, key: bytes):
        """删除键|Delete key"""
        self.db.delete(key)

    def write(self, batch):
        """写入批量|Write batch"""
        self.db.write(batch.batch)

    def write_batch(self):
        """创建批量写入对象|Create write batch object"""
        return self.WriteBatch(self._rocksdb_module)

    def close(self):
        """关闭数据库|Close database"""
        del self.db
        self.db = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    class WriteBatch:
        """批量写入对象（兼容plyvel）|Write batch object (plyvel compatible)"""

        def __init__(self, rocksdb_module):
            self._rocksdb = rocksdb_module
            self.batch = self._rocksdb.WriteBatch()

        def put(self, key: bytes, value: bytes):
            """添加到批量|Add to batch"""
            self.batch.put(key, value)

        def delete(self, key: bytes):
            """添加删除到批量|Add delete to batch"""
            self.batch.delete(key)

        def write(self):
            """执行批量写入|Execute batch write (需要传入db对象)"""
            # 这个方法不会被调用，实际调用的是RocksDBAdapter.write(batch)
            pass

        def clear(self):
            """清空批量|Clear batch"""
            self.batch = self._rocksdb.WriteBatch()


def DB(db_path: str, create_if_missing: bool = True,
        error_if_exists: bool = False, max_memory_gb: float = None, **kwargs):
    """兼容plyvel.DB的工厂函数|Factory function compatible with plyvel.DB

    Args:
        db_path: 数据库路径|Database path
        create_if_missing: 如果不存在则创建|Create if not exists
        error_if_exists: 如果已存在则报错|Error if exists
        max_memory_gb: 最大内存使用量（GB），None表示自动检测系统内存|Max memory usage in GB, None means auto-detect system memory
        **kwargs: 额外参数（忽略）|Additional parameters (ignored)

    Returns:
        RocksDBAdapter实例|RocksDBAdapter instance
    """
    return RocksDBAdapter(db_path, create_if_missing, error_if_exists,
                         max_memory_gb=max_memory_gb)


class KMCMatrixBuilder:
    """KMC丰度矩阵构建器|KMC Abundance Matrix Builder"""

    def __init__(self, config: KMCConfig, logger: Optional[KMCLogger] = None):
        """初始化矩阵构建器|Initialize matrix builder"""
        self.config = config
        self.config.validate()

        if logger is None:
            self.logger_manager = KMCLogger(self.config.output_path)
            self.logger = self.logger_manager.get_logger()
        else:
            self.logger_manager = logger
            self.logger = logger.get_logger()

        # HDF5文件路径|HDF5 file path
        self.matrix_file = self.config.output_path / 'abundance_matrix.h5'
        self.kmer_dict_file = self.config.output_path / 'kmer_dictionary.h5'
        self.metadata_file = self.config.output_path / 'metadata.json'

    def build_kmer_dictionary(self, sample_names: List[str]) -> Dict[str, int]:
        """构建全局k-mer字典|Build global k-mer dictionary"""
        self.logger.info("开始构建k-mer字典|Starting to build k-mer dictionary")

        # 首先合并所有样本数据库|First union all sample databases
        from .kmer_counter import KMCCounter
        counter = KMCCounter(self.config, self.logger_manager)

        self.logger.info("合并所有样本数据库|Union all sample databases")
        if not counter.union_databases(sample_names, 'global_kmers'):
            self.logger.error("合并数据库失败|Union databases failed")
            return {}

        # 导出全局k-mer列表|Export global k-mer list
        global_db = self.config.get_global_db_path()
        kmer_list_file = self.config.output_path / 'global_kmers.txt'

        cmd = f"{self.config.get_kmc_tools_bin()} transform {global_db} dump {kmer_list_file}"
        result = os.system(cmd)

        if result != 0:
            self.logger.error("导出k-mer列表失败|Export k-mer list failed")
            return {}

        # 读取k-mer列表|Read k-mer list
        kmer_dict = {}
        with open(kmer_list_file, 'r') as f:
            for i, line in enumerate(f):
                kmer = line.strip().split()[0]  # 格式: kmer count
                kmer_dict[kmer] = i

        self.logger.info(f"k-mer字典构建完成|K-mer dictionary built: {len(kmer_dict)} 个k-mer")

        # 保存k-mer字典到HDF5|Save k-mer dictionary to HDF5
        self._save_kmer_dictionary(kmer_dict)

        return kmer_dict

    def _save_kmer_dictionary(self, kmer_dict: Dict[str, int]):
        """保存k-mer字典到HDF5|Save k-mer dictionary to HDF5"""
        self.logger.info("保存k-mer字典|Saving k-mer dictionary")

        with h5py.File(self.kmer_dict_file, 'w') as f:
            # 创建k-mer字符串数据集|Create k-mer string dataset
            kmer_list = sorted(kmer_dict.keys(), key=lambda x: kmer_dict[x])
            kmer_array = np.array(kmer_list, dtype='S')

            f.create_dataset('kmer', data=kmer_array, compression='gzip', maxshape=(None,))
            f.attrs['n_kmers'] = len(kmer_dict)
            f.attrs['kmer_size'] = self.config.kmer_size

        self.logger.info(f"k-mer字典已保存|K-mer dictionary saved: {self.kmer_dict_file}")

    def _save_kmer_dictionary_from_file(self, kmer_file: str, n_kmers: int):
        """从文件保存k-mer字典到HDF5（流式处理）|Save k-mer dictionary from file to HDF5 (streaming)"""
        self.logger.info("从文件保存k-mer字典到HDF5|Saving k-mer dictionary from file to HDF5")

        kmer_list = []
        with open(kmer_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if parts:
                    kmer_list.append(parts[0])

        with h5py.File(self.kmer_dict_file, 'w') as f:
            kmer_array = np.array(kmer_list, dtype='S')
            f.create_dataset('kmer', data=kmer_array, compression='gzip', maxshape=(None,))
            f.attrs['n_kmers'] = n_kmers
            f.attrs['kmer_size'] = self.config.kmer_size

        self.logger.info(f"k-mer字典已保存|K-mer dictionary saved: {self.kmer_dict_file}")

    def _load_kmer_dictionary(self) -> Dict[str, int]:
        """加载k-mer字典|Load k-mer dictionary"""
        if not os.path.exists(self.kmer_dict_file):
            return {}

        kmer_dict = {}
        with h5py.File(self.kmer_dict_file, 'r') as f:
            kmer_array = f['kmer'][:]
            for i, kmer in enumerate(kmer_array):
                kmer_dict[kmer.decode('utf-8')] = i

        return kmer_dict

    def scan_existing_databases(self) -> List[str]:
        """扫描已存在的KMC数据库|Scan existing KMC databases"""
        self.logger.info("扫描已存在的KMC数据库|Scanning existing KMC databases")

        db_path = self.config.kmc_db_path
        if not db_path.exists():
            self.logger.error(f"数据库目录不存在|Database directory does not exist: {db_path}")
            return []

        # 查找所有 .kmc_pre 文件|Find all .kmc_pre files
        kmc_files = list(db_path.glob("*.kmc_pre"))

        if not kmc_files:
            self.logger.error(f"未找到KMC数据库文件|No KMC database files found in: {db_path}")
            return []

        # 提取样本名称，并过滤掉临时文件|Extract sample names and filter out temp files
        sample_names = []
        for kmc_file in kmc_files:
            sample_name = kmc_file.stem  # .kmc_pre文件去掉.kmc_pre后缀就是样本名
            # 过滤临时文件（以.temp_union_开头）|Filter out temporary files (starting with .temp_union_)
            if not sample_name.startswith('.temp_union_'):
                sample_names.append(sample_name)
                self.logger.debug(f"找到数据库|Found database: {sample_name}")

        self.logger.info(f"找到 {len(sample_names)} 个KMC数据库|Found {len(sample_names)} KMC databases")

        return sample_names

    def _load_kmer_size_from_count_metadata(self):
        """从count步骤保存的metadata读取kmer_size|Load kmer_size from count metadata"""
        import json

        metadata_file = self.config.output_path / 'kmer_metadata.json'

        if metadata_file.exists():
            try:
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)

                if 'kmer_size' in metadata and metadata['kmer_size'] > 0:
                    old_kmer_size = self.config.kmer_size
                    self.config.kmer_size = metadata['kmer_size']
                    self.logger.info(f"从count metadata读取k-mer大小|Read k-mer size from count metadata: "
                                   f"{old_kmer_size} → {self.config.kmer_size}")
            except Exception as e:
                self.logger.warning(f"读取count metadata失败，使用默认kmer_size|"
                                  f"Failed to read count metadata, using default kmer_size: {e}")
        else:
            self.logger.info(f"未找到count metadata，使用默认kmer_size|"
                           f"Count metadata not found, using default kmer_size: {self.config.kmer_size}")

    def build_matrix(self, sample_names: Optional[List[str]] = None) -> bool:
        """构建丰度矩阵|Build abundance matrix

        使用union+流式处理（速度快且内存友好）|Use union+streaming (fast and memory-friendly)
        支持断点续传|Support resuming from partial results
        """
        # 首先尝试从count步骤保存的metadata读取kmer_size|First try to read kmer_size from count metadata
        self._load_kmer_size_from_count_metadata()

        # 如果没有提供样本名，自动扫描|If no sample names provided, auto-scan
        if sample_names is None:
            self.logger.info("未指定样本名称，自动扫描数据库|No sample names specified, auto-scanning databases")
            sample_names = self.scan_existing_databases()
            if not sample_names:
                self.logger.error("未找到任何数据库，请先运行count命令|No databases found, please run count command first")
                return False
        else:
            self.logger.info(f"使用指定的 {len(sample_names)} 个样本|Using {len(sample_names)} specified samples")

        self.logger.info(f"开始构建丰度矩阵|Starting to build abundance matrix: {len(sample_names)} 样本|samples")

        # 检查断点续传|Check for resuming
        existing_samples = self._get_existing_samples()
        if existing_samples:
            self.logger.info(f"发现已存在的矩阵包含 {len(existing_samples)} 个样本|"
                           f"Found existing matrix with {len(existing_samples)} samples")
            self.logger.info(f"已处理样本|Already processed: {existing_samples}")

            # 过滤掉已处理的样本|Filter out already processed samples
            remaining_samples = [s for s in sample_names if s not in existing_samples]
            if not remaining_samples:
                self.logger.info("所有样本已处理，无需重新构建|All samples already processed, no need to rebuild")
                return True

            self.logger.info(f"需要处理剩余 {len(remaining_samples)} 个样本|"
                           f"Need to process {len(remaining_samples)} remaining samples")
            sample_names = remaining_samples

        n_samples = len(sample_names)

        # 检查是否已有k-mer文件（断点续传）|Check if k-mer files already exist (resume)
        kmer_with_id_file = self.config.output_path / 'global_kmers_with_id.txt'
        global_kmers_files_exist = kmer_with_id_file.exists()

        if global_kmers_files_exist:
            self.logger.info(f"检测到已存在的k-mer文件，跳过union步骤|Found existing k-mer file, skipping union: {kmer_with_id_file}")
            kmer_file = str(kmer_with_id_file)

            # 读取元数据获取k-mer总数|Read metadata to get total k-mer count
            import json
            metadata_file = self.config.output_path / 'kmer_metadata.json'
            if metadata_file.exists():
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                    n_kmers = metadata['n_kmers']

                # 如果n_kmers为0，说明元数据不完整，需要统计实际行数
                # If n_kmers is 0, metadata is incomplete, need to count actual lines
                if n_kmers == 0:
                    self.logger.info(f"元数据中n_kmers为0，统计k-mer数量（可能需要几分钟）|"
                                   f"n_kmers is 0 in metadata, counting k-mers (may take a few minutes)")
                    n_kmers = sum(1 for _ in open(kmer_file))
                    self.logger.info(f"统计k-mer数量：{n_kmers:,}|Counted k-mers: {n_kmers:,}")
                    # 更新元数据|Update metadata
                    metadata['n_kmers'] = n_kmers
                    with open(metadata_file, 'w') as f:
                        json.dump(metadata, f, indent=2)
                    self.logger.info(f"元数据已更新|Metadata updated")
                else:
                    self.logger.info(f"从元数据加载 {n_kmers:,} 个k-mer|Loaded {n_kmers:,} k-mers from metadata")
            else:
                # 如果没有元数据，需要统计行数并创建元数据|If no metadata, count lines and create metadata
                self.logger.info(f"元数据文件不存在，统计k-mer数量（可能需要几分钟）|Metadata file not found, counting k-mers (may take a few minutes)")
                n_kmers = sum(1 for _ in open(kmer_file))
                self.logger.info(f"统计k-mer数量：{n_kmers:,}|Counted k-mers: {n_kmers:,}")

                # 创建元数据文件|Create metadata file
                with open(metadata_file, 'w') as f:
                    json.dump({
                        'n_kmers': n_kmers,
                        'kmer_file': kmer_file,
                        'kmer_size': self.config.kmer_size,  # 保存k-mer大小|Save k-mer size
                        'index_type': None,  # 索引尚未构建
                        'index_path': None
                    }, f, indent=2)
                self.logger.info(f"元数据文件已创建|Metadata file created: {metadata_file}")

            # 初始化HDF5矩阵文件|Initialize HDF5 matrix file
            self._initialize_matrix(n_kmers, n_samples, sample_names)
        elif not existing_samples:
            self.logger.info("全新构建：使用union生成全局k-mer列表|Fresh build: using union to generate global k-mer list")
            kmer_file = self._collect_all_kmers(sample_names)
            if not kmer_file:
                self.logger.error("收集k-mer失败|Failed to collect k-mers")
                return False

            # 读取元数据获取k-mer总数和k-mer大小|Read metadata to get total k-mer count and k-mer size
            import json
            with open(self.config.output_path / 'kmer_metadata.json', 'r') as f:
                metadata = json.load(f)
                n_kmers = metadata['n_kmers']
                # 从元数据读取k-mer大小并更新config（matrix模式可能没有指定-k参数）|Read k-mer size from metadata and update config
                if 'kmer_size' in metadata:
                    self.config.kmer_size = metadata['kmer_size']
                    self.logger.info(f"从元数据读取k-mer大小|Read k-mer size from metadata: {self.config.kmer_size}")
                else:
                    self.logger.warning(f"元数据中缺少kmer_size字段，使用默认值|kmer_size not in metadata, using default: {self.config.kmer_size}")

            self.logger.info(f"共收集到 {n_kmers} 个唯一k-mer|Collected {n_kmers} unique k-mers")

            # 初始化HDF5矩阵文件|Initialize HDF5 matrix file
            self._initialize_matrix(n_kmers, n_samples, sample_names)
        else:
            # 从现有矩阵加载k-mer信息|Load k-mer info from existing matrix
            kmer_metadata_file = self.config.output_path / 'kmer_metadata.json'
            if not kmer_metadata_file.exists():
                self.logger.error("未找到k-mer元数据文件，无法断点续传|K-mer metadata file not found, cannot resume")
                return False

            with open(kmer_metadata_file, 'r') as f:
                import json
                metadata = json.load(f)
                kmer_file = metadata['kmer_file']
                n_kmers = metadata['n_kmers']
                # 从元数据读取k-mer大小并更新config|Read k-mer size from metadata and update config
                if 'kmer_size' in metadata:
                    self.config.kmer_size = metadata['kmer_size']

            self.logger.info(f"从元数据加载 {n_kmers} 个k-mer|Loaded {n_kmers} k-mers from metadata")

            # 扩展矩阵|Extend matrix
            self._extend_matrix_file(n_kmers, sample_names, existing_samples)

        # 获取k-mer文件大小|Get k-mer file size
        kmer_file_size_gb = os.path.getsize(kmer_file) / (1024**3)
        self.logger.info(f"k-mer文件大小|K-mer file size: {kmer_file_size_gb:.2f}GB")

        # 构建或加载索引（自动选择内存或数据库）|Build or load index (auto select memory or DB)
        self.logger.info("加载或构建索引|Loading or building index")
        kmer_index, index_type = self._load_or_create_index(kmer_file, kmer_file_size_gb)

        # 查询每个样本的丰度|Query abundance for each sample
        for i, sample_name in enumerate(sample_names, 1):
            self.logger.info(f"处理样本|Processing sample {i}/{len(sample_names)}: {sample_name}")

            abundances = self._query_sample_abundances(sample_name, kmer_index, index_type)

            # 保存到矩阵|Save to matrix
            sample_idx = len(existing_samples) + i - 1 if existing_samples else i - 1
            self._save_sample_abundances(sample_name, sample_idx, abundances)

        # 关闭数据库索引（如果是RocksDB）|Close database index (if RocksDB)
        if index_type == 'rocksdb':
            kmer_index.close()
            self.logger.info("RocksDB索引已关闭|RocksDB index closed")

        # 保存k-mer列表到HDF5（用于导出）|Save k-mer list to HDF5 (for export)
        self._save_kmer_dictionary_from_file(kmer_file, n_kmers)

        self.logger.info(f"丰度矩阵构建完成|Abundance matrix built: {self.matrix_file}")

        return True

    def _get_existing_samples(self) -> List[str]:
        """获取已处理的样本列表|Get list of already processed samples"""
        if not self.matrix_file.exists():
            return []

        try:
            with h5py.File(self.matrix_file, 'r') as f:
                return [s.decode('utf-8') for s in f['sample_names'][:]]
        except Exception as e:
            self.logger.warning(f"读取现有矩阵失败|Failed to read existing matrix: {e}")
            return []

    def _load_existing_kmers(self) -> List[str]:
        """从现有矩阵加载k-mer列表|Load k-mer list from existing matrix"""
        if not self.kmer_dict_file.exists():
            return []

        try:
            with h5py.File(self.kmer_dict_file, 'r') as f:
                return [k.decode('utf-8') for k in f['kmer'][:]]
        except Exception as e:
            self.logger.warning(f"读取k-mer字典失败|Failed to read k-mer dictionary: {e}")
            return []

    def _collect_all_kmers(self, sample_names: List[str]) -> Optional[str]:
        """收集所有样本的唯一k-mer（使用union+流式处理，速度快且内存友好）|Collect all unique k-mers using union+streaming"""
        self.logger.info("收集唯一k-mer列表（union+流式处理）|Collecting unique k-mer list (union+streaming)")

        # 使用union生成全局k-mer列表（速度快）|Use union to generate global k-mer list (fast)
        from .kmer_counter import KMCCounter
        counter = KMCCounter(self.config, self.logger_manager)

        self.logger.info("合并数据库生成全局k-mer列表|Union databases to generate global k-mer list")
        if not counter.union_databases(sample_names, 'global_kmers'):
            self.logger.error("合并数据库失败|Union databases failed")
            return None

        # 导出k-mer列表|Export k-mer list
        global_db = self.config.output_path / 'global_kmers'
        kmer_list_file = self.config.output_path / 'global_kmers.txt'

        self.logger.info("导出全局k-mer列表|Exporting global k-mer list")
        cmd = f"{self.config.get_kmc_tools_bin()} transform {global_db} dump {kmer_list_file}"
        result = os.system(cmd)

        if result != 0:
            self.logger.error("导出k-mer列表失败|Export k-mer list failed")
            return None

        # 流式添加kmer_id列（不加载到内存）|Stream and add kmer_id column (not load to memory)
        kmer_with_id_file = self.config.output_path / 'global_kmers_with_id.txt'
        self._add_kmer_id_column(kmer_list_file, kmer_with_id_file)

        # 统计k-mer数量|Count k-mers
        n_kmers = sum(1 for _ in open(kmer_list_file))

        # 创建元数据文件|Create metadata file
        import json
        metadata_file = self.config.output_path / 'kmer_metadata.json'
        metadata_content = {
            'n_kmers': n_kmers,
            'kmer_file': str(kmer_list_file),
            'kmer_size': self.config.kmer_size,  # 保存k-mer大小|Save k-mer size
            'index_type': None,  # 索引尚未构建
            'index_path': None
        }
        with open(metadata_file, 'w') as f:
            json.dump(metadata_content, f, indent=2)
        self.logger.info(f"元数据文件已创建|Metadata file created: {metadata_file}")
        self.logger.info(f"元数据内容|Metadata content: kmer_size={metadata_content['kmer_size']}, n_kmers={n_kmers}")

        # 保留中间文件用于断点续传|Keep intermediate files for resuming
        # - global_kmers.txt: 原始导出文件
        # - global_kmers.kmc_pre/suf: union后的KMC数据库
        # - global_kmers_with_id.txt: 添加ID后的文件（用于后续处理）
        self.logger.info("保留中间文件用于断点续传|Keeping intermediate files for resuming")

        self.logger.info(f"k-mer列表已生成并添加ID列|K-mer list generated with ID column: {kmer_with_id_file}")

        return str(kmer_with_id_file)

    def _add_kmer_id_column(self, input_file: Path, output_file: Path):
        """流式添加kmer_id列（不加载到内存）|Stream and add kmer_id column (not load to memory)

        输入格式|Input format:  kmer count
        输出格式|Output format: kmer count kmer_id
        """
        self.logger.info(f"添加kmer_id列（流式处理）|Adding kmer_id column (streaming): {input_file}")

        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            for kmer_id, line in enumerate(fin):
                fout.write(line.strip() + f'\t{kmer_id}\n')

        # 记录总k-mer数量|Log total k-mer count
        kmer_count = kmer_id + 1
        self.logger.info(f"k-mer ID列添加完成|K-mer ID column added: {kmer_count} k-mers")

        # 保存元数据|Save metadata
        with open(self.config.output_path / 'kmer_metadata.json', 'w') as f:
            import json
            json.dump({
                'n_kmers': kmer_count,
                'kmer_file': str(output_file),
                'kmer_size': self.config.kmer_size  # 保存k-mer大小|Save k-mer size
            }, f)

    def _should_use_db_index(self, kmer_file_size_gb: float) -> bool:
        """判断是否使用数据库索引|Decide whether to use database index

        Args:
            kmer_file_size_gb: k-mer文件大小（GB）|K-mer file size in GB

        Returns:
            True表示使用数据库，False表示使用内存|True for DB, False for memory
        """
        index_mode = self.config.index_mode.lower()

        if index_mode == 'memory':
            self.logger.info("强制使用内存索引|Forced to use memory index")
            return False
        elif index_mode == 'db':
            self.logger.info("强制使用数据库索引|Forced to use database index")
            return True
        else:  # auto
            threshold = self.config.index_threshold_gb
            if kmer_file_size_gb >= threshold:
                self.logger.info(f"文件大小 {kmer_file_size_gb:.2f}GB >= 阈值 {threshold}GB，使用数据库索引|"
                               f"File size {kmer_file_size_gb:.2f}GB >= threshold {threshold}GB, using database index")
                return True
            else:
                self.logger.info(f"文件大小 {kmer_file_size_gb:.2f}GB < 阈值 {threshold}GB，使用内存索引|"
                               f"File size {kmer_file_size_gb:.2f}GB < threshold {threshold}GB, using memory index")
                return False

    def _build_rocksdb_index(self, kmer_file: str) -> str:
        """构建RocksDB索引|Build RocksDB index

        Args:
            kmer_file: k-mer文件路径（包含ID列）|K-mer file path (with ID column)

        Returns:
            RocksDB数据库路径|RocksDB database path
        """
        # 使用RocksDB（通过兼容层）|Use RocksDB (via compatibility layer)
        db_path = self.config.output_path / 'kmer_index.rdb'

        self.logger.info(f"构建RocksDB索引|Building RocksDB index: {db_path}")

        # 创建数据库（使用配置的最大内存）|Create database (use configured max memory)
        db = DB(str(db_path), create_if_missing=True, error_if_exists=False,
                max_memory_gb=self.config.max_memory)  # 使用配置的最大内存|Use configured max memory

        # 记录实际配置|Log actual configuration
        config = db.config_info
        self.logger.info(
            f"最大内存|Max memory: {config['system_memory_gb']:.0f} GB\n"
            f"RocksDB配置|RocksDB config:\n"
            f"  块缓存|Block cache: {config['block_cache_size_gb']:.0f} GB\n"
            f"  MemTable: {config['write_buffer_size_gb']:.0f} GB × {config['max_write_buffer_number']} = {config['write_buffer_size_gb'] * config['max_write_buffer_number']:.0f} GB\n"
            f"  压缩缓存|Compressed cache: {config['block_cache_compressed_size_gb']:.0f} GB\n"
            f"  后台线程|Background threads: {config['max_background_compactions']} (压缩) + {config['max_background_flushes']} (刷盘)\n"
            f"  布隆过滤器|Bloom filter: 10位|10-bit\n"
            f"  压缩算法|Compression: LZ4"
        )

        # 批量写入（控制内存占用）|Batch write (control memory usage)
        batch_size = 100000000  # 1亿条记录/100M records (~5-10GB memory)
        current_batch = []
        total_kmers = 0

        with open(kmer_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    kmer = parts[0]
                    # 第2列是kmer_id，第3列是count（从文件内容确认）|Column 2 is kmer_id, column 3 is count (confirmed from file)
                    kmer_id = parts[1].encode()
                    current_batch.append((kmer.encode(), kmer_id))
                    total_kmers += 1

                    # 批量写入|Batch write
                    if len(current_batch) >= batch_size:
                        wb = db.write_batch()
                        for k, v in current_batch:
                            wb.put(k, v)
                        db.write(wb)  # 调用RocksDBAdapter.write()|Call RocksDBAdapter.write()
                        current_batch = []

                        # 进度日志|Progress log
                        if total_kmers % 10000000 == 0:  # 每1000万条记录记录一次|Log every 10M
                            self.logger.info(f"  已处理|Processed: {total_kmers:,} k-mers")

            # 写入剩余数据|Write remaining data
            if current_batch:
                wb = db.write_batch()
                for k, v in current_batch:
                    wb.put(k, v)
                db.write(wb)  # 调用RocksDBAdapter.write()|Call RocksDBAdapter.write()

        db.close()

        self.logger.info(f"RocksDB索引构建完成：{total_kmers:,} k-mers|RocksDB index built: {total_kmers:,} k-mers")

        return str(db_path)

    def _load_or_create_index(self, kmer_file: str, kmer_file_size_gb: float):
        """加载或创建索引（自动选择）|Load or create index (auto select)

        Args:
            kmer_file: k-mer文件路径|K-mer file path
            kmer_file_size_gb: k-mer文件大小（GB）|K-mer file size in GB

        Returns:
            (索引对象, 索引类型)|Index object, index type
        """
        # 判断使用哪种索引|Decide which index to use
        use_db = self._should_use_db_index(kmer_file_size_gb)

        if use_db:
            # 使用RocksDB|Use RocksDB
            db_path = self.config.output_path / 'kmer_index.rdb'

            # 检查数据库是否已存在|Check if database already exists
            if db_path.exists() and any(db_path.iterdir()):  # 目录存在且非空
                self.logger.info(f"检测到已存在的RocksDB索引，直接加载|Found existing RocksDB index, loading: {db_path}")
                try:
                    db = DB(str(db_path), create_if_missing=False,
                           max_memory_gb=self.config.max_memory)  # 使用配置的最大内存|Use configured max memory
                    # 更新元数据|Update metadata
                    metadata_file = self.config.output_path / 'kmer_metadata.json'
                    import json
                    if metadata_file.exists():
                        with open(metadata_file, 'r') as f:
                            metadata = json.load(f)
                        metadata['index_type'] = 'rocksdb'
                        metadata['index_path'] = str(db_path)
                        with open(metadata_file, 'w') as f:
                            json.dump(metadata, f, indent=2)
                    return db, 'rocksdb'
                except Exception as e:
                    self.logger.warning(f"打开现有RocksDB索引失败，将重建: {e}")

            # 需要创建新的索引|Need to create new index
            db_path_str = self._build_rocksdb_index(kmer_file)
            db = DB(db_path_str, create_if_missing=False,
                   max_memory_gb=self.config.max_memory)  # 使用配置的最大内存|Use configured max memory

            # 更新元数据|Update metadata
            metadata_file = self.config.output_path / 'kmer_metadata.json'
            import json
            metadata = {}
            if metadata_file.exists():
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
            metadata['index_type'] = 'rocksdb'
            metadata['index_path'] = db_path_str
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)

            return db, 'rocksdb'
        else:
            # 使用内存索引|Use memory index
            return self._build_memory_index(kmer_file), 'memory'

    def _build_memory_index(self, kmer_file: str) -> Dict[str, int]:
        """构建内存索引|Build in-memory index

        注意：这个方法会占用大量内存，建议用于小规模数据
        Note: This method uses lots of memory, recommended for small data only
        """
        self.logger.warning("构建内存索引（可能占用大量内存）|Building in-memory index (may use lots of memory)")

        kmer_to_id = {}
        with open(kmer_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    kmer = parts[0]
                    # 第2列是kmer_id，第3列是count（从文件内容确认）|Column 2 is kmer_id, column 3 is count (confirmed from file)
                    kmer_id = int(parts[1])
                    kmer_to_id[kmer] = kmer_id

        self.logger.info(f"内存索引构建完成：{len(kmer_to_id):,} k-mers|Memory index built: {len(kmer_to_id):,} k-mers")
        return kmer_to_id

    def _build_kmer_lookup_index(self, kmer_file: str) -> Dict[str, int]:
        """构建k-mer查找索引（用于快速查询）|Build k-mer lookup index (for fast query）

        注意：这个方法会占用内存，建议用于小规模数据或使用数据库索引
        Note: This method uses memory, recommended for small data or use database index
        """
        self.logger.warning("构建内存索引（可能占用大量内存）|Building in-memory index (may use lots of memory)")

        kmer_to_id = {}
        with open(kmer_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    kmer = parts[0]
                    # 第2列是kmer_id，第3列是count（从文件内容确认）|Column 2 is kmer_id, column 3 is count (confirmed from file)
                    kmer_id = int(parts[1])
                    kmer_to_id[kmer] = kmer_id

        self.logger.info(f"索引构建完成：{len(kmer_to_id)} 个k-mer|Index built: {len(kmer_to_id)} k-mers")
        return kmer_to_id

    def _save_kmer_set(self, kmer_set: set, filename: str):
        """保存k-mer集合到临时文件|Save k-mer set to temp file"""
        temp_file = self.config.output_path / filename
        with open(temp_file, 'w') as f:
            for kmer in sorted(kmer_set):
                f.write(f"{kmer}\n")

    def _setup_matrix_file(self, all_kmers: List[str], new_samples: List[str],
                           existing_samples: List[str]):
        """设置矩阵文件（新建或扩展）|Setup matrix file (create or extend)"""
        n_kmers = len(all_kmers)
        n_existing = len(existing_samples)
        n_new = len(new_samples)
        n_total = n_existing + n_new

        if not existing_samples:
            # 全新创建|Create new
            self.logger.info("创建新矩阵文件|Creating new matrix file")
            self._initialize_matrix(n_kmers, n_new, new_samples)
        else:
            # 扩展现有矩阵|Extend existing matrix
            self.logger.info(f"扩展现有矩阵: {n_existing} -> {n_total} 个样本|"
                           f"Extend existing matrix: {n_existing} -> {n_total} samples")
            self._extend_matrix_file(n_kmers, new_samples, existing_samples)

    def _extend_matrix_file(self, n_kmers: int, new_samples: List[str],
                            existing_samples: List[str]):
        """扩展矩阵文件|Extend matrix file"""
        try:
            with h5py.File(self.matrix_file, 'a') as f:
                storage_type = f.attrs['storage_type']
                # 兼容字符串和bytes|Compatible with both str and bytes
                if isinstance(storage_type, bytes):
                    storage_type_str = storage_type.decode('utf-8')
                else:
                    storage_type_str = storage_type

                if storage_type_str == 'dense':
                    # 扩展密集矩阵|Extend dense matrix
                    old_shape = f['abundance'].shape
                    new_shape = (n_kmers, len(existing_samples) + len(new_samples))

                    # 创建新的数据集|Create new dataset
                    new_data = np.zeros(new_shape, dtype=np.uint16)
                    new_data[:old_shape[0], :old_shape[1]] = f['abundance'][:]

                    del f['abundance']
                    f.create_dataset('abundance', data=new_data, compression='gzip', chunks=True)
                else:
                    # 稀疏矩阵：不需要扩展结构，只需更新样本列表|Sparse: no need to extend structure
                    pass

                # 更新样本列表|Update sample names
                all_samples = existing_samples + new_samples
                del f['sample_names']
                f.create_dataset('sample_names', data=np.array(all_samples, dtype='S'))

                # 更新元数据|Update metadata
                f.attrs['n_kmers'] = n_kmers
                f.attrs['n_samples'] = len(all_samples)

        except Exception as e:
            self.logger.error(f"扩展矩阵文件失败|Failed to extend matrix file: {e}")
            raise

    def _initialize_matrix(self, n_kmers: int, n_samples: int, sample_names: List[str]):
        """初始化HDF5矩阵|Initialize HDF5 matrix"""
        self.logger.info("初始化矩阵文件|Initializing matrix file")

        # 删除旧矩阵文件（确保全新构建）|Delete old matrix file (ensure fresh build)
        if self.matrix_file.exists():
            self.logger.info(f"删除旧矩阵文件|Deleting old matrix file: {self.matrix_file}")
            self.matrix_file.unlink()

        # 删除旧的k-mer字典文件|Delete old k-mer dictionary file
        if self.kmer_dict_file.exists():
            self.logger.info(f"删除旧k-mer字典文件|Deleting old k-mer dictionary file: {self.kmer_dict_file}")
            self.kmer_dict_file.unlink()

        with h5py.File(self.matrix_file, 'w') as f:
            # 创建稀疏存储的数据集|Create sparse storage dataset
            # 使用COO格式: (kmer_id, sample_id, abundance)|Use COO format
            if self.config.sparse_storage:
                # 预估非零元素数量|Estimate non-zero elements
                # 假设稀疏度为5%|Assume 5% sparsity
                estimated_nnz = int(n_kmers * n_samples * 0.05)

                f.create_dataset('kmer_id', shape=(estimated_nnz,), dtype=np.uint64,
                                maxshape=(None,), compression='gzip')
                f.create_dataset('sample_id', shape=(estimated_nnz,), dtype=np.uint16,
                                maxshape=(None,), compression='gzip')
                f.create_dataset('abundance', shape=(estimated_nnz,), dtype=np.uint16,
                                maxshape=(None,), compression='gzip')
                f.attrs['storage_type'] = b'sparse'
            else:
                # 密集矩阵|Dense matrix
                f.create_dataset('abundance', shape=(n_kmers, n_samples), dtype=np.uint16,
                                compression='gzip', chunks=True)
                f.attrs['storage_type'] = b'dense'

            # 元数据|Metadata
            f.attrs['n_kmers'] = n_kmers
            f.attrs['n_samples'] = n_samples
            f.attrs['kmer_size'] = self.config.kmer_size
            f.create_dataset('sample_names', data=np.array(sample_names, dtype='S'))

        self.logger.info("矩阵文件初始化完成|Matrix file initialized")

    def _query_sample_abundances(self, sample_name: str, kmer_index, index_type: str) -> Dict[int, int]:
        """查询单个样本的丰度|Query single sample abundances

        Args:
            sample_name: 样本名称|Sample name
            kmer_index: 索引对象（字典或RocksDB）|Index object (dict or RocksDB)
            index_type: 索引类型 ('memory' 或 'rocksdb')|Index type
        """
        self.logger.info(f"查询样本丰度|Querying sample abundances: {sample_name}")

        # 使用kmc_tools transform dump查询|Use kmc_tools transform dump to query
        db_path = self.config.kmc_db_path / sample_name
        dump_file = self.config.dump_path / f'{sample_name}_dump.txt'

        # 检查是否已有dump文件（断点续传）|Check if dump file already exists (resume)
        if dump_file.exists():
            self.logger.info(f"使用已有dump文件|Using existing dump file: {dump_file}")
        else:
            # 生成dump文件|Generate dump file
            cmd = f"{self.config.get_kmc_tools_bin()} transform {db_path} dump {dump_file}"
            result = os.system(cmd)

            if result != 0:
                self.logger.error(f"查询样本|Query sample {sample_name} 失败|failed")
                return {}

        # 读取dump文件并构建丰度字典|Read dump file and build abundance dictionary
        # 优化：先收集所有k-mer，然后批量查询|Optimize: collect all k-mers, then batch query
        kmer_count_pairs = []
        with open(dump_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    kmer = parts[0]
                    count = int(parts[1])
                    kmer_count_pairs.append((kmer, count))

        self.logger.info(f"从dump文件读取|Read from dump: {len(kmer_count_pairs):,} k-mers")

        # 批量查询|Batch query
        abundances = {}

        if index_type == 'memory':
            # 内存字典查找（很快，无需优化）|Memory dict lookup (fast, no need to optimize)
            for kmer, count in kmer_count_pairs:
                if kmer in kmer_index:
                    kmer_id = kmer_index[kmer]
                    abundances[kmer_id] = count
        else:  # rocksdb
            # RocksDB批量查询优化（使用多线程）|RocksDB batch query optimization (using multithreading)
            # RocksDB是线程安全的，多线程可以共享同一个数据库连接|RocksDB is thread-safe, multiple threads can share the same db connection
            if len(kmer_count_pairs) > 100000:  # 超过10万条使用并行|Use parallel for >100k
                abundances = self._batch_query_rocksdb_parallel(kmer_index, kmer_count_pairs)
            else:
                # 少量数据直接查询|Direct query for small data
                for kmer, count in kmer_count_pairs:
                    kid_bytes = kmer_index.get(kmer.encode())
                    if kid_bytes:
                        kmer_id = int(kid_bytes)
                        abundances[kmer_id] = count

        # 根据配置决定是否删除dump文件|Delete or keep dump file based on configuration
        if not self.config.keep_dump:
            os.remove(dump_file)
            self.logger.debug(f"已删除dump文件|Dump file deleted: {dump_file}")
        else:
            self.logger.debug(f"保留dump文件|Dump file kept: {dump_file}")

        self.logger.info(f"查询完成|Query completed: {len(abundances):,} 个k-mer")

        return abundances

    def _batch_query_rocksdb_parallel(self, db, kmer_count_pairs, num_processes=None):
        """并行批量查询RocksDB（使用多线程）|Parallel batch query RocksDB (using multithreading)

        RocksDB是线程安全的，可以使用多线程共享同一个数据库连接。
        比多进程更高效，因为不需要序列化数据和多次打开数据库。

        Args:
            db: RocksDB数据库对象|RocksDB database object
            kmer_count_pairs: (kmer, count) 列表|(kmer, count) list
            num_processes: 线程数|Number of threads (None表示使用config.threads)|None means use config.threads

        Returns:
            {kmer_id: count} 字典|{kmer_id: count} dict
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed

        if num_processes is None:
            # 使用配置的线程数|Use configured threads
            num_processes = self.config.threads

        self.logger.info(f"使用 {num_processes} 线程并行查询|Using {num_processes} threads for parallel query")

        # 将数据分块|Split data into chunks
        chunk_size = len(kmer_count_pairs) // num_processes + 1
        chunks = []
        for i in range(0, len(kmer_count_pairs), chunk_size):
            chunks.append(kmer_count_pairs[i:i+chunk_size])

        self.logger.info(f"数据分块|Data chunks: {len(chunks)}, 每块约|~{chunk_size:,} 条")

        # 多线程并行查询|Multithreaded parallel query
        # RocksDB是线程安全的，多个线程可以共享同一个db对象|RocksDB is thread-safe, multiple threads can share the same db object
        chunk_args = [(chunk, db) for chunk in chunks]

        abundances = {}
        with ThreadPoolExecutor(max_workers=num_processes) as executor:
            # 使用submit+as_completed以获得更好的错误处理|Use submit+as_completed for better error handling
            future_to_chunk = {executor.submit(_query_thread_worker, args): args[0] for args in chunk_args}

            for future in as_completed(future_to_chunk):
                chunk = future_to_chunk[future]
                try:
                    result = future.result()
                    abundances.update(result)
                except Exception as e:
                    self.logger.error(f"查询块出错|Error querying chunk: {e}")
                    raise

        return abundances

    def _save_sample_abundances(self, sample_name: str, sample_id: int,
                                abundances: Dict[int, int]):
        """保存样本丰度|Save sample abundances"""
        self.logger.info(f"保存样本丰度|Saving sample abundances: {sample_name}")

        with h5py.File(self.matrix_file, 'a') as f:
            storage_type = f.attrs['storage_type']
            # 兼容字符串和bytes|Compatible with both str and bytes
            if isinstance(storage_type, bytes):
                storage_type_str = storage_type.decode('utf-8')
            else:
                storage_type_str = storage_type
            self.logger.info(f"存储类型|Storage type: {storage_type_str} (类型|type: {type(storage_type)})")

            if storage_type_str == 'sparse':
                # 稀疏存储|Sparse storage
                self.logger.info(f"使用稀疏存储|Using sparse storage")
                # 需要追加数据|Need to append data
                current_size = f['kmer_id'].shape[0]
                self.logger.info(f"当前数据集大小|Current dataset size: {current_size}, 新增|new: {len(abundances)}")

                # 如果空间不足，扩展数据集|If not enough space, extend dataset
                new_size = current_size + len(abundances)
                for ds_name in ['kmer_id', 'sample_id', 'abundance']:
                    f[ds_name].resize((new_size,))
                self.logger.info(f"数据集已扩展|Dataset resized to: {new_size}")

                # 写入数据（批量写入，比逐个写入快100倍）|Write data (batch write, 100x faster)
                kmer_ids = np.array(list(abundances.keys()), dtype=np.uint64)
                sample_ids = np.full(len(abundances), sample_id, dtype=np.uint16)
                abunds = np.array(list(abundances.values()), dtype=np.uint16)

                f['kmer_id'][current_size:new_size] = kmer_ids
                f['sample_id'][current_size:new_size] = sample_ids
                f['abundance'][current_size:new_size] = abunds

            else:
                # 密集矩阵|Dense matrix
                self.logger.info(f"使用密集矩阵存储|Using dense matrix storage")
                self.logger.info(f"abundance 数据集形状|abundance dataset shape: {f['abundance'].shape}")
                for kmer_id, abundance in abundances.items():
                    f['abundance'][kmer_id, sample_id] = abundance

        self.logger.info(f"样本丰度已保存|Sample abundances saved: {sample_name}")

    def add_samples(self, new_samples: List[str], new_sample_names: List[str]) -> bool:
        """添加新样本到矩阵|Add new samples to matrix"""
        self.logger.info(f"添加 {len(new_samples)} 个新样本|Adding {len(new_samples)} new samples")

        # 1. 为新样本建立KMC数据库|Build KMC databases for new samples
        from .kmer_counter import KMCCounter
        counter = KMCCounter(self.config, self.logger_manager)

        # 更新配置|Update config
        self.config.input_files = new_samples
        self.config.sample_names = new_sample_names

        results = counter.count_samples(new_samples, new_sample_names)

        if not all(results.values()):
            self.logger.error("部分样本统计失败|Some samples counting failed")
            return False

        # 2. 加载现有k-mer字典|Load existing k-mer dictionary
        old_kmer_dict = self._load_kmer_dictionary()
        old_n_kmers = len(old_kmer_dict)

        # 3. 构建新的全局k-mer字典|Build new global k-mer dictionary
        # 获取现有样本名称|Get existing sample names
        with h5py.File(self.matrix_file, 'r') as f:
            existing_samples = [s.decode('utf-8') for s in f['sample_names'][:]]

        all_samples = existing_samples + new_sample_names

        # 重新构建字典|Rebuild dictionary
        new_kmer_dict = self.build_kmer_dictionary(all_samples)

        new_n_kmers = len(new_kmer_dict)
        added_kmers = new_n_kmers - old_n_kmers

        self.logger.info(f"新增k-mer数量|New k-mers added: {format_number(added_kmers)}")

        # 4. 更新矩阵|Update matrix
        if added_kmers > 0:
            # 需要扩展矩阵|Need to extend matrix
            self._extend_matrix(new_n_kmers, all_samples)

        # 5. 查询新样本的丰度|Query new sample abundances
        for i, sample_name in enumerate(new_sample_names, len(existing_samples)):
            self.logger.info(f"处理新样本|Processing new sample: {sample_name}")

            abundances = self._query_sample_abundances(sample_name, new_kmer_dict)
            self._save_sample_abundances(sample_name, i, abundances)

        # 6. 如果有新k-mer，查询旧样本的丰度|If new k-mers, query old sample abundances
        if added_kmers > 0:
            self._query_old_samples_for_new_kmers(existing_samples, old_kmer_dict, new_kmer_dict)

        self.logger.info("新样本添加完成|New samples added")

        return True

    def add_samples_incremental(self, new_samples: List[str], new_sample_names: List[str]) -> bool:
        """增量式添加新样本到矩阵|Incrementally add new samples to matrix

        优化版：不重新处理旧样本，只添加新样本和新 k-mer
        Optimized: Do not reprocess old samples, only add new samples and new k-mers

        Args:
            new_samples: 新样本的FASTQ文件列表|New sample FASTQ files
            new_sample_names: 新样本名称列表|New sample names

        Returns:
            是否成功|Success or not
        """
        self.logger.info(f"增量式添加 {len(new_samples)} 个新样本|Incrementally adding {len(new_samples)} new samples")

        # 1. 为新样本建立KMC数据库|Build KMC databases for new samples
        from .kmer_counter import KMCCounter
        counter = KMCCounter(self.config, self.logger_manager)

        # 为每个新样本建立数据库|Build database for each new sample
        # 跳过已经存在的数据库|Skip already existing databases
        results = {}
        for sample_input, sample_name in zip(new_samples, new_sample_names):
            # 检查KMC数据库是否已存在|Check if KMC database already exists
            # KMC数据库文件格式: {sample_name}.kmc_pre 和 {sample_name}.kmc_suf
            kmc_pre = self.config.kmc_db_path / f'{sample_name}.kmc_pre'
            kmc_suf = self.config.kmc_db_path / f'{sample_name}.kmc_suf'

            if kmc_pre.exists() and kmc_suf.exists():
                self.logger.info(f"样本数据库已存在，跳过统计|Sample database already exists, skipping: {sample_name}")
                results[sample_name] = True
            else:
                self.logger.info(f"统计新样本|Counting new sample: {sample_name}")
                success = counter.count_sample(sample_input, sample_name)
                results[sample_name] = success

        if not all(results.values()):
            self.logger.error("部分样本统计失败|Some samples counting failed")
            return False

        # 2. 加载现有 k-mer 信息|Load existing k-mer info
        import json
        metadata_file = self.config.output_path / 'kmer_metadata.json'

        if not metadata_file.exists():
            self.logger.error("元数据文件不存在，无法增量添加|Metadata file not exists, cannot incremental add")
            return False

        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        old_n_kmers = metadata['n_kmers']
        self.logger.info(f"现有 k-mer 数量|Existing k-mers: {old_n_kmers:,}")

        # 3. 增量式 Union：global_kmers + 新样本|Incremental union: global_kmers + new samples
        self.logger.info("步骤1/6: 增量式 Union 全局数据库|Step 1/6: Incremental union global database")

        # 获取现有样本名称|Get existing sample names
        with h5py.File(self.matrix_file, 'r') as f:
            existing_samples = [s.decode('utf-8') for s in f['sample_names'][:]]

        self.logger.info(f"现有样本|Existing samples: {existing_samples}")
        self.logger.info(f"新样本|New samples: {new_sample_names}")

        # 增量式 Union：global_kmers + 第一个新样本|Incremental union: global_kmers + first new sample
        global_db = self.config.get_global_db_path()
        global_db_path = Path(global_db)  # 转换为Path对象|Convert to Path object

        # 备份旧的全局数据库|Backup old global database
        # KMC数据库由两个文件组成：.kmc_pre 和 .kmc_suf
        backup_global_db = self.config.output_path / 'global_kmers_backup'
        backup_global_db.mkdir(parents=True, exist_ok=True)

        # 复制KMC数据库文件|Copy KMC database files
        for suffix in ['.kmc_pre', '.kmc_suf']:
            src = self.config.output_path / f'global_kmers{suffix}'
            dst = backup_global_db / f'global_kmers{suffix}'
            if src.exists():
                shutil.copy2(src, dst)

        try:
            # 逐个添加新样本到全局数据库|Add new samples one by one to global database
            current_global = 'global_kmers'
            for i, new_sample in enumerate(new_sample_names):
                self.logger.info(f"  合并新样本|Merging new sample: {new_sample}")

                temp_output = f'global_kmers_temp_{i}'
                if not counter.union_two_databases(current_global, new_sample, temp_output):
                    self.logger.error(f"合并样本失败|Union sample failed: {new_sample}")
                    # 恢复备份|Restore backup
                    if backup_global_db.exists():
                        for suffix in ['.kmc_pre', '.kmc_suf']:
                            src = backup_global_db / f'global_kmers{suffix}'
                            dst = self.config.output_path / f'global_kmers{suffix}'
                            if src.exists():
                                shutil.copy2(src, dst)
                    return False

                # 更新current_global为刚生成的数据库|Update current_global to the newly generated database
                current_global = temp_output

            # 将最终的临时数据库重命名为global_kmers|Rename final temp database to global_kmers
            if current_global != 'global_kmers':
                for suffix in ['.kmc_pre', '.kmc_suf']:
                    src = self.config.output_path / f'{current_global}{suffix}'
                    dst = self.config.output_path / f'global_kmers{suffix}'
                    if src.exists():
                        if dst.exists():
                            os.remove(dst)
                        shutil.move(str(src), str(dst))

            # 导出更新后的全局 k-mer 列表|Export updated global k-mer list
            kmer_list_file = self.config.output_path / 'global_kmers.txt'
            cmd = f"{self.config.get_kmc_tools_bin()} transform {global_db} dump {kmer_list_file}"
            result = os.system(cmd)
            if result != 0:
                self.logger.error("导出k-mer列表失败|Export k-mer list failed")
                return False

            # 统计新增的 k-mer|Count new k-mers
            new_n_kmers = 0
            with open(kmer_list_file, 'r') as f:
                for _ in f:
                    new_n_kmers += 1

            added_kmers = new_n_kmers - old_n_kmers
            self.logger.info(f"新增 k-mer 数量|New k-mers added: {added_kmers:,}")
            self.logger.info(f"更新后 k-mer 总数|Total k-mers after update: {new_n_kmers:,}")

        except Exception as e:
            self.logger.error(f"增量式 Union 失败|Incremental union failed: {e}")
            # 恢复备份|Restore backup
            if backup_global_db.exists():
                for suffix in ['.kmc_pre', '.kmc_suf']:
                    src = backup_global_db / f'global_kmers{suffix}'
                    dst = self.config.output_path / f'global_kmers{suffix}'
                    if src.exists():
                        shutil.copy2(src, dst)
            return False

        # 4. 导出新增的 k-mer（只导出新增的）|Export only new k-mers
        self.logger.info("步骤2/6: 导出新增 k-mer|Step 2/6: Export new k-mers")
        new_kmers_file = self._export_new_kmers_only(kmer_list_file, old_n_kmers, added_kmers)

        # 5. 更新 k-mer 字典（追加而不是重建）|Update k-mer dictionary (append, not rebuild)
        self.logger.info("步骤3/6: 更新 k-mer 字典|Step 3/6: Update k-mer dictionary")
        self._append_kmers_to_dictionary(new_kmers_file, new_n_kmers)

        # 6. 更新 RocksDB 索引（追加而不是重建）|Update RocksDB index (append, not rebuild)
        self.logger.info("步骤4/6: 更新 RocksDB 索引|Step 4/6: Update RocksDB index")
        self._append_to_rocksdb(new_kmers_file, old_n_kmers)

        # 7. 查询新样本的丰度|Query new sample abundances
        self.logger.info("步骤5/6: 查询新样本丰度|Step 5/6: Query new sample abundances")

        # 加载或创建索引（使用完整的全局 k-mer 列表，而不是只有新增的）|Load or create index (use complete global k-mer list, not only new ones)
        kmer_file_for_index = kmer_list_file  # 使用完整的 k-mer 列表|Use complete k-mer list
        kmer_file_size_gb = os.path.getsize(kmer_file_for_index) / (1024**3)
        kmer_index, index_type = self._load_or_create_index(kmer_file_for_index, kmer_file_size_gb)

        for i, sample_name in enumerate(new_sample_names, len(existing_samples)):
            self.logger.info(f"处理新样本|Processing new sample: {sample_name}")

            abundances = self._query_sample_abundances(sample_name, kmer_index, index_type)
            self._save_sample_abundances(sample_name, i, abundances)

        # 关闭 RocksDB（如果是数据库索引）|Close RocksDB if database index
        if index_type == 'rocksdb':
            kmer_index.close()
            self.logger.info("RocksDB索引已关闭|RocksDB index closed")

        # 8. 如果有新增 k-mer，查询旧样本中的这些新 k-mer|If new k-mers, query them in old samples
        if added_kmers > 0:
            self.logger.info("步骤6/6: 查询旧样本中的新 k-mer|Step 6/6: Query new k-mers in old samples")
            self._query_old_samples_for_new_kmers_incremental(existing_samples, new_kmers_file, old_n_kmers)

        # 9. 更新元数据|Update metadata
        self.logger.info("更新元数据|Updating metadata")
        with open(metadata_file, 'w') as f:
            json.dump({
                'n_kmers': new_n_kmers,
                'kmer_file': str(kmer_list_file),
                'index_type': metadata.get('index_type'),
                'index_path': metadata.get('index_path')
            }, f, indent=2)

        # 10. 重新生成 global_kmers_with_id.txt（用于后续）|Regenerate global_kmers_with_id.txt
        self.logger.info("重新生成带 ID 的 k-mer 文件|Regenerating k-mer file with IDs")
        self._add_kmer_id_column(kmer_list_file, self.config.output_path / 'global_kmers_with_id.txt')

        # 11. 更新 HDF5 文件的样本列表和元数据|Update HDF5 file sample list and metadata
        all_samples = existing_samples + new_sample_names
        with h5py.File(self.matrix_file, 'a') as f:
            # 更新样本名称|Update sample names
            if 'sample_names' in f:
                del f['sample_names']
            f.create_dataset('sample_names', data=np.array(all_samples, dtype='S'))

            # 更新样本数量|Update sample count
            f.attrs['n_samples'] = len(all_samples)

            self.logger.info(f"HDF5 文件已更新：样本数量 {len(all_samples)}|"
                            f"HDF5 file updated: {len(all_samples)} samples")

        self.logger.info("增量式添加完成|Incremental add completed")

        return True

    def _extend_matrix(self, new_n_kmers: int, all_samples: List[str]):
        """扩展矩阵|Extend matrix"""
        self.logger.info(f"扩展矩阵|Extending matrix: {new_n_kmers} k-mers")

        # 备份旧矩阵|Backup old matrix
        backup_file = self.matrix_file.with_suffix('.h5.bak')
        os.rename(self.matrix_file, backup_file)

        try:
            with h5py.File(backup_file, 'r') as f_old:
                storage_type = f_old.attrs['storage_type']
                # 兼容字符串和bytes|Compatible with both str and bytes
                if isinstance(storage_type, bytes):
                    storage_type_str = storage_type.decode('utf-8')
                else:
                    storage_type_str = storage_type

                with h5py.File(self.matrix_file, 'w') as f_new:
                    # 复制旧数据|Copy old data
                    if storage_type_str == 'sparse':
                        # 稀疏存储|Sparse storage
                        old_size = f_old['kmer_id'].shape[0]
                        f_new.create_dataset('kmer_id', data=f_old['kmer_id'][:],
                                            maxshape=(None,), compression='gzip')
                        f_new.create_dataset('sample_id', data=f_old['sample_id'][:],
                                            maxshape=(None,), compression='gzip')
                        f_new.create_dataset('abundance', data=f_old['abundance'][:],
                                            maxshape=(None,), compression='gzip')
                    else:
                        # 密集矩阵|Dense matrix
                        old_shape = f_old['abundance'].shape
                        new_shape = (new_n_kmers, old_shape[1])

                        f_new.create_dataset('abundance', shape=new_shape, dtype=np.uint16,
                                            compression='gzip', chunks=True)
                        f_new['abundance'][:old_shape[0], :] = f_old['abundance'][:]

                    # 更新元数据|Update metadata
                    for attr_name in f_old.attrs:
                        if attr_name != 'n_kmers':
                            f_new.attrs[attr_name] = f_old.attrs[attr_name]

                    f_new.attrs['n_kmers'] = new_n_kmers
                    f_new.create_dataset('sample_names', data=np.array(all_samples, dtype='S'))

        except Exception as e:
            self.logger.error(f"扩展矩阵失败|Extend matrix failed: {e}")
            # 恢复备份|Restore backup
            os.rename(backup_file, self.matrix_file)
            return False

        # 删除备份|Delete backup
        os.remove(backup_file)

        return True

    def _query_old_samples_for_new_kmers(self, old_samples: List[str],
                                         old_kmer_dict: Dict[str, int],
                                         new_kmer_dict: Dict[str, int]):
        """查询旧样本中新增k-mer的丰度|Query new k-mer abundances in old samples"""
        # 获取新增的k-mer|Get new k-mers
        new_kmers = set(new_kmer_dict.keys()) - set(old_kmer_dict.keys())

        if not new_kmers:
            return

        self.logger.info(f"查询旧样本中 {len(new_kmers)} 个新k-mer|"
                        f"Querying {len(new_kmers)} new k-mers in old samples")

        # 构建新k-mer的子字典|Build sub-dictionary of new k-mers
        new_kmer_subdict = {k: new_kmer_dict[k] for k in new_kmers}

        for sample_name in old_samples:
            self.logger.info(f"查询样本|Querying sample: {sample_name}")

            abundances = self._query_sample_abundances(sample_name, new_kmer_subdict)
            sample_id = old_samples.index(sample_name)

            self._save_sample_abundances(sample_name, sample_id, abundances)

    def get_matrix_statistics(self) -> Dict:
        """获取矩阵统计信息|Get matrix statistics"""
        if not os.path.exists(self.matrix_file):
            return {}

        stats = {}

        with h5py.File(self.matrix_file, 'r') as f:
            stats['n_kmers'] = int(f.attrs['n_kmers'])
            stats['n_samples'] = int(f.attrs['n_samples'])
            stats['kmer_size'] = int(f.attrs['kmer_size'])
            # storage_type 可能是bytes或str，需要处理|storage_type may be bytes or str
            storage_type = f.attrs['storage_type']
            if isinstance(storage_type, bytes):
                stats['storage_type'] = storage_type.decode('utf-8')
            else:
                stats['storage_type'] = storage_type
            stats['sample_names'] = [s.decode('utf-8') if isinstance(s, bytes) else s for s in f['sample_names'][:]]

            if stats['storage_type'] == 'sparse':
                # 过滤掉丰度为0的无效记录（预分配但未使用的空间）
                # Filter out invalid zero-abundance records (pre-allocated but unused space)
                abundances = f['abundance'][:]
                valid_entries = np.sum(abundances > 0)
                stats['nnz_entries'] = valid_entries
                stats['total_entries'] = f['kmer_id'].shape[0]  # 保留总条目数用于调试

                # 防止除以0|Prevent division by zero
                total_cells = stats['n_kmers'] * stats['n_samples']
                if total_cells > 0:
                    sparsity = valid_entries / total_cells * 100
                    stats['sparsity_percent'] = f"{sparsity:.2f}%"
                else:
                    stats['sparsity_percent'] = "N/A (n_kmers or n_samples is 0)"
                    self.logger.warning(f"警告：n_kmers或n_samples为0，无法计算稀疏度|"
                                       f"Warning: n_kmers or n_samples is 0, cannot calculate sparsity")

        return stats

    def _export_new_kmers_only(self, kmer_list_file: Path, old_n_kmers: int, added_kmers: int) -> Path:
        """只导出新增的 k-mer|Export only new k-mers

        Args:
            kmer_list_file: 完整的 k-mer 列表文件|Complete k-mer list file
            old_n_kmers: 旧的 k-mer 数量|Old k-mer count
            added_kmers: 新增的 k-mer 数量|Added k-mer count

        Returns:
            新 k-mer 文件路径|New k-mer file path
        """
        new_kmers_file = self.config.output_path / 'new_kmers.txt'

        self.logger.info(f"导出新增的 {added_kmers} 个 k-mer|Exporting {added_kmers} new k-mers")

        with open(kmer_list_file, 'r') as f_in, open(new_kmers_file, 'w') as f_out:
            # 跳过前 old_n_kmers 个 k-mer（只写入新增的）|Skip first old_n_kmers k-mers
            for i, line in enumerate(f_in):
                if i >= old_n_kmers:
                    f_out.write(line)

        self.logger.info(f"新 k-mer 文件已保存|New k-mer file saved: {new_kmers_file}")

        return new_kmers_file

    def _append_kmers_to_dictionary(self, new_kmers_file: Path, total_n_kmers: int):
        """追加 k-mer 到字典（不重建）|Append k-mers to dictionary (do not rebuild)

        Args:
            new_kmers_file: 新 k-mer 文件|New k-mer file
            total_n_kmers: k-mer 总数|Total k-mer count
        """
        self.logger.info("追加 k-mer 到字典|Appending k-mers to dictionary")

        # 检查并修复kmer dataset的maxshape|Check and fix kmer dataset maxshape
        old_kmers_backup = None
        with h5py.File(self.kmer_dict_file, 'a') as f_check:
            if 'kmer' in f_check:
                maxshape = f_check['kmer'].maxshape
                if maxshape[0] is not None and maxshape[0] < total_n_kmers:
                    self.logger.warning(f"kmer dataset的maxshape不足（{maxshape[0]} < {total_n_kmers}），需要扩展|kmer dataset maxshape too small, need to expand")
                    self.logger.info(f"备份旧k-mer数据|Backing up old k-mer data: {f_check.attrs['n_kmers']:,}")
                    # 备份旧数据|Backup old data
                    old_kmers_backup = f_check['kmer'][:]
                    del f_check['kmer']

        # 读取现有 k-mer 数量|Read existing k-mer count
        with h5py.File(self.kmer_dict_file, 'a') as f:
            existing_n_kmers = f.attrs['n_kmers']
            self.logger.info(f"现有 k-mer 数量|Existing k-mers: {existing_n_kmers:,}")

            # 如果kmer dataset被删除了，重建它|If kmer dataset was deleted, rebuild it
            if 'kmer' not in f:
                self.logger.info(f"创建新的kmer dataset（包含旧数据+新数据）|Creating new kmer dataset (old + new data)")

                # 先写入旧的k-mer数据|First write old k-mer data
                if old_kmers_backup is not None:
                    self.logger.info(f"恢复 {len(old_kmers_backup):,} 个旧k-mer|Restoring {len(old_kmers_backup):,} old k-mers")
                    f.create_dataset('kmer', data=old_kmers_backup, compression='gzip', maxshape=(None,))
                else:
                    # 如果没有备份数据，从文件重建（这种情况不应该发生）|If no backup, rebuild from file (should not happen)
                    self.logger.warning(f"没有备份数据，从global_kmers.txt重建|No backup data, rebuilding from global_kmers.txt")
                    kmer_list_file = self.config.output_path / 'global_kmers.txt'
                    if kmer_list_file.exists():
                        kmer_list = []
                        with open(kmer_list_file, 'r') as f_in:
                            for i, line in enumerate(f_in):
                                if i < existing_n_kmers:
                                    kmer = line.strip().split()[0]
                                    kmer_list.append(kmer)
                                else:
                                    break
                        kmer_array = np.array(kmer_list, dtype='S')
                        f.create_dataset('kmer', data=kmer_array, compression='gzip', maxshape=(None,))

            # 读取新 k-mer|Read new k-mers
            new_kmers = []
            with open(new_kmers_file, 'r') as f_in:
                for line in f_in:
                    kmer = line.strip().split()[0]
                    new_kmers.append(kmer)

            # 追加到数据集|Append to dataset
            if 'kmer' in f:
                # 调整数据集大小|Resize dataset
                f['kmer'].resize((total_n_kmers,))

                # 追加新 k-mer|Append new k-mers
                kmer_array = np.array(new_kmers, dtype='S')
                f['kmer'][existing_n_kmers:total_n_kmers] = kmer_array

                # 更新属性|Update attributes
                f.attrs['n_kmers'] = total_n_kmers

                self.logger.info(f"已追加 {len(new_kmers)} 个 k-mer|Appended {len(new_kmers)} k-mers")
            else:
                self.logger.error("k-mer 数据集不存在，无法追加|k-mer dataset not exists, cannot append")

    def _append_to_rocksdb(self, new_kmers_file: Path, old_n_kmers: int):
        """追加到 RocksDB（不重建）|Append to RocksDB (do not rebuild)

        Args:
            new_kmers_file: 新 k-mer 文件|New k-mer file
            old_n_kmers: 旧的 k-mer 数量|Old k-mer count
        """
        self.logger.info("追加到 RocksDB|Appending to RocksDB")

        # 打开现有 RocksDB|Open existing RocksDB
        db_path = self.config.output_path / 'kmer_index.rdb'
        db = DB(str(db_path), create_if_missing=False,
               max_memory_gb=self.config.max_memory)  # 使用配置的最大内存|Use configured max memory

        # 使用大批次以提高写入速度|Use large batch size for better write performance
        batch_size = 100000000  # 每1亿条提交一次|Commit every 100M entries
        batch = db.write_batch()

        count = 0
        batch_count = 0

        with open(new_kmers_file, 'r') as f:
            for line in f:
                kmer = line.strip().split()[0]
                kmer_id = old_n_kmers + count
                batch.put(kmer.encode(), str(kmer_id).encode())
                count += 1
                batch_count += 1

                # 每1亿条提交一次|Commit every 100M entries
                if batch_count >= batch_size:
                    db.write(batch)  # 调用RocksDBAdapter.write()|Call RocksDBAdapter.write()
                    self.logger.info(f"  已写入|Written: {count:,}")
                    batch = db.write_batch()
                    batch_count = 0

        # 写入剩余记录|Write remaining records
        if batch_count > 0:
            db.write(batch)  # 调用RocksDBAdapter.write()|Call RocksDBAdapter.write()

        db.close()

        self.logger.info(f"RocksDB 追加完成|RocksDB append completed: {count:,} entries")

    def _query_old_samples_for_new_kmers_incremental(self, old_samples: List[str],
                                                       new_kmers_file: Path, old_n_kmers: int):
        """增量式查询旧样本中的新 k-mer|Incrementally query new k-mers in old samples

        只查询新 k-mer，不查询旧 k-mer
        Only query new k-mers, not old k-mers

        Args:
            old_samples: 旧样本列表|Old sample list
            new_kmers_file: 新 k-mer 文件|New k-mer file
            old_n_kmers: 旧的 k-mer 数量（用于计算偏移）|Old k-mer count (for offset)
        """
        self.logger.info(f"查询 {len(old_samples)} 个旧样本中的新 k-mer|"
                        f"Querying new k-mers in {len(old_samples)} old samples")

        # 构建新 k-mer 字典（带偏移的 ID）|Build new k-mer dict (with offset IDs)
        new_kmer_dict = {}
        with open(new_kmers_file, 'r') as f:
            for i, line in enumerate(f):
                kmer = line.strip().split()[0]
                kmer_id = old_n_kmers + i  # 添加偏移|Add offset
                new_kmer_dict[kmer] = kmer_id

        self.logger.info(f"新 k-mer 字典大小|New k-mer dict size: {len(new_kmer_dict)}")

        # 加载或创建索引|Load or create index
        kmer_file_size_gb = os.path.getsize(new_kmers_file) / (1024**3)
        kmer_index, index_type = self._load_or_create_index(new_kmers_file, kmer_file_size_gb)

        # 查询每个旧样本|Query each old sample
        for sample_name in old_samples:
            self.logger.info(f"查询旧样本|Querying old sample: {sample_name}")

            abundances = self._query_sample_abundances(sample_name, kmer_index, index_type)

            # 只保存新 k-mer 的丰度|Only save new k-mer abundances
            sample_id = old_samples.index(sample_name)
            self._save_sample_abundances(sample_name, sample_id, abundances)

        # 关闭 RocksDB（如果是数据库索引）|Close RocksDB if database index
        if index_type == 'rocksdb':
            kmer_index.close()
            self.logger.info("RocksDB索引已关闭|RocksDB index closed")

    def _add_kmer_id_column(self, input_file: Path, output_file: Path):
        """给 k-mer 文件添加 ID 列|Add ID column to k-mer file

        Args:
            input_file: 输入文件（只有 k-mer 和 count）|Input file (k-mer and count only)
            output_file: 输出文件（添加 ID 列）|Output file (with ID column)
        """
        self.logger.info(f"添加 ID 列|Adding ID column: {input_file} → {output_file}")

        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            for i, line in enumerate(f_in):
                parts = line.strip().split()
                kmer = parts[0]
                count = parts[1] if len(parts) > 1 else '0'
                f_out.write(f"{kmer}\t{i}\t{count}\n")

        self.logger.info(f"ID 列已添加|ID column added: {output_file}")


# 模块级函数用于多线程查询|Module-level function for multithreading
def _query_thread_worker(args):
    """多线程查询工作函数|Multithreading query worker function

    Args:
        args: (chunk, db) 元组，chunk是(kmer, count)列表，db是共享的RocksDB对象|(chunk, db) tuple where chunk is (kmer, count) list and db is shared RocksDB object

    Returns:
        {kmer_id: count} 字典|{kmer_id: count} dict
    """
    chunk, db = args

    result = {}
    for kmer, count in chunk:
        kid_bytes = db.get(kmer.encode())
        if kid_bytes:
            kmer_id = int(kid_bytes)
            result[kmer_id] = count

    return result

