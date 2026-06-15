"""
RocksDB查询模块|RocksDB Query Module
"""

import os
from typing import Dict, List, Tuple


class RocksDBQuerier:
    """RocksDB查询器|RocksDB Querier"""

    def __init__(self, db_path: str,
                 bloom_bits: int = 15,
                 header_db_key: str = "kmer_header",
                 logger=None):
        """
        初始化查询器|Initialize querier

        Args:
            db_path: RocksDB数据库路径|RocksDB database path
            bloom_bits: Bloom filter位数|Bloom filter bits per key
            header_db_key: 数据库中的header key|Header key in database
            logger: 日志器|Logger
        """
        self.db_path = db_path
        self.bloom_bits = bloom_bits
        self.header_db_key = header_db_key
        self.logger = logger
        self.db = None
        self.header_fields = []

    def open_db(self) -> bool:
        """
        打开数据库|Open database

        Returns:
            bool: 是否成功|Success
        """
        try:
            import rocksdb
        except ImportError:
            error_msg = "未安装python-rocksdb|python-rocksdb not installed. 请先安装|Please install: pip install python-rocksdb"
            if self.logger:
                self.logger.error(error_msg)
            else:
                print(error_msg)
            return False

        if not os.path.exists(self.db_path):
            error_msg = f"数据库不存在|Database does not exist: {self.db_path}"
            if self.logger:
                self.logger.error(error_msg)
            else:
                print(error_msg)
            return False

        try:
            # 获取选项|Get options
            opts = self._get_rocksdb_options()
            self.db = rocksdb.DB(self.db_path, opts, read_only=True)

            # 读取header|Read header
            self._read_header()

            if self.logger:
                self.logger.info(f"数据库打开成功|Database opened successfully: {self.db_path}")
                self.logger.info(f"样本数量|Sample count: {len(self.header_fields)}")

            return True

        except Exception as e:
            error_msg = f"打开数据库失败|Failed to open database: {e}"
            if self.logger:
                self.logger.error(error_msg)
            else:
                print(error_msg)
            return False

    def _read_header(self):
        """读取数据库header|Read database header"""
        try:
            header_value_bytes = self.db.get(self.header_db_key.encode('utf-8'))
            if header_value_bytes:
                header_str = header_value_bytes.decode('utf-8').strip()
                self.header_fields = header_str.split('\t')
                if self.logger:
                    self.logger.debug(f"读取header成功|Header read successfully: {len(self.header_fields)} 个字段|fields")
            else:
                if self.logger:
                    self.logger.warning(f"Header key '{self.header_db_key}' 未找到|not found in database")
        except Exception as e:
            if self.logger:
                self.logger.error(f"读取header失败|Failed to read header: {e}")

    def query_kmers(self, kmer_list: List[str],
                    include_rc: bool = True) -> Dict[str, Dict]:
        """
        查询k-mer列表|Query k-mer list

        Args:
            kmer_list: k-mer列表|K-mer list
            include_rc: 是否包含反向互补|Whether to include reverse complement

        Returns:
            dict: 查询结果|Query results
            格式|Format: {
                "kmer1": {"found": True, "abundance": "01001..."},
                "kmer2": {"found": False, "abundance": None}
            }
        """
        if not self.db:
            error_msg = "数据库未打开|Database not opened"
            if self.logger:
                self.logger.error(error_msg)
            else:
                print(error_msg)
            return {}

        results = {}

        try:
            # 准备查询的k-mer|Prepare kmers to query
            kmers_to_query = set()
            original_to_rc = {}  # 原始k-mer到反向互补的映射|Original k-mer to RC mapping

            for kmer in kmer_list:
                kmer = kmer.strip().upper()
                if not kmer:
                    continue

                # 检查是否为有效DNA序列|Check if valid DNA sequence
                if not all(c in "ATCGN" for c in kmer):
                    if self.logger:
                        self.logger.warning(f"跳过无效k-mer|Skipping invalid k-mer: {kmer}")
                    continue

                kmers_to_query.add(kmer)

                if include_rc:
                    rc_kmer = self._reverse_complement(kmer)
                    original_to_rc[kmer] = rc_kmer
                    kmers_to_query.add(rc_kmer)

            # 批量查询|Batch query
            kmers_bytes = [k.encode('utf-8') for k in kmers_to_query]
            chunk_size = 100000

            retrieved_values = {}
            for i in range(0, len(kmers_bytes), chunk_size):
                chunk = kmers_bytes[i:i + chunk_size]
                chunk_results = self.db.multi_get(chunk)
                retrieved_values.update(chunk_results)

                if self.logger:
                    self.logger.debug(f"查询进度|Query progress: {min(i + chunk_size, len(kmers_bytes))}/{len(kmers_bytes)}")

            # 处理结果|Process results
            for original_kmer in kmer_list:
                original_kmer = original_kmer.strip().upper()
                if not original_kmer:
                    continue

                # 尝试获取原始k-mer的值|Try to get value for original k-mer
                value_bytes = retrieved_values.get(original_kmer.encode('utf-8'))

                # 如果没找到，尝试反向互补|If not found, try reverse complement
                if not value_bytes and include_rc and original_kmer in original_to_rc:
                    rc_kmer = original_to_rc[original_kmer]
                    value_bytes = retrieved_values.get(rc_kmer.encode('utf-8'))

                if value_bytes:
                    value_str = value_bytes.decode('utf-8')
                    results[original_kmer] = {
                        "found": True,
                        "abundance": value_str
                    }
                else:
                    results[original_kmer] = {
                        "found": False,
                        "abundance": None
                    }

            if self.logger:
                found_count = sum(1 for r in results.values() if r["found"])
                self.logger.info(f"查询完成|Query completed: {found_count}/{len(results)} 个k-mer找到|kmers found")

            return results

        except Exception as e:
            error_msg = f"查询k-mer时出错|Error querying kmers: {e}"
            if self.logger:
                self.logger.error(error_msg)
            else:
                print(error_msg)
            return {}

    def write_results(self, results: Dict[str, Dict], output_file: str,
                     delimiter: str = '\t'):
        """
        写入查询结果到文件|Write query results to file

        Args:
            results: 查询结果|Query results
            output_file: 输出文件路径|Output file path
            delimiter: 分隔符|Delimiter
        """
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                # 写入表头|Write header
                if self.header_fields:
                    header_line = f"kmer{delimiter}{delimiter.join(self.header_fields)}\n"
                    f.write(header_line)
                else:
                    f.write(f"kmer{delimiter}abundance\n")

                # 写入数据|Write data
                for kmer, result in results.items():
                    if result["found"]:
                        abundance = result["abundance"]

                        # 如果有header，拆分丰度值|If header exists, split abundance values
                        if self.header_fields and len(abundance) == len(self.header_fields):
                            values = delimiter.join(abundance)
                            f.write(f"{kmer}{delimiter}{values}\n")
                        elif self.header_fields and len(abundance) != len(self.header_fields):
                            # 长度不匹配，作为整体输出|Length mismatch, output as whole
                            f.write(f"{kmer}{delimiter}{abundance}\n")
                        else:
                            f.write(f"{kmer}{delimiter}{abundance}\n")
                    else:
                        # 未找到，写入占位符|Not found, write placeholders
                        if self.header_fields:
                            placeholders = delimiter.join(['-'] * len(self.header_fields))
                            f.write(f"{kmer}{delimiter}{placeholders}\n")
                        else:
                            f.write(f"{kmer}{delimiter}-\n")

            if self.logger:
                self.logger.info(f"查询结果已写入|Query results written to: {output_file}")

        except Exception as e:
            error_msg = f"写入结果文件失败|Failed to write results file: {e}"
            if self.logger:
                self.logger.error(error_msg)
            else:
                print(error_msg)

    def _reverse_complement(self, sequence: str) -> str:
        """
        计算反向互补序列|Calculate reverse complement

        Args:
            sequence: DNA序列|DNA sequence

        Returns:
            str: 反向互补序列|Reverse complement sequence
        """
        complement_map = str.maketrans("ATCGN", "TAGCN")
        complement_seq = sequence.translate(complement_map)
        return complement_seq[::-1]

    def _get_rocksdb_options(self):
        """
        获取RocksDB选项|Get RocksDB options

        Returns:
            rocksdb.Options: RocksDB选项|RocksDB options
        """
        import rocksdb

        opts = rocksdb.Options()
        try:
            opts.compression = rocksdb.CompressionType.zstd_compression
        except AttributeError:
            try:
                opts.compression = rocksdb.CompressionType.snappy_compression
            except AttributeError:
                opts.compression = rocksdb.CompressionType.no_compression

        opts.table_factory = rocksdb.BlockBasedTableFactory(
            filter_policy=rocksdb.BloomFilterPolicy(bits_per_key=self.bloom_bits),
            block_cache=rocksdb.LRUCache(128 * (1024**2)),
        )

        return opts

    def close(self):
        """关闭数据库|Close database"""
        if self.db:
            self.db = None
            if self.logger:
                self.logger.info("数据库已关闭|Database closed")
