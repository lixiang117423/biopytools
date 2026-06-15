"""Kmer查询器|Kmer Query

从RocksDB数据库查询kmer
Query kmers from RocksDB database
"""

import os
from typing import Dict, Set


class KmerQuery:
    """Kmer查询器|Kmer Query

    从RocksDB数据库查询kmer及其样本存在情况
    Query kmers and their sample presence from RocksDB database
    """

    def __init__(self, config, logger):
        """初始化Kmer查询器|Initialize Kmer query

        Args:
            config: KmerPAVConfig配置对象|KmerPAVConfig object
            logger: 日志器|Logger instance
        """
        self.config = config
        self.logger = logger

    def reverse_complement(self, dna_seq: str) -> str:
        """计算反向互补序列|Calculate reverse complement sequence

        Args:
            dna_seq: DNA序列|DNA sequence

        Returns:
            str: 反向互补序列|Reverse complement sequence
        """
        complement_map = str.maketrans("ATCGN", "TAGCN")
        complement_seq = dna_seq.translate(complement_map)
        return complement_seq[::-1]

    def query_kmers_from_rocksdb(self) -> bool:
        """从RocksDB查询kmer|Query kmers from RocksDB

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info("从RocksDB查询kmer|Querying kmers from RocksDB")

        # 检查rocksdb是否可用|Check if rocksdb is available
        try:
            import rocksdb
        except ImportError:
            self.logger.error("rocksdb Python库未安装|rocksdb Python library not installed")
            self.logger.error("请安装: pip install python-rocksdb|Please install: pip install python-rocksdb")
            return False

        # 检查数据库是否存在|Check if database exists
        if not os.path.exists(self.config.rocksdb_path):
            self.logger.error(f"RocksDB数据库不存在|RocksDB database not found: {self.config.rocksdb_path}")
            return False

        # 检查kmer输入文件|Check kmer input file
        if not os.path.exists(self.config.kmer_input):
            self.logger.error(f"Kmer输入文件不存在|Kmer input file not found: {self.config.kmer_input}")
            return False

        # 打开数据库|Open database
        try:
            # 简化的RocksDB配置|Simplified RocksDB options
            opts = rocksdb.Options()
            opts.create_if_missing = False

            db = rocksdb.DB(self.config.rocksdb_path, opts, read_only=True)
            self.logger.info(f"已打开RocksDB数据库|RocksDB database opened: {self.config.rocksdb_path}")

        except Exception as e:
            self.logger.error(f"打开RocksDB失败|Failed to open RocksDB: {e}")
            return False

        # 读取kmer文件|Read kmer file
        kmers_to_query = set()
        original_kmer_to_data = {}

        try:
            with open(self.config.kmer_input, 'r') as f:
                for line in f:
                    original_kmer = line.strip().upper()
                    if not original_kmer or not all(c in "ATCGN" for c in original_kmer):
                        continue

                    rc_kmer = self.reverse_complement(original_kmer)
                    kmers_to_query.add(original_kmer)
                    kmers_to_query.add(rc_kmer)

                    original_kmer_to_data[original_kmer] = {
                        "rc_kmer": rc_kmer,
                        "value_str": None,
                        "status": "PENDING"
                    }

        except Exception as e:
            self.logger.error(f"读取kmer文件失败|Failed to read kmer file: {e}")
            return False

        self.logger.info(f"待查询kmer数量|Kmers to query: {len(kmers_to_query)}")

        # 查询kmer|Query kmers
        try:
            # 批量查询|Batch query
            chunk_size = 100000
            kmers_bytes = [k.encode('utf-8') for k in kmers_to_query]

            retrieved_values = {}
            for i in range(0, len(kmers_bytes), chunk_size):
                chunk = kmers_bytes[i:i+chunk_size]
                retrieved_values.update(db.multi_get(chunk))

            # 处理查询结果|Process query results
            found_count = 0
            not_found_count = 0

            for original_kmer, data in original_kmer_to_data.items():
                original_bytes = original_kmer.encode('utf-8')
                rc_bytes = data["rc_kmer"].encode('utf-8')

                # 尝试获取原始或反向互补kmer的值
                # Try to get value from original or reverse complement kmer
                value_bytes = retrieved_values.get(original_bytes)
                if not value_bytes:
                    value_bytes = retrieved_values.get(rc_bytes)

                if value_bytes:
                    data["value_str"] = value_bytes.decode('utf-8')
                    data["status"] = "FOUND"
                    found_count += 1
                else:
                    data["status"] = "NOT_FOUND"
                    not_found_count += 1

            self.logger.info(f"查询完成|Query completed")
            self.logger.info(f"找到kmer|Found: {found_count}")
            self.logger.info(f"未找到kmer|Not found: {not_found_count}")

            # 输出结果|Output results
            self._write_query_results(original_kmer_to_data)

            return True

        except Exception as e:
            self.logger.error(f"查询kmer失败|Failed to query kmers: {e}")
            return False

    def _write_query_results(self, kmer_data: Dict) -> bool:
        """写入查询结果|Write query results

        Args:
            kmer_data: kmer查询结果数据|Kmer query result data

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info(f"写入查询结果|Writing query results to: {self.config.query_output}")

        try:
            with open(self.config.query_output, 'w') as f:
                # 写入表头|Write header
                f.write("ID\t")
                f.write("\t".join([str(i) for i in range(100)]))  # 简化处理
                f.write("\n")

                # 写入数据|Write data
                for kmer, data in kmer_data.items():
                    if data["status"] == "FOUND":
                        f.write(f"{kmer}\t{data['value_str']}\n")
                    else:
                        f.write(f"{kmer}\t")
                        f.write("\t".join(["-"] * 100))  # 简化处理
                        f.write("\n")

            self.logger.info("查询结果已写入|Query results written")
            return True

        except Exception as e:
            self.logger.error(f"写入查询结果失败|Failed to write query results: {e}")
            return False
