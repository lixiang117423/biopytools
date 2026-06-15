"""
RocksDB导入模块|RocksDB Import Module
"""

import os
import gzip
import shutil
import traceback


class RocksDBImporter:
    """RocksDB导入器|RocksDB Importer"""

    def __init__(self, matrix_file: str, db_path: str,
                 input_delimiter: str = '\t',
                 batch_size: int = 20000,
                 bloom_bits: int = 15,
                 force_overwrite: bool = True,
                 header_file: str = None,
                 header_db_key: str = "kmer_header",
                 logger=None):
        """
        初始化导入器|Initialize importer

        Args:
            matrix_file: 输入矩阵文件（可以是gzip压缩）|Input matrix file (can be gzipped)
            db_path: RocksDB输出路径|RocksDB output path
            input_delimiter: 输入文件分隔符|Input file delimiter
            batch_size: 批量写入大小|Batch write size
            bloom_bits: Bloom filter位数|Bloom filter bits per key
            force_overwrite: 强制覆盖|Force overwrite
            header_file: Header文件路径|Header file path
            header_db_key: 数据库中的header key|Header key in database
            logger: 日志器|Logger
        """
        self.matrix_file = matrix_file
        self.db_path = db_path
        self.input_delimiter = input_delimiter
        self.batch_size = batch_size
        self.bloom_bits = bloom_bits
        self.force_overwrite = force_overwrite
        self.header_file = header_file
        self.header_db_key = header_db_key
        self.logger = logger

    def import_to_rocksdb(self) -> bool:
        """
        执行导入|Execute import

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

        if self.logger:
            self.logger.info(f"开始导入到RocksDB|Starting import to RocksDB: {self.db_path}")

        # 检查并删除已存在的数据库|Check and remove existing database
        if os.path.exists(self.db_path):
            if self.force_overwrite:
                if self.logger:
                    self.logger.warning(f"删除已存在的数据库|Removing existing database: {self.db_path}")
                shutil.rmtree(self.db_path)
            else:
                error_msg = f"数据库已存在|Database already exists: {self.db_path}. 使用force_overwrite=True来覆盖|Use force_overwrite=True to overwrite"
                if self.logger:
                    self.logger.error(error_msg)
                else:
                    print(error_msg)
                return False

        # 获取RocksDB选项|Get RocksDB options
        opts = self._get_rocksdb_options()

        # 尝试打开数据库|Try to open database
        db = None
        compression_preferences = [
            ("ZSTD", rocksdb.CompressionType.zstd_compression),
            ("Snappy", rocksdb.CompressionType.snappy_compression),
            ("LZ4", rocksdb.CompressionType.lz4_compression),
            ("None", rocksdb.CompressionType.no_compression)
        ]

        for comp_name, comp_type in compression_preferences:
            try:
                if self.logger:
                    self.logger.debug(f"尝试使用|Trying {comp_name} 压缩|compression...")
                opts.compression = comp_type
                db = rocksdb.DB(self.db_path, opts)
                if self.logger:
                    self.logger.info(f"使用|Using {comp_name} 压缩打开数据库|compression opened database")
                break
            except (rocksdb.errors.InvalidArgument, AttributeError) as e:
                if "not linked with the binary" in str(e) or "has no attribute" in str(e):
                    if os.path.exists(self.db_path):
                        shutil.rmtree(self.db_path)
                    continue
                else:
                    if self.logger:
                        self.logger.error(f"打开数据库时出错|Error opening database: {e}")
                    return False

        if not db:
            error_msg = "无法使用任何压缩类型打开RocksDB|Failed to open RocksDB with any compression type"
            if self.logger:
                self.logger.error(error_msg)
            else:
                print(error_msg)
            return False

        try:
            # 写入header|Write header
            batch = rocksdb.WriteBatch()
            count_in_batch = 0
            imported_records = 0

            if self.header_file and os.path.exists(self.header_file):
                header_fields = []
                with open(self.header_file, 'r') as hf:
                    for line in hf:
                        field = line.strip()
                        if field:
                            header_fields.append(field)

                if header_fields:
                    header_value = '\t'.join(header_fields)
                    batch.put(
                        self.header_db_key.encode('utf-8'),
                        header_value.encode('utf-8')
                    )
                    count_in_batch += 1
                    if self.logger:
                        self.logger.info(f"添加header|Added header: {len(header_fields)} 个样本|samples")

            # 读取并导入矩阵数据|Read and import matrix data
            open_func = gzip.open if self.matrix_file.endswith('.gz') else open

            with open_func(self.matrix_file, 'rt', encoding='utf-8') as infile:
                for line_idx, line in enumerate(infile):
                    line = line.rstrip('\n')
                    if not line:
                        continue

                    parts = line.split(self.input_delimiter, 1)
                    key_str = parts[0]

                    if not key_str:
                        continue

                    value_str = parts[1] if len(parts) > 1 else ""

                    key_bytes = key_str.encode('utf-8')
                    value_bytes = value_str.encode('utf-8')

                    batch.put(key_bytes, value_bytes)
                    count_in_batch += 1

                    # 批量写入|Batch write
                    if count_in_batch >= self.batch_size:
                        db.write(batch)
                        imported_records += batch.count()
                        if self.logger and (line_idx + 1) % (self.batch_size * 100) == 0:
                            self.logger.debug(f"已处理|Processed {line_idx + 1} 行|lines, 导入|imported {imported_records} 条记录|records")
                        batch.clear()
                        count_in_batch = 0

                # 写入剩余记录|Write remaining records
                if batch.count() > 0:
                    db.write(batch)
                    imported_records += batch.count()

            if self.logger:
                self.logger.info("-" * 30)
                self.logger.info("导入完成|Import completed")
                self.logger.info(f"总记录数|Total records: {imported_records}")

            return True

        except Exception as e:
            error_msg = f"导入时出错|Error during import: {e}"
            if self.logger:
                self.logger.error(error_msg)
                traceback.print_exc()
            else:
                print(error_msg)
            return False

    def _get_rocksdb_options(self):
        """
        获取RocksDB选项|Get RocksDB options

        Returns:
            rocksdb.Options: RocksDB选项|RocksDB options
        """
        import rocksdb

        opts = rocksdb.Options()
        opts.create_if_missing = True

        opts.table_factory = rocksdb.BlockBasedTableFactory(
            filter_policy=rocksdb.BloomFilterPolicy(bits_per_key=self.bloom_bits),
            block_cache=rocksdb.LRUCache(128 * (1024**2)),
        )

        cpu_cores = os.cpu_count()
        if cpu_cores and cpu_cores > 1:
            if cpu_cores <= 2:
                opts.max_background_compactions = 1
                opts.max_background_flushes = 1
            else:
                num_total_bg_threads = max(2, cpu_cores // 2)
                if num_total_bg_threads < 4 and cpu_cores > 2:
                    num_total_bg_threads = min(4, cpu_cores)
                opts.max_background_flushes = max(1, num_total_bg_threads // 4 if num_total_bg_threads > 3 else 1)
                opts.max_background_compactions = max(1, num_total_bg_threads - opts.max_background_flushes)
        else:
            opts.max_background_compactions = 1
            opts.max_background_flushes = 1

        opts.write_buffer_size = 128 * (1024**2)
        opts.max_write_buffer_number = 4
        opts.target_file_size_base = 128 * (1024**2)
        opts.max_bytes_for_level_base = 512 * (1024**2)
        opts.level0_file_num_compaction_trigger = 10
        opts.level0_slowdown_writes_trigger = 20
        opts.level0_stop_writes_trigger = 30

        return opts
