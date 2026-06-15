"""
解析器模块|Parser Module
"""

from .utils import print_interval_progress, format_number, open_input
from .parallel import ChunkProcessor, convert_chunk_to_binary


class KmerAbundanceParser:
    """Kmer丰度解析器|Kmer Abundance Parser"""

    def __init__(self, config, logger):
        """
        初始化解析器|Initialize parser

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger

    def convert_to_binary_format(self):
        """
        Pass 1: 转换丰度矩阵为0/1格式|Convert abundance matrix to binary format
        使用多线程分块处理|Use multithreaded chunk processing

        Returns:
            int: 总行数|Total lines processed
        """
        self.logger.info("=" * 100)
        self.logger.info(" Step 1: 转换丰度矩阵为0/1格式|Step 1: Converting abundance matrix to binary format")
        self.logger.info("=" * 100)
        self.logger.info(f"样本数|Number of samples: {self.config.num_samples}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")

        # 判断是否使用分桶模式|Check if use bucket mode
        use_bucket_mode = self.config.should_use_bucket_mode()

        if use_bucket_mode:
            self.logger.info("检测到超大数据集|Large dataset detected")
            self.logger.info("启用分桶处理模式|Enabling bucket processing mode")
            return self._convert_with_bucket_mode()
        else:
            self.logger.info("使用标准处理模式|Using standard processing mode")
            return self._convert_with_standard_mode()

    def _convert_with_standard_mode(self):
        """
        标准模式：多线程分块处理|Standard mode: multithreaded chunk processing

        Returns:
            int: 总行数|Total lines processed
        """
        temp_binary_file = self.config.get_temp_file_path('binary_matrix.txt')

        try:
            # 使用多线程分块处理|Use multithreaded chunk processing
            if self.config.threads > 1:
                processor = ChunkProcessor(self.config, self.logger)
                total_lines = processor.process_file_in_chunks(
                    input_file=self.config.input_matrix,
                    output_file=temp_binary_file,
                    process_func=convert_chunk_to_binary,
                    skip_header=False
                )
            else:
                # 单线程处理（fallback）|Single-threaded processing (fallback)
                total_lines = self._convert_single_thread(temp_binary_file)

            self.logger.info(f"总行数|Total lines processed: {format_number(total_lines)}")
            self.logger.info(f"临时文件|Temporary file: {temp_binary_file}")

        except Exception as e:
            self.logger.error(f"转换失败|Conversion failed: {e}")
            raise

        return total_lines

    def _convert_with_bucket_mode(self):
        """
        分桶模式：一次扫描完成转换+分桶+统计|Bucket mode: single pass conversion+bucketing+counting

        Returns:
            int: 总行数|Total lines processed
        """
        from .parallel_optimized import BucketProcessor

        bucket_dir = self.config.get_temp_file_path('buckets')

        try:
            processor = BucketProcessor(self.config, self.logger)

            # 一次扫描：转换 + 分桶|Single pass: convert + bucket
            total_lines, bucket_counts = processor.process_and_bucket(
                input_file=self.config.input_matrix,
                bucket_dir=bucket_dir,
                skip_header=True
            )

            self.logger.info(f"总行数|Total lines processed: {format_number(total_lines)}")
            self.logger.info(f"分桶目录|Bucket directory: {bucket_dir}")

        except Exception as e:
            self.logger.error(f"分桶转换失败|Bucket conversion failed: {e}")
            raise

        return total_lines

    def _convert_single_thread(self, output_file: str) -> int:
        """
        单线程转换（fallback方法）|Single-threaded conversion (fallback method)

        Args:
            output_file: 输出文件|Output file path

        Returns:
            int: 总行数|Total lines processed
        """
        total_lines = 0

        with open_input(self.config.input_matrix) as f_in, \
             open(output_file, 'w') as f_out:

            # 跳过表头并写入新表头|Skip header and write new header
            header = f_in.readline().strip()
            f_out.write(header + '\n')

            for line in f_in:
                total_lines += 1

                # 进度显示|Progress display (每100万行|Every 1M lines)
                if total_lines % 1000000 == 0:
                    print_interval_progress(total_lines, 1000000)

                line = line.rstrip('\n')

                # 尝试制表符分隔|Try tab delimiter first
                parts = line.split('\t')

                # 如果制表符分隔结果少于2列，尝试空格分隔|If tab gives <2 columns, try space
                if len(parts) < 2:
                    parts = line.split()

                if len(parts) < 2:
                    continue

                kmer_id = parts[0]
                abundances = parts[1:]

                # 转换丰度为0/1|Convert abundances to binary
                binary_values = ['1' if int(abd) > 0 else '0' for abd in abundances]

                # 写入临时文件|Write to temporary file
                f_out.write(f"{kmer_id}\t{' '.join(binary_values)}\n")

            print()  # 换行|New line

        return total_lines
