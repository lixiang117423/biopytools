"""
优化版多线程处理模块|Optimized Multithreading Processing Module

适用于大规模数据处理（100亿行+）
Optimized for large-scale data processing (10B+ rows)

核心优化|Core Optimization:
- 一次扫描完成：转换 + 分桶 + 统计
- Single pass: conversion + bucketing + frequency counting
- 智能桶数量计算
- Auto-tuning bucket count based on dataset size and memory limit
"""

import os
import hashlib
from pathlib import Path
from collections import defaultdict
from typing import Callable, Tuple, Dict
from concurrent.futures import ThreadPoolExecutor, as_completed
from .utils import open_input


class BucketProcessor:
    """分桶处理器|Bucket-based Processor for Large-scale Data"""

    def __init__(self, config, logger):
        """
        初始化分桶处理器|Initialize bucket processor

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger

        # 自动计算最优参数|Auto-calculate optimal parameters
        self.num_buckets = self._calculate_optimal_buckets()
        self.bucket_memory_limit = self._calculate_bucket_memory_limit()

        self.logger.info(f"自动计算参数|Auto-calculated parameters:")
        self.logger.info(f"  分桶数量|Number of buckets: {self.num_buckets}")
        self.logger.info(f"  每桶内存限制|Memory limit per bucket: {self.bucket_memory_limit / (1024**3):.1f}GB")

    def _calculate_optimal_buckets(self) -> int:
        """
        根据数据集大小和内存限制自动计算最优桶数量|Auto-calculate optimal bucket count

        策略|Strategy:
        1. 基于输入文件大小估算
        2. 确保每个桶不超过内存限制
        3. 桶数量为2的幂次（方便哈希分配）
        4. 最少64个，最多512个

        Returns:
            int: 最优桶数量|Optimal bucket count (power of 2)
        """
        # 获取输入文件大小|Get input file size
        input_size = Path(self.config.input_matrix).stat().st_size
        size_gb = input_size / (1024**3)

        self.logger.info(f"输入文件大小|Input file size: {size_gb:.2f}GB")

        # 每个桶的目标大小|Target size per bucket
        # 根据内存限制和并发数调整|Adjust based on memory limit and concurrency
        max_memory_gb = 280  # 最大内存|Maximum memory
        threads = self.config.threads

        # 每个线程可用的内存|Memory available per thread
        memory_per_thread_gb = max_memory_gb / threads

        # 每个桶的目标大小（保守估计：每GB数据需要2GB内存排序）|Target bucket size
        # 但由于使用外部排序，可以适当放宽|But we can relax this for external sort
        target_bucket_size_gb = min(
            memory_per_thread_gb * 2,  # 每个桶不超过线程可用内存的2倍
            100  # 最大100GB per bucket
        )

        # 计算初始桶数量|Calculate initial bucket count
        initial_buckets = max(
            64,  # 最少64个桶|Minimum 64 buckets
            min(512, int(size_gb / target_bucket_size_gb))  # 最多512个
        )

        # 调整为2的幂次（方便哈希分配）|Adjust to power of 2
        num_buckets = 1
        while num_buckets < initial_buckets:
            num_buckets *= 2

        self.logger.info(f"计算桶数量|Calculated buckets:")
        self.logger.info(f"  文件大小|File size: {size_gb:.2f}GB")
        self.logger.info(f"  每桶目标大小|Target per bucket: {target_bucket_size_gb:.1f}GB")
        self.logger.info(f"  初始桶数|Initial buckets: {initial_buckets}")
        self.logger.info(f"  最终桶数（2的幂次）|Final buckets (power of 2): {num_buckets}")

        return num_buckets

    def _calculate_bucket_memory_limit(self) -> int:
        """
        计算每个桶的内存限制|Calculate memory limit per bucket

        Returns:
            int: 内存限制（字节）|Memory limit in bytes
        """
        max_memory = 280 * 1024**3  # 280GB
        threads = self.config.threads

        # 每个线程的内存限制|Memory limit per thread
        memory_per_thread = max_memory // threads

        # 但不能太大，避免单个桶排序时内存溢出|But not too large
        max_bucket_memory = 50 * 1024**3  # 最大50GB per bucket

        return min(memory_per_thread, max_bucket_memory)

    def _get_bucket_id(self, kmer_seq: str) -> int:
        """
        根据kmer序列计算分桶ID|Calculate bucket ID from kmer sequence

        使用kmer序列的hash值来分桶，确保相同kmer进入同一桶
        Use hash of kmer sequence to bucket, ensuring same kmer goes to same bucket

        Args:
            kmer_seq: kmer序列|Kmer sequence

        Returns:
            int: 桶ID (0 到 num_buckets-1)|Bucket ID
        """
        # 使用kmer的hash值来分桶|Use hash of kmer to bucket
        hash_val = hashlib.md5(kmer_seq.encode()).digest()
        bucket_id = int.from_bytes(hash_val[:4], byteorder='big') % self.num_buckets

        return bucket_id

    def process_and_bucket(
        self,
        input_file: str,
        bucket_dir: str,
        skip_header: bool = True
    ) -> Tuple[int, Dict[str, int]]:
        """
        一次扫描：转换数据 + 分桶 + 统计频次
        Single pass: convert data + bucketing + frequency counting

        Args:
            input_file: 输入文件|Input file path
            bucket_dir: 分桶目录|Bucket directory path
            skip_header: 是否跳过输入文件表头|Whether to skip input file header

        Returns:
            Tuple[int, Dict[str, int]]: (总行数, 每个桶的行数)|Total lines, lines per bucket
        """
        self.logger.info("=" * 100)
        self.logger.info(" 分桶处理模式|Bucket Processing Mode")
        self.logger.info("=" * 100)
        self.logger.info(f"分桶数量|Number of buckets: {self.num_buckets}")
        self.logger.info(f"块大小|Chunk size: {self.config.chunk_size} lines")

        # 创建分桶目录|Create bucket directory
        os.makedirs(bucket_dir, exist_ok=True)

        # 打开所有分桶文件|Open all bucket files
        bucket_files = []
        for i in range(self.num_buckets):
            bucket_path = os.path.join(bucket_dir, f'bucket_{i:04d}.txt.gz')
            # 使用快速压缩级别|Use fast compression level
            f = gzip.open(bucket_path, 'wt', compresslevel=1)
            bucket_files.append(f)

        # 统计信息|Statistics
        total_lines = 0
        bucket_counts = [0] * self.num_buckets
        header_line = None

        try:
            with open_input(input_file) as f:
                # 跳过header|Skip header
                if skip_header:
                    header_line = f.readline()
                    if header_line:
                        header_line = header_line.strip()

                # 逐行处理|Process line by line
                chunk_lines = 0
                for line in f:
                    total_lines += 1
                    chunk_lines += 1

                    # 进度显示|Progress display
                    if chunk_lines >= self.config.chunk_size:
                        from .utils import format_number
                        self.logger.info(f"已处理|Processed: {format_number(total_lines)} lines")
                        chunk_lines = 0

                    # 解析行|Parse line
                    line = line.rstrip('\n')

                    # 尝试制表符分隔|Try tab delimiter first
                    parts = line.split('\t')
                    if len(parts) < 2:
                        parts = line.split()

                    if len(parts) < 2:
                        continue

                    kmer_id = parts[0]
                    abundances = parts[1:]

                    # 数据验证|Data validation
                    try:
                        int(abundances[0])
                    except (ValueError, IndexError):
                        continue

                    # 转换为0/1|Convert to binary
                    binary_values = ['1' if int(abd) > 0 else '0' for abd in abundances]
                    processed_line = f"{kmer_id}\t{' '.join(binary_values)}"

                    # 分桶|Bucket assignment
                    bucket_id = self._get_bucket_id(kmer_id)

                    # 写入对应桶|Write to corresponding bucket
                    bucket_files[bucket_id].write(processed_line + '\n')
                    bucket_counts[bucket_id] += 1

        finally:
            # 关闭所有桶文件|Close all bucket files
            for f in bucket_files:
                f.close()

        # 输出统计|Output statistics
        from .utils import format_number
        self.logger.info(f"总行数|Total lines: {format_number(total_lines)}")
        self.logger.info(f"分桶统计|Bucket distribution:")

        # 显示每个桶的大小|Show size of each bucket
        total_in_buckets = sum(bucket_counts)
        avg_count = total_in_buckets / self.num_buckets if self.num_buckets > 0 else 0

        self.logger.info(f"平均每桶|Average per bucket: {format_number(int(avg_count))} lines")

        # 显示非空桶|Show non-empty buckets
        non_empty = sum(1 for c in bucket_counts if c > 0)
        self.logger.info(f"非空桶数|Non-empty buckets: {non_empty}/{self.num_buckets}")

        bucket_counts_dict = {f'bucket_{i:04d}.txt.gz': count for i, count in enumerate(bucket_counts)}

        return total_lines, bucket_counts_dict

    def process_buckets_in_parallel(
        self,
        bucket_dir: str,
        output_vcf: str,
        output_txt: str,
        min_agg_count: int
    ) -> int:
        """
        并行处理所有桶：统计频次 + 过滤 + 排序 + 写VCF
        Process all buckets in parallel: frequency counting + filtering + sorting + VCF writing

        Args:
            bucket_dir: 分桶目录|Bucket directory
            output_vcf: 输出VCF文件|Output VCF file
            output_txt: 输出TXT文件|Output TXT file
            min_agg_count: 最小聚合频次|Minimum aggregated count

        Returns:
            int: 保留的行数|Number of kept lines
        """
        import glob

        # 获取所有桶文件|Get all bucket files
        bucket_files = sorted(glob.glob(os.path.join(bucket_dir, 'bucket_*.txt.gz')))

        self.logger.info("=" * 100)
        self.logger.info(f" 并行处理桶文件|Processing Buckets in Parallel")
        self.logger.info("=" * 100)
        self.logger.info(f"找到|Found {len(bucket_files)} 个桶文件|bucket files")

        # 限制并发数以避免内存溢出|Limit concurrency to avoid memory overflow
        # 计算安全并发数：总内存限制(280GB) / 每桶排序内存限制
        # Calculate safe concurrent count: total memory limit (280GB) / per-bucket sort memory limit
        safe_concurrent_by_memory = int(280 * 1024**3 / self.bucket_memory_limit)

        # 最终并发数：取线程数、内存限制、桶数量中的最小值
        # Final concurrent count: min(threads, memory_limit, bucket_count)
        max_concurrent = min(
            self.config.threads,
            safe_concurrent_by_memory,
            len(bucket_files)
        )

        self.logger.info(f"内存限制计算|Memory limit calculation:")
        self.logger.info(f"  总内存限制|Total memory limit: 280GB")
        self.logger.info(f"  每桶排序内存|Per-bucket sort memory: {self.bucket_memory_limit / (1024**3):.1f}GB")
        self.logger.info(f"  内存允许并发|Memory-allowed concurrency: {safe_concurrent_by_memory}")
        self.logger.info(f"  实际并发数|Actual concurrent threads: {max_concurrent}")

        # 使用线程池并行处理|Use thread pool for parallel processing
        kept_lines = 0
        temp_vcf_parts = []
        temp_txt_parts = []

        with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            # 提交所有任务|Submit all tasks
            futures = {}
            for idx, bucket_file in enumerate(bucket_files):
                bucket_name = os.path.basename(bucket_file)
                temp_vcf = os.path.join(bucket_dir, f'{bucket_name}.vcf')
                temp_txt = os.path.join(bucket_dir, f'{bucket_name}.txt')

                future = executor.submit(
                    self._process_single_bucket,
                    bucket_file,
                    temp_vcf,
                    temp_txt,
                    min_agg_count,
                    bucket_name,
                    idx
                )
                futures[future] = (temp_vcf, temp_txt, bucket_name)

            # 收集结果|Collect results
            completed = 0
            for future in as_completed(futures):
                temp_vcf, temp_txt, bucket_name = futures[future]
                try:
                    count = future.result()
                    kept_lines += count
                    temp_vcf_parts.append(temp_vcf)
                    temp_txt_parts.append(temp_txt)
                    completed += 1

                    if completed % 10 == 0:
                        from .utils import format_number
                        self.logger.info(f"已完成|Completed: {completed}/{len(bucket_files)} 桶|buckets")

                except Exception as e:
                    self.logger.error(f"处理桶失败|Error processing bucket {bucket_name}: {e}")

        # 合并所有部分文件|Merge all partial files
        self.logger.info("=" * 100)
        self.logger.info(" 合并输出文件|Merging Output Files")
        self.logger.info("=" * 100)

        # 合并VCF|Merge VCF
        self._merge_vcf_files(temp_vcf_parts, output_vcf)

        # 合并TXT|Merge TXT
        self._merge_txt_files(temp_txt_parts, output_txt)

        # 删除临时文件|Delete temporary files
        self.logger.info("清理临时文件|Cleaning up temporary files...")
        for temp_file in temp_vcf_parts + temp_txt_parts:
            try:
                os.remove(temp_file)
            except:
                pass

        return kept_lines

    def _merge_vcf_files(self, temp_vcf_parts: list, output_vcf: str):
        """
        合并VCF文件|Merge VCF files

        Args:
            temp_vcf_parts: 临时VCF文件列表|List of temporary VCF files
            output_vcf: 输出VCF文件|Output VCF file
        """
        with open(output_vcf, 'w') as f_out:
            # 写入header（从第一个文件读取）|Write header (read from first file)
            if temp_vcf_parts:
                with open(temp_vcf_parts[0], 'r') as f_first:
                    # 写入VCF header lines（##开头）|Write VCF header lines
                    for line in f_first:
                        f_out.write(line)
                        if line.startswith('#CHROM'):
                            break

                # 追加其他文件的数据行|Append data lines from other files
                for idx, temp_vcf in enumerate(temp_vcf_parts):
                    with open(temp_vcf, 'r') as f_in:
                        # 跳过header|Skip header
                        for line in f_in:
                            if line.startswith('#'):
                                continue
                            f_out.write(line)

                    if (idx + 1) % 50 == 0:
                        self.logger.info(f"已合并|Merged: {idx + 1}/{len(temp_vcf_parts)} VCF parts")

        self.logger.info(f"VCF文件合并完成|VCF merging completed: {output_vcf}")

    def _merge_txt_files(self, temp_txt_parts: list, output_txt: str):
        """
        合并TXT文件|Merge TXT files

        Args:
            temp_txt_parts: 临时TXT文件列表|List of temporary TXT files
            output_txt: 输出TXT文件|Output TXT file
        """
        with open(output_txt, 'w') as f_out:
            # 写入header（从第一个文件读取）|Write header
            if temp_txt_parts:
                with open(temp_txt_parts[0], 'r') as f_first:
                    header = f_first.readline()
                    f_out.write(header)

                # 追加数据|Append data
                for idx, temp_txt in enumerate(temp_txt_parts):
                    with open(temp_txt, 'r') as f_in:
                        # 跳过header|Skip header
                        next(f_in)
                        f_out.write(f_in.read())

                    if (idx + 1) % 50 == 0:
                        self.logger.info(f"已合并|Merged: {idx + 1}/{len(temp_txt_parts)} TXT parts")

        self.logger.info(f"TXT文件合并完成|TXT merging completed: {output_txt}")

    def _process_single_bucket(
        self,
        bucket_file: str,
        output_vcf: str,
        output_txt: str,
        min_agg_count: int,
        bucket_name: str,
        bucket_idx: int
    ) -> int:
        """
        处理单个桶文件|Process single bucket file

        步骤|Steps:
        1. 读取桶文件并统计频次|Read bucket and count frequencies
        2. 过滤并写入临时文件|Filter and write temp file
        3. 排序|Sort
        4. 写入VCF和TXT|Write VCF and TXT

        Args:
            bucket_file: 桶文件路径|Bucket file path
            output_vcf: 输出VCF|Output VCF
            output_txt: 输出TXT|Output TXT
            min_agg_count: 最小聚合频次|Minimum aggregated count
            bucket_name: 桶名称|Bucket name
            bucket_idx: 桶索引|Bucket index

        Returns:
            int: 保留的行数|Number of kept lines
        """
        import subprocess

        try:
            # 步骤1：统计频次|Step 1: Count frequencies
            kmer_freq = defaultdict(int)

            with gzip.open(bucket_file, 'rt') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue

                    kmer_id = parts[0]
                    binary_str = ' '.join(parts[1:])
                    count = binary_str.count('1')

                    kmer_freq[kmer_id] += count

            # 步骤2：过滤并写入临时文件|Step 2: Filter and write temp file
            temp_filtered = bucket_file + '.filtered'

            kept_count = 0
            with gzip.open(bucket_file, 'rt') as f_in, \
                 open(temp_filtered, 'w') as f_out:

                for line in f_in:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue

                    kmer_id = parts[0]
                    agg_freq = kmer_freq.get(kmer_id, 0)

                    if min_agg_count <= agg_freq <= self.config.num_samples:
                        f_out.write(line)
                        kept_count += 1

            # 步骤3：排序|Step 3: Sort
            temp_sorted = temp_filtered + '.sorted'

            # 使用外部排序，控制内存使用|Use external sort with memory limit
            sort_memory_gb = self.bucket_memory_limit / (1024**3)

            subprocess.run(
                ['sort', '-T', os.path.dirname(temp_filtered),
                 '-S', f'{sort_memory_gb:.1f}G',
                 '-o', temp_sorted, temp_filtered],
                check=True,
                capture_output=True
            )

            # 步骤4：写VCF和TXT|Step 4: Write VCF and TXT
            with open(temp_sorted, 'r') as f_in, \
                 open(output_vcf, 'w') as f_vcf, \
                 open(output_txt, 'w') as f_txt:

                # 写入VCF header|Write VCF header
                f_vcf.write("##fileformat=VCFv4.2\n")
                f_vcf.write("##source=kmer2vcf\n")
                f_vcf.write(f"##filter_agg_range={min_agg_count}_to_{self.config.num_samples}\n")
                f_vcf.write("##INFO=<ID=KMER,Number=1,Type=String,Description=\"Kmer sequence\">\n")
                f_vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                f_vcf.write("##contig=<ID=Chr_kmer,length=2147483647>\n")

                samples = self.config.samples
                vcf_header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples
                f_vcf.write("\t".join(vcf_header) + "\n")

                # 写入TXT header|Write TXT header
                txt_header = ["KmerID"] + samples
                f_txt.write("\t".join(txt_header) + "\n")

                # 写入数据行|Write data rows
                pos = 1
                for line in f_in:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue

                    kmer_id = parts[0]
                    abundances = parts[1:]

                    # 写入TXT|Write TXT
                    f_txt.write(line + '\n')

                    # 写入VCF|Write VCF
                    genotypes = []
                    for abd in abundances:
                        if int(abd) > 0:
                            genotypes.append("0/0")
                        else:
                            genotypes.append("1/1")

                    vcf_row = [
                        "Chr_kmer",
                        str(pos),
                        ".",
                        kmer_id,
                        ".",
                        ".",
                        "PASS",
                        f"KMER={kmer_id}",
                        "GT"
                    ] + genotypes

                    f_vcf.write("\t".join(vcf_row) + "\n")
                    pos += 1

            # 清理临时文件|Clean up temp files
            try:
                os.remove(temp_filtered)
                os.remove(temp_sorted)
            except:
                pass

            return kept_count

        except Exception as e:
            self.logger.error(f"处理桶时出错|Error processing bucket {bucket_name}: {e}")
            raise
