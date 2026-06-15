"""
核心计算逻辑模块|Core Calculation Logic Module
"""

import gzip
import os
from typing import Dict, Set, Tuple, Optional
from multiprocessing import Pool
from .utils import format_number

# 全局变量用于worker进程共享目标kmer数据|Global variables for sharing target kmer data across workers
_global_target_kmers_dict = None
_global_target_kmers_set = None


def _parse_kmer_line(line: str) -> Tuple[Optional[str], list]:
    """
    解析kmer矩阵文件的每一行，兼容不同分隔符格式|Parse kmer matrix line, compatible with different delimiter formats

    支持的格式|Supported formats:
    1. 混合格式（制表符+空格）|Mixed format (tab+space): KMER\tval1 val2 val3 ...
    2. 纯制表符分隔|Tab-delimited: KMER\tval1\tval2\t...
    3. 纯空格分隔|Space-delimited: KMER val1 val2 ...

    Args:
        line: 数据行|Data line

    Returns:
        Tuple[Optional[str], list]: (kmer序列, 丰度列表) 或 (None, []) 如果解析失败

    Examples:
        >>> _parse_kmer_line("AAAA\t1 2 3")  # 混合格式|Mixed format
        ('AAAA', ['1', '2', '3'])
        >>> _parse_kmer_line("AAAA\t1\t2\t3")  # 纯制表符|Pure tab
        ('AAAA', ['1', '2', '3'])
        >>> _parse_kmer_line("AAAA 1 2 3")  # 纯空格|Pure space
        ('AAAA', ['1', '2', '3'])
    """
    line = line.strip()
    if not line:
        return None, []

    # 先尝试用制表符分割第一列|Try tab delimiter for first column
    parts_tab = line.split('\t', 1)  # 只分割一次|Split only once
    if len(parts_tab) == 2:
        kmer_seq = parts_tab[0].upper()
        # 检查第一列是否是纯ATCGN序列|Check if first column is pure ATCGN sequence
        if kmer_seq and all(c in 'ATCGN' for c in kmer_seq):
            # 用空格分割第二部分|Split second part with spaces
            abundances = parts_tab[1].split()
            if abundances:
                return kmer_seq, abundances

    # 如果制表符分割失败，使用空格/空白字符分割|If tab split failed, use space/whitespace split
    parts_space = line.split()
    if len(parts_space) >= 2:
        # 第一部分应该是kmer序列|First part should be kmer sequence
        kmer_seq = parts_space[0].upper()
        if kmer_seq and all(c in 'ATCGN' for c in kmer_seq):
            return kmer_seq, parts_space[1:]

    # 如果都无法解析，返回None|If parsing failed, return None
    return None, []


def _worker_init(target_kmers_dict: Dict[str, str], target_kmers_set: Set[str]):
    """
    Worker初始化函数：在fork后设置全局变量|Worker initializer: set global variables after fork

    利用Linux的copy-on-write机制，多个worker共享同一份只读数据|Using Linux copy-on-write mechanism,
    multiple workers share the same read-only data

    Args:
        target_kmers_dict: 目标kmer字典|Target kmer dictionary
        target_kmers_set: 目标kmer集合|Target kmer set
    """
    global _global_target_kmers_dict, _global_target_kmers_set
    _global_target_kmers_dict = target_kmers_dict
    _global_target_kmers_set = target_kmers_set


class KmerIntersectCalculator:
    """Kmer交集计算器|Kmer Intersection Calculator"""

    def __init__(self, config, logger):
        """
        初始化计算器|Initialize calculator

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger

    def load_target_kmers(self) -> tuple:
        """
        从fasta文件加载目标kmer|Load target kmers from fasta file

        修改：保留所有重复的kmer|Modified: Keep all duplicate kmers

        Returns:
            tuple: (kmer序列 -> kmer ID的映射（存储第一个ID）, kmer ID列表（保持fasta原始顺序，包括重复）, kmer ID -> 序列的映射)
                  (Kmer sequence to kmer ID mapping (stores first ID), Kmer ID list preserving fasta order including duplicates, kmer ID to sequence mapping)
        """
        self.logger.info("=" * 100)
        self.logger.info(" 加载目标kmer|Loading Target Kmers")
        self.logger.info("=" * 100)
        self.logger.info(f"fasta文件|Fasta file: {self.config.kmer_fasta}")

        kmer_to_id = {}  # sequence -> first ID (for query)
        kmer_id_list = []  # 保持fasta中的原始顺序，包括重复|Preserve original order from fasta including duplicates
        kmer_id_to_seq = {}  # ID -> sequence (for write results)
        seen_sequences = set()  # 用于检测重复的kmer序列|For detecting duplicate kmer sequences
        total_kmers = 0
        duplicate_count = 0

        # 判断是否为压缩文件|Check if file is compressed
        open_func = gzip.open if self.config.kmer_fasta.endswith('.gz') else open

        with open_func(self.config.kmer_fasta, 'rt') as f:
            current_id = None

            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith('>'):
                    # fasta header行|Fasta header line
                    current_id = line[1:].split()[0]  # 取第一个字段作为ID|Take first field as ID
                else:
                    # 序列行|Sequence line
                    sequence = line.upper()

                    if current_id and sequence:
                        total_kmers += 1

                        # 检查是否重复|Check for duplicates
                        if sequence in seen_sequences:
                            duplicate_count += 1
                            # 仍然添加到列表中|Still add to list
                            kmer_id_list.append(current_id)
                            kmer_id_to_seq[current_id] = sequence
                        else:
                            # 第一次出现|First occurrence
                            kmer_to_id[sequence] = current_id
                            kmer_id_list.append(current_id)  # 按fasta顺序添加ID列表|Add ID in fasta order
                            kmer_id_to_seq[current_id] = sequence
                            seen_sequences.add(sequence)

        self.logger.info(f"总kmer数|Total kmers: {format_number(total_kmers)}")
        self.logger.info(f"唯一kmer数|Unique kmers: {format_number(len(kmer_to_id))}")
        if duplicate_count > 0:
            self.logger.info(f"重复kmer数|Duplicate kmers: {format_number(duplicate_count)} (已保留|preserved)")

        return kmer_to_id, kmer_id_list, kmer_id_to_seq

    def extract_kmer_abundance(
        self,
        target_kmers: Dict[str, str],
        sample_names: list
    ) -> Tuple[Dict[str, dict], int, list]:
        """
        从kmer矩阵中提取目标kmer的丰度|Extract target kmer abundance from kmer matrix

        自动选择单进程或多进程模式|Automatically choose single or multi-process mode

        Args:
            target_kmers: 目标kmer字典 (序列 -> ID)|Target kmer dictionary (sequence -> ID)
            sample_names: 样本名列表|List of sample names

        Returns:
            Tuple[Dict[str, dict], int, list]: (kmer ID -> 丰度信息, 找到的kmer数, 样本名列表)
        """
        # 判断是否使用多进程模式|Check if use multi-process mode
        # 阈值：目标kmer数 > 1000万 或 文件大小 > 10GB
        use_parallel = len(target_kmers) > 10_000_000

        if use_parallel and self.config.threads > 1:
            # 使用多进程模式|Use multi-process mode
            return self.extract_kmer_abundance_parallel(target_kmers, sample_names)
        else:
            # 使用单进程模式|Use single-process mode
            return self.extract_kmer_abundance_single(target_kmers, sample_names)

    def extract_kmer_abundance_single(
        self,
        target_kmers: Dict[str, str],
        sample_names: list
    ) -> Tuple[Dict[str, dict], int, list]:
        """
        单进程提取kmer丰度|Extract kmer abundance in single-process mode

        Args:
            target_kmers: 目标kmer字典 (序列 -> ID)|Target kmer dictionary (sequence -> ID)
            sample_names: 样本名列表|List of sample names

        Returns:
            Tuple[Dict[str, dict], int, list]: (kmer ID -> 丰度信息, 找到的kmer数, 样本名列表)
        """
        self.logger.info("=" * 100)
        self.logger.info(" 提取kmer丰度（单进程模式）|Extracting Kmer Abundance (Single-Process Mode)")
        self.logger.info("=" * 100)
        self.logger.info(f"kmer矩阵文件|Kmer matrix file: {self.config.kmer_matrix}")
        self.logger.info(f"目标kmer数|Target kmers: {format_number(len(target_kmers))}")

        # 结果字典|Result dictionary
        # kmer_id -> {'sequence': str, 'found_as': str, 'abundances': list}
        results = {}

        # 用于快速查找的集合|Set for fast lookup
        target_set = set(target_kmers.keys())

        # 统计信息|Statistics
        found_forward = 0
        found_reverse = 0
        not_found = 0
        total_lines = 0

        # 判断是否为压缩文件|Check if file is compressed
        open_func = gzip.open if self.config.kmer_matrix.endswith('.gz') else open

        with open_func(self.config.kmer_matrix, 'rt') as f:
            # 跳过header行|Skip header line
            header_line = f.readline()

            if not sample_names:
                # 如果没有提供样本名，从header解析|If sample names not provided, parse from header
                # 尝试空格分割（优先）|Try space delimiter first (preferred)
                parts_space = header_line.strip().split()

                if len(parts_space) >= 2 and parts_space[0] in ['KMER', 'KmerID']:
                    # 格式：KMER Sample1 Sample2 Sample3 ...（纯空格分隔）|Format: KMER Sample1 Sample2... (space-only)
                    sample_names = parts_space[1:]
                else:
                    # 尝试制表符分割（兼容旧格式）|Try tab delimiter (backward compatibility)
                    parts_tab = header_line.strip().split('\t', 1)
                    if len(parts_tab) == 2:
                        if parts_tab[0] == 'KMER' or parts_tab[0] == 'KmerID':
                            # 格式：KMER<TAB>Sample1 Sample2 ...（混合格式）|Format: KMER<TAB>Sample1 Sample2... (mixed)
                            sample_names = parts_tab[1].split()
                    elif len(parts_tab) >= 4 and parts_tab[0] == 'KmerID':
                        # 兼容旧格式：KmerID Sequence FoundAs Sample1 Sample2 ...|Old format: KmerID Sequence FoundAs Sample1...
                        sample_names = parts_tab[3:]

            # 逐行处理|Process line by line
            chunk_lines = 0
            for line in f:
                total_lines += 1
                chunk_lines += 1

                # 进度显示|Progress display
                if chunk_lines >= self.config.chunk_size:
                    self.logger.info(f"已扫描|Scanned: {format_number(total_lines)} lines, "
                                   f"找到|Found: {format_number(found_forward + found_reverse)}")
                    chunk_lines = 0

                # 解析kmer行（兼容不同分隔符格式）|Parse kmer line (compatible with different delimiters)
                kmer_seq, abundances = _parse_kmer_line(line)
                if kmer_seq is None:
                    continue

                # 正向查找|Forward search
                if kmer_seq in target_set:
                    kmer_id = target_kmers[kmer_seq]
                    results[kmer_id] = {
                        'sequence': kmer_seq,
                        'found_as': 'forward',
                        'abundances': abundances
                    }
                    found_forward += 1
                    target_set.discard(kmer_seq)  # 从待查找集合中移除|Remove from target set
                    continue

                # 反向互补查找|Reverse complement search
                if self.config.use_reverse_complement:
                    from .utils import reverse_complement
                    rev_comp_seq = reverse_complement(kmer_seq)

                    if rev_comp_seq in target_set:
                        kmer_id = target_kmers[rev_comp_seq]
                        results[kmer_id] = {
                            'sequence': rev_comp_seq,  # 保留原始查询序列|Keep original query sequence
                            'found_as': 'reverse',
                            'abundances': abundances
                        }
                        found_reverse += 1
                        target_set.discard(rev_comp_seq)  # 从待查找集合中移除|Remove from target set

        not_found = len(target_set)

        self.logger.info(f"扫描完成|Scan completed")
        self.logger.info(f"总扫描行数|Total lines scanned: {format_number(total_lines)}")
        self.logger.info(f"正向找到|Found forward: {format_number(found_forward)}")
        self.logger.info(f"反向互补找到|Found reverse complement: {format_number(found_reverse)}")
        self.logger.info(f"未找到|Not found: {format_number(not_found)}")
        self.logger.info(f"总共找到|Total found: {format_number(len(results))}")

        return results, len(results), sample_names

    def _process_chunk_worker(self, args):
        """
        Worker函数：处理单个文件块（多进程）|Worker function: process single file chunk (multiprocess)

        使用全局变量共享目标kmer数据，避免重复加载|Uses global variables to share target kmer data,
        avoiding duplicate loading

        Args:
            args: (chunk_id, start_line, end_line, matrix_file, use_reverse_complement)

        Returns:
            tuple: (chunk_id, results_dict, found_count, stats)
        """
        global _global_target_kmers_dict, _global_target_kmers_set
        chunk_id, start_line, end_line, matrix_file, use_reverse_complement = args

        # 使用全局共享的目标kmer数据|Use globally shared target kmer data
        # 不创建副本，直接使用全局变量（只读操作）|Do not create copy, use global variable directly (read-only)
        target_kmers = _global_target_kmers_dict
        target_set = _global_target_kmers_set

        # 使用局部set记录已找到的kmer，避免修改全局set|Use local set to track found kmers, avoid modifying global set
        found_kmers = set()
        results = {}
        found_forward = 0
        found_reverse = 0
        total_lines = 0

        # 处理分配的行范围|Process assigned line range
        open_func_matrix = gzip.open if matrix_file.endswith('.gz') else open

        try:
            with open_func_matrix(matrix_file, 'rt') as f:
                # 跳到起始位置|Skip to start position
                for _ in range(start_line + 1):  # +1 跳过header行|Skip header line
                    next(f)

                # 处理行范围|Process line range
                for i, line in enumerate(f):
                    if i >= (end_line - start_line):
                        break

                    total_lines += 1

                    # 解析kmer行（兼容不同分隔符格式）|Parse kmer line (compatible with different delimiters)
                    kmer_seq, abundances = _parse_kmer_line(line)
                    if kmer_seq is None:
                        continue

                    # 正向查找|Forward search
                    if kmer_seq in target_set and kmer_seq not in found_kmers:
                        kmer_id = target_kmers[kmer_seq]
                        results[kmer_id] = {
                            'sequence': kmer_seq,
                            'found_as': 'forward',
                            'abundances': abundances
                        }
                        found_forward += 1
                        found_kmers.add(kmer_seq)
                        continue

                    # 反向互补查找|Reverse complement search
                    if use_reverse_complement:
                        from .utils import reverse_complement
                        rev_comp_seq = reverse_complement(kmer_seq)

                        if rev_comp_seq in target_set and rev_comp_seq not in found_kmers:
                            kmer_id = target_kmers[rev_comp_seq]
                            results[kmer_id] = {
                                'sequence': rev_comp_seq,
                                'found_as': 'reverse',
                                'abundances': abundances
                            }
                            found_reverse += 1
                            found_kmers.add(rev_comp_seq)

        except Exception as e:
            return (chunk_id, {}, 0, {'error': str(e)})

        stats = {
            'found_forward': found_forward,
            'found_reverse': found_reverse,
            'total_lines': total_lines
        }

        return (chunk_id, results, len(results), stats)

    def extract_kmer_abundance_parallel(
        self,
        target_kmers: Dict[str, str],
        sample_names: list
    ) -> Tuple[Dict[str, dict], int, list]:
        """
        多进程并行提取kmer丰度|Extract kmer abundance using multi-processing

        Args:
            target_kmers: 目标kmer字典|Target kmer dictionary
            sample_names: 样本名列表|List of sample names

        Returns:
            Tuple[Dict[str, dict], int, list]: (kmer ID -> 丰度信息, 找到的kmer数, 样本名列表)
        """
        self.logger.info("=" * 100)
        self.logger.info(" 提取kmer丰度（多进程模式）|Extracting Kmer Abundance (Multi-Process Mode)")
        self.logger.info("=" * 100)
        self.logger.info(f"kmer矩阵文件|Kmer matrix file: {self.config.kmer_matrix}")
        self.logger.info(f"目标kmer数|Target kmers: {format_number(len(target_kmers))}")
        self.logger.info(f"进程数|Number of processes: {self.config.threads}")

        # 读取header获取sample names|Read header to get sample names
        if not sample_names:
            open_func = gzip.open if self.config.kmer_matrix.endswith('.gz') else open
            with open_func(self.config.kmer_matrix, 'rt') as f:
                header_line = f.readline()
                # 尝试空格分割（优先）|Try space delimiter first (preferred)
                parts_space = header_line.strip().split()

                if len(parts_space) >= 2 and parts_space[0] in ['KMER', 'KmerID']:
                    # 格式：KMER Sample1 Sample2 Sample3 ...（纯空格分隔）|Format: KMER Sample1 Sample2... (space-only)
                    sample_names = parts_space[1:]
                else:
                    # 尝试制表符分割（兼容旧格式）|Try tab delimiter (backward compatibility)
                    parts_tab = header_line.strip().split('\t', 1)
                    if len(parts_tab) == 2:
                        if parts_tab[0] == 'KMER' or parts_tab[0] == 'KmerID':
                            # 格式：KMER<TAB>Sample1 Sample2 ...（混合格式）|Format: KMER<TAB>Sample1 Sample2... (mixed)
                            sample_names = parts_tab[1].split()
                    elif len(parts_tab) >= 4 and parts_tab[0] == 'KmerID':
                        # 兼容旧格式：KmerID Sequence FoundAs Sample1 Sample2 ...|Old format: KmerID Sequence FoundAs Sample1...
                        sample_names = parts_tab[3:]

        # 步骤1：估算文件总行数|Step 1: Estimate total lines in file
        self.logger.info("估算文件行数|Estimating file lines...")
        total_lines = self._count_lines(self.config.kmer_matrix) - 1  # 减去header行|Subtract header line
        self.logger.info(f"文件总行数|Total lines in file: {format_number(total_lines)}")

        # 步骤2：计算分块（智能调整进程数以避免内存溢出）|Step 2: Calculate chunks (smart process count adjustment to avoid OOM)
        # 估算每个worker的内存开销（考虑found_kmers局部set）|Estimate memory per worker (considering local found_kmers set)
        # 假设每个找到的kmer占用约200字节（set overhead）|Assume ~200 bytes per found kmer (set overhead)
        # 预期找到率约80%，所以每个worker最多需要存储约0.8*chunk_size个kmer|Expected find rate ~80%, so each worker stores ~0.8*chunk_size kmers
        # 安全限制：每个worker的found_kmers set不应超过500M条目（约100GB内存）|Safety limit: found_kmers set per worker should not exceed 500M entries (~100GB)
        max_safe_processes = self.config.threads

        # 基于目标kmer数量自动调整进程数|Auto-adjust process count based on target kmer count
        kmer_count = len(target_kmers)
        if kmer_count > 50_000_000:  # 超过50M kmer
            # 使用保守的进程数以避免内存溢出|Use conservative process count to avoid OOM
            recommended_processes = min(8, self.config.threads)
            if recommended_processes < self.config.threads:
                self.logger.warning(f"目标kmer数过大|Large target kmer count: {format_number(kmer_count)}")
                self.logger.warning(f"为避免内存溢出，将进程数从|To avoid OOM, reducing processes from {self.config.threads} 到|to {recommended_processes}")
            max_safe_processes = recommended_processes

        num_processes = min(max_safe_processes, total_lines)
        chunk_size = total_lines // num_processes

        self.logger.info(f"分块信息|Chunking info:")
        self.logger.info(f"  进程数|Processes: {num_processes}")
        self.logger.info(f"  每进程行数|Lines per process: {format_number(chunk_size)}")

        # 步骤3：准备任务参数（不再传递fasta文件，使用共享数据）|Step 3: Prepare task arguments (no fasta file, using shared data)
        tasks = []
        for i in range(num_processes):
            start_line = i * chunk_size
            end_line = (i + 1) * chunk_size if i < num_processes - 1 else total_lines

            tasks.append((
                i,
                start_line,
                end_line,
                self.config.kmer_matrix,
                self.config.use_reverse_complement
            ))

        # 步骤4：多进程执行（使用initializer共享目标kmer数据）|Step 4: Execute with multiprocessing (using initializer to share target kmer data)
        self.logger.info("开始并行处理（使用共享内存）|Starting parallel processing (using shared memory)...")

        all_results = {}
        total_found_forward = 0
        total_found_reverse = 0
        total_scanned = 0

        # 准备共享数据|Prepare shared data
        target_set = set(target_kmers.keys())

        with Pool(
            processes=num_processes,
            initializer=_worker_init,
            initargs=(target_kmers, target_set)
        ) as pool:
            # 使用imap_unordered以便实时获取进度|Use imap_unordered for real-time progress
            for i, (chunk_id, chunk_results, found_count, stats) in enumerate(
                pool.imap_unordered(self._process_chunk_worker, tasks)
            ):
                if 'error' in stats:
                    self.logger.error(f"Chunk {chunk_id} 出错|Error: {stats['error']}")
                    continue

                # 合并结果|Merge results
                all_results.update(chunk_results)
                total_found_forward += stats['found_forward']
                total_found_reverse += stats['found_reverse']
                total_scanned += stats['total_lines']

                # 进度显示|Progress display
                completed = i + 1
                self.logger.info(f"完成chunk|Completed chunk {completed}/{num_processes}, "
                               f"找到|Found: {format_number(found_count)} "
                               f"(总|Total: {format_number(len(all_results))})")

        not_found = len(target_kmers) - len(all_results)

        self.logger.info(f"扫描完成|Scan completed")
        self.logger.info(f"总扫描行数|Total lines scanned: {format_number(total_scanned)}")
        self.logger.info(f"正向找到|Found forward: {format_number(total_found_forward)}")
        self.logger.info(f"反向互补找到|Found reverse complement: {format_number(total_found_reverse)}")
        self.logger.info(f"未找到|Not found: {format_number(not_found)}")
        self.logger.info(f"总共找到|Total found: {format_number(len(all_results))}")

        return all_results, len(all_results), sample_names

    def _count_lines(self, file_path: str) -> int:
        """
        快速统计文件行数|Fast count lines in file

        Args:
            file_path: 文件路径|File path

        Returns:
            int: 行数|Number of lines
        """
        count = 0
        open_func = gzip.open if file_path.endswith('.gz') else open

        with open_func(file_path, 'rb') as f:
            # 跳过header|Skip header
            f.readline()
            count = 0

            # 按块读取以加速|Read in chunks for speed
            while True:
                chunk = f.read(1024 * 1024)  # 1MB chunks
                if not chunk:
                    break
                count += chunk.count(b'\n')

        return count

    def write_results(
        self,
        results: Dict[str, dict],
        sample_names: list,
        target_kmers: Dict[str, str],
        kmer_id_list: list = None,
        kmer_id_to_seq: Dict[str, str] = None
    ):
        """
        写入结果到输出文件|Write results to output file

        修改：支持重复的kmer ID|Modified: Support duplicate kmer IDs

        Args:
            results: 结果字典|Results dictionary
            sample_names: 样本名列表|List of sample names
            target_kmers: 目标kmer字典 (序列 -> ID)|Target kmer dictionary (sequence -> ID)
            kmer_id_list: kmer ID列表（保持fasta原始顺序，包括重复）|Kmer ID list preserving fasta order including duplicates
            kmer_id_to_seq: kmer ID到序列的映射|Kmer ID to sequence mapping
        """
        self.logger.info("=" * 100)
        self.logger.info(" 写入结果文件|Writing Results File")
        self.logger.info("=" * 100)
        self.logger.info(f"输出文件|Output file: {self.config.output_file}")

        written_count = 0

        # 判断是否为压缩文件|Check if output should be compressed
        open_func = gzip.open if self.config.output_file.endswith('.gz') else open

        with open_func(self.config.output_file, 'wt') as f_out:
            # 写入header行|Write header line
            if self.config.output_format == 'csv':
                sep = ','
            else:
                sep = '\t'

            header = ['KmerID', 'Sequence', 'FoundAs'] + sample_names
            f_out.write(sep.join(header) + '\n')

            # 创建sequence -> result的映射，用于处理重复ID|Create sequence -> result mapping for handling duplicate IDs
            seq_to_result = {}
            for kmer_id, result in results.items():
                seq_to_result[result['sequence']] = result

            # 创建反向映射：kmer_id -> sequence，用于快速查找|Create reverse mapping: kmer_id -> sequence for fast lookup
            id_to_seq = {v: k for k, v in target_kmers.items()}

            # 按照fasta文件中的原始顺序写入|Write in the original order from fasta file
            if kmer_id_list:
                # 使用提供的kmer_id_list（保持fasta原始顺序，包括重复）|Use provided kmer_id_list (preserves fasta order including duplicates)
                for kmer_id in kmer_id_list:
                    # 从kmer_id_to_seq获取序列|Get sequence from kmer_id_to_seq
                    kmer_seq = kmer_id_to_seq.get(kmer_id, '')

                    if not kmer_seq:
                        # 回退到旧方法|Fallback to old method
                        kmer_seq = id_to_seq.get(kmer_id, '')

                    if kmer_seq in seq_to_result:
                        # 找到了该序列的结果|Found result for this sequence
                        result = seq_to_result[kmer_seq]
                        row = [kmer_id, result['sequence'], result['found_as']] + result['abundances']
                        f_out.write(sep.join(row) + '\n')
                        written_count += 1
                    elif self.config.keep_not_found:
                        # 未找到但保留|Not found but keep
                        row = [kmer_id, kmer_seq, 'not_found'] + ['NA'] * len(sample_names)
                        f_out.write(sep.join(row) + '\n')
                        written_count += 1
                    # else: 未找到且不保留，跳过|Not found and skip
            else:
                # 回退到按字典序排序（如果没有提供kmer_id_list）|Fallback to lexicographic order (if kmer_id_list not provided)
                for kmer_seq, kmer_id in sorted(target_kmers.items(), key=lambda x: x[1]):
                    if kmer_id in results:
                        # 找到了|Found
                        result = results[kmer_id]
                        row = [kmer_id, result['sequence'], result['found_as']] + result['abundances']
                    elif self.config.keep_not_found:
                        # 未找到但保留|Not found but keep
                        row = [kmer_id, kmer_seq, 'not_found'] + ['NA'] * len(sample_names)
                    else:
                        # 未找到且不保留|Not found and skip
                        continue

                    f_out.write(sep.join(row) + '\n')
                    written_count += 1

        self.logger.info(f"写入完成|Write completed")
        self.logger.info(f"输出行数|Output rows: {format_number(written_count)}")

        # 写入窗口统计文件|Write window statistics file
        if kmer_id_list and self.config.window_size > 0:
            self._write_window_stats(results, kmer_id_list, sample_names)

    def _write_window_stats(
        self,
        results: Dict[str, dict],
        kmer_id_list: list,
        sample_names: list
    ):
        """
        写入窗口统计文件|Write window statistics file

        按照指定的窗口大小统计每个窗口内每个样品的kmer丰度总和
        Calculate the sum of kmer abundances for each sample in each window

        Args:
            results: 结果字典|Results dictionary
            kmer_id_list: kmer ID列表（保持fasta原始顺序）|Kmer ID list preserving fasta order
            sample_names: 样本名列表|List of sample names
        """
        window_stats_file = self.config.output_file + '.window_stats.txt'

        self.logger.info("=" * 100)
        self.logger.info(" 写入窗口统计文件|Writing Window Statistics File")
        self.logger.info("=" * 100)
        self.logger.info(f"窗口统计文件|Window stats file: {window_stats_file}")
        self.logger.info(f"窗口大小|Window size: {format_number(self.config.window_size)}")

        window_size = self.config.window_size
        total_kmers = len(kmer_id_list)

        # 计算窗口数|Calculate number of windows
        num_windows = (total_kmers + window_size - 1) // window_size

        with open(window_stats_file, 'wt') as f_out:
            # 写入header|Write header
            header = ['Window', 'Start', 'End'] + sample_names
            f_out.write('\t'.join(header) + '\n')

            # 统计每个窗口|Statistics for each window
            for window_idx in range(num_windows):
                start_idx = window_idx * window_size
                end_idx = min((window_idx + 1) * window_size, total_kmers)

                # 获取当前窗口的kmer IDs|Get kmer IDs for current window
                window_kmers = kmer_id_list[start_idx:end_idx]

                # 初始化每个样品的丰度总和|Initialize abundance sum for each sample
                sample_sums = [0] * len(sample_names)

                # 累加每个kmer在各个样品中的丰度|Accumulate abundances for each kmer across samples
                for kmer_id in window_kmers:
                    if kmer_id in results:
                        # 该kmer在矩阵中找到，累加各样品丰度
                        abundances = results[kmer_id]['abundances']
                        # 只处理前len(sample_names)个丰度值，避免索引越界
                        # Only process first len(sample_names) abundance values to avoid index out of range
                        for i, abundance in enumerate(abundances[:len(sample_names)]):
                            try:
                                sample_sums[i] += int(abundance)
                            except (ValueError, TypeError):
                                # 如果丰度值不是数字，跳过|Skip if abundance is not a number
                                pass

                        # 如果abundances比sample_names长，记录警告
                        # If abundances is longer than sample_names, log warning
                        if len(abundances) > len(sample_names):
                            self.logger.warning(f"Kmer {kmer_id} 有 {len(abundances)} 个丰度值，"
                                              f"但只有 {len(sample_names)} 个样品，"
                                              f"多余的 {len(abundances) - len(sample_names)} 个值被忽略")
                    # 如果kmer没找到，丰度视为0，不需要操作|If kmer not found, abundance is 0, no action needed

                # 输出（使用1-based索引）|Output (using 1-based indexing)
                window_num = window_idx + 1
                start_num = start_idx + 1
                end_num = end_idx

                row = [str(window_num), str(start_num), str(end_num)] + [str(s) for s in sample_sums]
                f_out.write('\t'.join(row) + '\n')

        self.logger.info(f"窗口统计文件写入完成|Window stats file written")
        self.logger.info(f"总窗口数|Total windows: {num_windows}")
        self.logger.info(f"总kmer数|Total kmers: {format_number(total_kmers)}")
        self.logger.info(f"样品数|Number of samples: {len(sample_names)}")
