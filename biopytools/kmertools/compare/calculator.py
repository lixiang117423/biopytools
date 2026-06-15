"""
核心计算逻辑模块|Core Calculation Logic Module
"""

import gzip
from typing import Dict, Set, List, Tuple
from .utils import format_number


class KmerCompareCalculator:
    """Kmer矩阵比较计算器|Kmer Matrix Comparison Calculator"""

    def __init__(self, config, logger):
        """
        初始化计算器|Initialize calculator

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger

    def extract_sequences(self, file_path: str) -> Set[str]:
        """
        提取文件中所有的kmer序列（流式处理）|Extract all kmer sequences from file (streaming)

        Args:
            file_path: 文件路径|File path

        Returns:
            Set[str]: kmer序列集合|Set of kmer sequences
        """
        self.logger.info(f"提取kmer序列|Extracting kmer sequences: {file_path}")

        sequences = set()
        open_func = gzip.open if file_path.endswith('.gz') else open

        with open_func(file_path, 'rt') as f:
            # 跳过header|Skip header
            next(f)

            line_count = 0
            for line in f:
                line_count += 1
                if line_count % 1000000 == 0:
                    self.logger.info(f"  已读取|Read: {format_number(line_count)} 行|lines")

                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    sequence = parts[1]
                    sequences.add(sequence)

        self.logger.info(f"  提取完成|Extracted: {format_number(len(sequences))} 序列|sequences")
        return sequences

    def find_unique_sequences(self, seqs1: Set[str], seqs2: Set[str]) -> Tuple[Set[str], Set[str]]:
        """
        找出两个集合中特有的序列|Find unique sequences in two sets

        Args:
            seqs1: 第一个文件的序列集合|Sequence set from first file
            seqs2: 第二个文件的序列集合|Sequence set from second file

        Returns:
            Tuple[Set[str], Set[str]]: (文件1特有序列, 文件2特有序列)
        """
        self.logger.info("查找特有序列|Finding unique sequences...")

        unique1 = seqs1 - seqs2
        unique2 = seqs2 - seqs1
        common = seqs1 & seqs2

        self.logger.info(f"文件1特有序列|Unique sequences in file1: {format_number(len(unique1))}")
        self.logger.info(f"文件2特有序列|Unique sequences in file2: {format_number(len(unique2))}")
        self.logger.info(f"共同序列|Common sequences: {format_number(len(common))}")

        return unique1, unique2

    def process_file_with_stats(self, file_path: str, unique_seqs: Set[str],
                                output_file: str, file_id: str) -> List[str]:
        """
        流式处理文件，标注unique并计算窗口统计|Stream process file, label unique and calculate window stats

        Args:
            file_path: 输入文件路径|Input file path
            unique_seqs: 特有序列集合|Set of unique sequences
            output_file: 输出文件路径|Output file path
            file_id: 文件标识|File identifier

        Returns:
            List[str]: 样品名列表|List of sample names
        """
        self.logger.info(f"处理文件|Processing file: {file_path}")
        self.logger.info("标注unique kmer并计算窗口统计|Labeling unique kmers and calculating window stats...")

        open_func = gzip.open if file_path.endswith('.gz') else open

        # 存储窗口统计信息|Store window statistics
        windows = []
        window_valid_counts = []

        sample_names = []
        line_count = 0
        window_idx = 0

        with open_func(file_path, 'rt') as f:
            # 读取header|Read header
            header_line = next(f).strip()
            header_parts = header_line.split('\t')

            if len(header_parts) >= 4:
                sample_names = header_parts[3:]

            # 初始化第一个窗口|Initialize first window
            windows.append([0] * len(sample_names))
            window_valid_counts.append(0)

            # 逐行处理|Process line by line
            for line in f:
                line_count += 1
                if line_count % 1000000 == 0:
                    self.logger.info(f"  已处理|Processed: {format_number(line_count)} 行|lines")

                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue

                sequence = parts[1]
                found_as = parts[2]
                abundances = parts[3:]

                # 检查是否为unique|Check if unique
                is_unique = sequence in unique_seqs

                # 更新当前窗口的统计|Update current window statistics
                if is_unique:
                    # 所有unique kmer都计入分母（包括not_found）
                    window_valid_counts[window_idx] += 1

                    # 只统计非not_found的kmer在各样品中的0/1值
                    # only_found kmer中的0/1值
                    if found_as != 'not_found':
                        # 对每个样品单独统计|Statistics for each sample separately
                        for sample_idx, abundance in enumerate(abundances):
                            if sample_idx >= len(sample_names):
                                break
                            try:
                                if int(abundance) == 1:
                                    windows[window_idx][sample_idx] += 1
                            except (ValueError, TypeError):
                                pass

                # 检查是否需要切换到下一个窗口|Check if need to move to next window
                if line_count % self.config.window_size == 0:
                    window_idx += 1
                    windows.append([0] * len(sample_names))
                    window_valid_counts.append(0)

        self.logger.info(f"  处理完成|Processed: {format_number(line_count)} 行|lines")
        self.logger.info(f"  窗口数|Number of windows: {len(windows)}")

        # 写入统计结果|Write statistics results
        self.write_window_stats(windows, window_valid_counts, sample_names,
                               output_file, file_id, line_count)

        return sample_names

    def write_window_stats(self, windows: List[List[int]], window_valid_counts: List[int],
                          sample_names: List[str], output_file: str, file_id: str, total_lines: int):
        """
        写入窗口统计结果|Write window statistics results

        Args:
            windows: 窗口统计信息|Window statistics
            window_valid_counts: 每个窗口的有效kmer数|Valid kmer count per window
            sample_names: 样品名列表|List of sample names
            output_file: 输出文件路径|Output file path
            file_id: 文件标识|File identifier
            total_lines: 总行数|Total lines
        """
        self.logger.info(f"写入统计结果|Writing statistics results: {output_file}")

        with open(output_file, 'w') as f_out:
            # 写入header|Write header
            header = ['File_ID', 'Start', 'End'] + [f'{name}_ratio' for name in sample_names]
            f_out.write('\t'.join(header) + '\n')

            # 写入每个窗口的统计|Write statistics for each window
            num_windows = len(windows)
            for window_idx in range(num_windows):
                start_line = window_idx * self.config.window_size + 1
                end_line = min((window_idx + 1) * self.config.window_size, total_lines)

                # 计算比例|Calculate ratios
                valid_count = window_valid_counts[window_idx]

                if valid_count > 0:
                    ratios = [f"{windows[window_idx][i] / valid_count:.4f}"
                             for i in range(len(sample_names))]
                else:
                    ratios = ['0.0000'] * len(sample_names)

                row = [file_id, str(start_line), str(end_line)] + ratios
                f_out.write('\t'.join(row) + '\n')

        self.logger.info(f"  写入完成|Write completed: {num_windows} 窗口|windows")
