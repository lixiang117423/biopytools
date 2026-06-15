"""
基因组装配统计分析核心模块|Genome Assembly Statistics Analysis Core Module
"""

from pathlib import Path
from typing import Dict, List, Tuple
import re

from .utils import parse_fasta, parse_fastq


class SequenceStats:
    """序列统计信息类|Sequence Statistics Information Class"""

    def __init__(self, header: str, sequence: str):
        self.header = header
        self.sequence = sequence
        self.length = len(sequence)

    def count_gaps(self) -> int:
        """计算gap数量|Count gaps"""
        # Gap是任何连续的N或n碱基|Gap is any consecutive run of Ns or ns
        gaps = re.findall(r'[Nn]+', self.sequence)
        return len(gaps)

    def count_n_bases(self) -> int:
        """计算N碱基数量|Count N bases"""
        return self.sequence.upper().count('N')


class AssemblyStatsAnalyzer:
    """基因组装配统计分析器|Genome Assembly Statistics Analyzer"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def analyze_file(self, file_path: str) -> Dict:
        """
        分析单个文件|Analyze single file

        Args:
            file_path: 文件路径|File path

        Returns:
            dict: 统计结果|Statistics results
        """
        self.logger.info(f"分析文件|Analyzing file: {file_path}")

        # 判断文件类型|Determine file type
        file_ext = Path(file_path).suffix.lower()

        # 解析序列|Parse sequences
        sequences = []
        if file_ext in ['.fa', '.fasta', '.fna']:
            parser = parse_fasta(file_path, self.config.min_length)
            for header, seq in parser:
                sequences.append(SequenceStats(header, seq))
        elif file_ext in ['.fq', '.fastq']:
            parser = parse_fastq(file_path, self.config.min_length)
            for header, seq in parser:
                sequences.append(SequenceStats(header, seq))
        else:
            # 尝试FASTA格式|Try FASTA format
            try:
                parser = parse_fasta(file_path, self.config.min_length)
                for header, seq in parser:
                    sequences.append(SequenceStats(header, seq))
            except:
                self.logger.error(f"无法解析文件格式|Cannot parse file format: {file_path}")
                return None

        if not sequences:
            self.logger.warning(f"未找到有效序列|No valid sequences found in: {file_path}")
            return None

        # 计算统计信息|Calculate statistics
        stats = self._calculate_statistics(sequences, file_path)
        return stats

    def _calculate_statistics(self, sequences: List[SequenceStats], file_path: str) -> Dict:
        """
        计算统计信息|Calculate statistics

        Args:
            sequences: 序列列表|List of sequences
            file_path: 文件路径|File path

        Returns:
            dict: 统计结果|Statistics results
        """
        # 按长度排序（降序）|Sort by length (descending)
        sequences.sort(key=lambda x: x.length, reverse=True)

        # 基本统计|Basic statistics
        n_seqs = len(sequences)
        total_length = sum(s.length for s in sequences)
        avg_length = total_length / n_seqs if n_seqs > 0 else 0
        largest = sequences[0].length if sequences else 0

        # 计算Nx|Calculate Nx
        nx_stats = {}
        for percentage in [50, 60, 70, 80, 90, 100]:
            nx_value, n_seqs_needed = self._calculate_nx(sequences, total_length, percentage)
            nx_stats[f'N{percentage}'] = nx_value
            nx_stats[f'N{percentage}_n'] = n_seqs_needed

        # Gap统计|Gap statistics
        total_gaps = sum(s.count_gaps() for s in sequences)
        total_n_count = sum(s.count_n_bases() for s in sequences)

        # 组装结果|Assemble results
        stats = {
            'file_name': Path(file_path).name,
            'sum': total_length,
            'n': n_seqs,
            'ave': avg_length,
            'largest': largest,
            'gaps': total_gaps,
            'n_count': total_n_count,
            **nx_stats
        }

        return stats

    def _calculate_nx(self, sequences: List[SequenceStats], total_length: int, percentage: int) -> Tuple[int, int]:
        """
        计算Nx值|Calculate Nx value

        Args:
            sequences: 序列列表（已排序）|List of sequences (sorted)
            total_length: 总长度|Total length
            percentage: 百分比|Percentage (50, 60, 70, 80, 90, 100)

        Returns:
            tuple: (Nx值, 序列数)|(Nx value, number of sequences)
        """
        target_length = total_length * percentage / 100.0
        cumulative_length = 0

        for i, seq in enumerate(sequences):
            cumulative_length += seq.length
            if cumulative_length >= target_length:
                return (seq.length, i + 1)

        return (0, len(sequences))

    def format_output(self, stats: Dict) -> str:
        """
        格式化输出结果|Format output results

        Args:
            stats: 统计结果|Statistics results

        Returns:
            str: 格式化输出|Formatted output
        """
        if stats is None:
            return ""

        file_name = stats['file_name']

        # 辅助函数：格式化大数字|Helper function: format large numbers
        def format_number(num: int) -> str:
            """格式化数字，以M为单位|Format number in M unit"""
            if num >= 1_000_000:
                return f"{num / 1_000_000:.2f}M"
            return str(num)

        # 竖着格式（键值对格式）|Vertical format (key-value pair format)
        if self.config.vertical_format:
            output = []
            output.append(f"File|文件\t{file_name}")
            output.append(f"Sum|总和\t{format_number(stats['sum'])}")
            output.append(f"N|序列数\t{stats['n']}")
            output.append(f"Average|平均长度\t{format_number(int(stats['ave']))}")
            output.append(f"Largest|最大长度\t{format_number(stats['largest'])}")

            # Nx统计|Nx statistics
            for percentage in [50, 60, 70, 80, 90, 100]:
                nx_key = f'N{percentage}'
                nx_n_key = f'N{percentage}_n'
                output.append(f"{nx_key}\t{format_number(stats[nx_key])}")
                output.append(f"{nx_key}_count|{nx_key}_序列数\t{stats[nx_n_key]}")

            # Gap和N统计|Gap and N statistics
            output.append(f"N_count|N碱基数\t{format_number(stats['n_count'])}")
            output.append(f"Gaps|Gap数\t{stats['gaps']}")

            return "\n".join(output)

        # 原有横着格式|Original horizontal format
        output = []

        # 文件名|File name
        if not self.config.grep_friendly:
            output.append(f"stats for {file_name}")

        # 基本统计|Basic statistics
        output.append(f"sum = {format_number(stats['sum'])}, n = {stats['n']}, ave = {format_number(int(stats['ave']))}, largest = {format_number(stats['largest'])}")

        # Nx统计|Nx statistics
        for percentage in [50, 60, 70, 80, 90, 100]:
            nx_key = f'N{percentage}'
            nx_n_key = f'N{percentage}_n'
            output.append(f"{nx_key} = {format_number(stats[nx_key])}, n = {stats[nx_n_key]}")

        # Gap和N统计|Gap and N statistics
        output.append(f"N_count = {format_number(stats['n_count'])}")
        output.append(f"Gaps = {stats['gaps']}")

        if self.config.tab_delimited:
            # Tab分隔输出|Tab-delimited output
            if not self.config.no_header:
                header = ["File", "Sum", "N", "Ave", "Largest"] + \
                        [f"N{x}" for x in [50, 60, 70, 80, 90, 100]] + \
                        ["N_count", "Gaps"]
                if not self.config.grep_friendly:
                    output = ["\t".join(header)]

            values = [
                file_name,
                format_number(stats['sum']),
                str(stats['n']),
                format_number(int(stats['ave'])),
                format_number(stats['largest'])
            ] + [format_number(stats[f'N{x}']) for x in [50, 60, 70, 80, 90, 100]] + \
            [format_number(stats['n_count']), str(stats['gaps'])]

            output.append("\t".join(values))
            return "\n".join(output)

        return "\n".join(output)
