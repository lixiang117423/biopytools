"""
基因组Gap统计主程序模块|Genome Gap Statistics Main Module
"""

import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Dict

from .config import GapStatConfig
from .utils import GapStatLogger


class GapStat:
    """基因组Gap统计主类|Genome Gap Statistics Main Class"""

    def __init__(self, **kwargs):
        """初始化|Initialize"""
        # 初始化配置|Initialize configuration
        self.config = GapStatConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        if self.config.output_file:
            output_dir = Path(self.config.output_file).parent
            log_file = output_dir / "gapstat.log"
        else:
            log_file = None
        self.logger_manager = GapStatLogger(log_file)
        self.logger = self.logger_manager.get_logger()

    def find_gaps(self) -> List[Tuple[str, int, int, int]]:
        """
        查找FASTA文件中的所有gap|Find all gaps in FASTA file

        Returns:
            List[Tuple]: [(序列名, 起始位置, 终止位置, gap长度), ...] (1-based坐标)
        """
        self.logger.info(f"解析FASTA文件|Parsing FASTA file: {self.config.fasta_file}")

        results = []
        current_name = None
        current_seq = []
        seq_count = 0

        with open(self.config.fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    # 处理上一条序列|Process previous sequence
                    if current_name:
                        seq = "".join(current_seq)
                        gaps = self._parse_gaps(current_name, seq)
                        results.extend(gaps)
                        seq_count += 1

                    # 提取序列名|Extract sequence name
                    current_name = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

        # 处理最后一条序列|Process last sequence
        if current_name:
            seq = "".join(current_seq)
            gaps = self._parse_gaps(current_name, seq)
            results.extend(gaps)
            seq_count += 1

        self.logger.info(f"解析完成|Parsing completed: {seq_count} 条序列|sequences, {len(results)} 个gap|gaps found")
        return results

    def _parse_gaps(self, name: str, seq: str) -> List[Tuple[str, int, int, int]]:
        """
        解析单条序列中的gap|Parse gaps in a single sequence

        Args:
            name: 序列名|Sequence name
            seq: 序列|Sequence

        Returns:
            List[Tuple]: [(序列名, 起始, 终止, 长度), ...] (1-based坐标)
        """
        gaps = []
        # 正则匹配至少min_n个连续N|Regex match at least min_n consecutive Ns
        pattern = f"N{{{self.config.min_n},}}"

        for match in re.finditer(pattern, seq, re.IGNORECASE):
            start = match.start() + 1  # 转为1-based|Convert to 1-based
            end = match.end()
            length = end - start + 1
            gaps.append((name, start, end, length))

        return gaps

    def calculate_statistics(self, gaps: List[Tuple[str, int, int, int]]) -> Dict:
        """
        计算gap统计信息|Calculate gap statistics

        Args:
            gaps: gap列表|List of gaps

        Returns:
            Dict: 统计信息字典|Statistics dictionary
        """
        self.logger.info("计算统计信息|Calculating statistics")

        stats = {
            'total_gaps': len(gaps),
            'total_length': sum(g[3] for g in gaps),
            'gaps_by_sequence': defaultdict(int),
            'gaps_by_length': defaultdict(int),
            'sequence_gap_lengths': defaultdict(list)
        }

        # 按序列统计|Statistics by sequence
        for seq_name, start, end, length in gaps:
            stats['gaps_by_sequence'][seq_name] += 1
            stats['gaps_by_length'][length] += 1
            stats['sequence_gap_lengths'][seq_name].append(length)

        # 每条序列的gap总长度|Total gap length per sequence
        stats['sequence_total_gap_length'] = {
            seq: sum(lengths) for seq, lengths in stats['sequence_gap_lengths'].items()
        }

        return stats

    def run(self) -> bool:
        """
        运行gap统计流程|Run gap statistics pipeline

        Returns:
            bool: 是否成功|Success status
        """
        try:
            self.logger.info("开始基因组Gap统计|Starting genome gap statistics")

            # 查找gap|Find gaps
            gaps = self.find_gaps()

            if not gaps:
                self.logger.info("未找到gap|No gaps found")
                # 仍然输出表头|Still output header
                self._write_results([])
                return True

            # 计算统计信息|Calculate statistics
            stats = self.calculate_statistics(gaps)

            # 写入结果|Write results
            self._write_results(gaps, stats)

            # 输出统计信息到stderr|Output statistics to stderr
            self._log_statistics(stats)

            self.logger.info("Gap统计完成|Gap statistics completed successfully")
            return True

        except Exception as e:
            self.logger.error(f"流程执行失败|Pipeline execution failed: {str(e)}")
            import traceback
            self.logger.debug(traceback.format_exc())
            return False

    def _write_results(self, gaps: List[Tuple[str, int, int, int]], stats: Dict = None):
        """写入结果到文件或终端|Write results to file or stdout"""
        # 确定输出目标|Determine output target
        if self.config.output_file:
            out_file = open(self.config.output_file, 'w')
            self.logger.info(f"写入结果到文件|Writing results to file: {self.config.output_file}")
        else:
            out_file = sys.stdout
            self.logger.info("写入结果到终端|Writing results to stdout")

        try:
            # 写入表头|Write header
            out_file.write("seq_name\tstart\tend\tgap_length\n")

            # 写入gap数据|Write gap data
            for seq_name, start, end, length in gaps:
                out_file.write(f"{seq_name}\t{start}\t{end}\t{length}\n")

            # 注意：统计信息只输出到stderr，不写入结果文件|Statistics only go to stderr, not to output file

        finally:
            if self.config.output_file:
                out_file.close()

    def _log_statistics(self, stats: Dict):
        """记录统计信息到stderr|Log statistics to stderr"""
        print("\n" + "=" * 60, file=sys.stderr)
        print("汇总|Summary", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        print(f"总gap数量|Total gaps: {stats['total_gaps']}", file=sys.stderr)
        print(f"总gap长度|Total gap length: {stats['total_length']} bp", file=sys.stderr)
        print(f"包含gap的序列数|Sequences with gaps: {len(stats['gaps_by_sequence'])}", file=sys.stderr)


def batch_process_gapstat(fasta_files: List[str], output_file: str, min_n: int = 1) -> bool:
    """
    批量处理多个FASTA文件的gap统计|Batch process gap statistics for multiple FASTA files

    Args:
        fasta_files: FASTA文件列表|List of FASTA files
        output_file: 输出文件路径|Output file path
        min_n: 最小连续N数量|Minimum consecutive N count

    Returns:
        bool: 是否成功|Success status
    """
    # 创建输出目录|Create output directory
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # 初始化日志|Initialize logging
    log_file = output_dir / "gapstat_batch.log"
    logger_manager = GapStatLogger(str(log_file))
    logger = logger_manager.get_logger()

    logger.info(f"开始批量Gap统计|Starting batch gap statistics")
    logger.info(f"文件数量|Number of files: {len(fasta_files)}")
    logger.info(f"输出文件|Output file: {output_file}")

    # 存储所有结果|Store all results
    all_gaps = []
    file_stats = {}

    # 处理每个文件|Process each file
    for fasta_file in fasta_files:
        # 提取样品名称|Extract sample name
        sample_name = Path(fasta_file).stem

        logger.info(f"处理文件|Processing file: {fasta_file} (样品|Sample: {sample_name})")

        try:
            # 创建GapStat实例|Create GapStat instance
            gapstat = GapStat(
                fasta_file=fasta_file,
                output_file=None,  # 不输出单个文件|Don't output single file
                min_n=min_n
            )

            # 查找gap|Find gaps
            gaps = gapstat.find_gaps()

            # 添加样品名称前缀|Add sample name prefix
            for seq_name, start, end, length in gaps:
                all_gaps.append((sample_name, seq_name, start, end, length))

            # 计算统计|Calculate statistics
            if gaps:
                stats = gapstat.calculate_statistics(gaps)
                file_stats[sample_name] = stats

        except Exception as e:
            logger.error(f"处理文件失败|Failed to process file {fasta_file}: {str(e)}")
            return False

    logger.info(f"解析完成|Parsing completed: {len(all_gaps)} 个gap|gaps total")

    # 写入结果|Write results
    try:
        logger.info(f"写入输出文件|Writing output file: {output_file}")

        with open(output_file, 'w', encoding='utf-8') as f:
            # 写入表头|Write header
            f.write("sample\tseq_name\tstart\tend\tgap_length\n")

            # 写入gap数据|Write gap data
            for sample, seq_name, start, end, length in all_gaps:
                f.write(f"{sample}\t{seq_name}\t{start}\t{end}\t{length}\n")

            # 写入统计信息|Write statistics
            f.write("\n" + "=" * 80 + "\n")
            f.write("统计信息|Statistics\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"总记录数|Total records: {len(all_gaps)}\n")
            f.write(f"样品数|Number of samples: {len(file_stats)}\n\n")

            total_gaps = len(all_gaps)
            total_length = sum(g[4] for g in all_gaps)
            f.write(f"总gap数量|Total gaps: {total_gaps}\n")
            f.write(f"总gap长度|Total gap length: {total_length} bp\n\n")

            f.write("各样品统计|Sample statistics:\n")
            for sample in sorted(file_stats.keys()):
                stats = file_stats[sample]
                f.write(f"\n  {sample}:\n")
                f.write(f"    Gap数量|Gap count: {stats['total_gaps']}\n")
                f.write(f"    Gap总长度|Total gap length: {stats['total_length']} bp\n")
                f.write(f"    序列数|Sequences: {len(stats['gaps_by_sequence'])}\n")

        logger.info("批量转换完成|Batch processing completed successfully")

        # 输出汇总到stderr|Output summary to stderr
        print("\n" + "=" * 60, file=sys.stderr)
        print("批量处理完成|Batch processing completed", file=sys.stderr)
        print(f"总gap数量|Total gaps: {total_gaps}", file=sys.stderr)
        print(f"总gap长度|Total gap length: {total_length} bp", file=sys.stderr)
        print("=" * 60, file=sys.stderr)

        return True

    except Exception as e:
        logger.error(f"写入输出文件失败|Failed to write output file: {str(e)}")
        import traceback
        logger.debug(traceback.format_exc())
        return False


def main():
    """主函数|Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='基因组Gap统计工具|Genome Gap Statistics Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('fasta',
                       help='输入FASTA文件|Input FASTA file')
    parser.add_argument('--min-n', type=int, default=1,
                       help='最少N的数量|Minimum consecutive N count')
    parser.add_argument('-o', '--output',
                       help='输出文件|Output file (default: stdout)')

    args = parser.parse_args()

    # 创建GapStat并运行|Create GapStat and run
    gapstat = GapStat(
        fasta_file=args.fasta,
        min_n=args.min_n,
        output_file=args.output
    )

    success = gapstat.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
